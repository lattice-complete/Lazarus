//! Core reduced statement, per-round proof message, and the final direct
//! verification (C: statement / proof / verify()).

use crate::comkey::{ComKey, KeyDomain};
use crate::params::{sis_secure, ComParams, OPNORM_T, SLACK};
use crate::transcript::Hash;
use crate::witness::Witness;
use crate::Rq;
use algebra::ring::PolyVec;

/// Upper-triangular (row <= col) pair index, row-major (C: triangularidx).
pub fn pair_idx(i: usize, j: usize, r: usize) -> usize {
    debug_assert!(i <= j && j < r);
    i * r - i * (i + 1) / 2 + i + (j - i)
}

pub fn npairs(r: usize) -> usize {
    (r * r + r) / 2
}

/// The reduced statement after one round (C: statement).
///
/// Relations claimed about a witness (z-digits s_0..s_{f-1}, aux v):
///   B v_t = u1 (with C v_g folded in),  D v_h = u2,
///   A z = sum_i c_i t_i,   <z,z> = sum c_i c_j g_ij (quadratic only),
///   <phi, z> = sum c_i c_j h_ij,   sum a_ij g_ij + sum h_ii = b,
///   sum ||parts||^2 <= betasq,
/// where z = sum_i 2^(i b) s_i and t/g/h are reconstructed from v.
#[derive(Debug, Clone)]
pub struct Statement {
    pub r: usize,
    pub n: usize,
    pub m: usize,
    pub tail: bool,
    pub quadratic: bool,
    pub cpp: ComParams,
    /// Outer commitment 1 (kappa1) — or raw inner commitments (r*kappa)
    /// followed by raw quadratic garbage (npairs) in tail mode.
    pub u1: PolyVec<64>,
    /// Outer commitment 2 (kappa1) — or raw linear-garbage values (2r-1)
    /// in tail mode.
    pub u2: PolyVec<64>,
    /// Amortization challenges (r).
    pub c: Vec<Rq>,
    /// Sparse quadratic form over the r chunks (row <= col).
    pub a: Vec<(usize, usize, Rq)>,
    /// Folded linear form (rank n).
    pub phi: PolyVec<64>,
    pub b: Rq,
    pub betasq: u128,
    pub h: Hash,
}

impl Statement {
    /// Offsets of the t / g / h blocks inside the aux vector v.
    pub fn aux_offsets(&self) -> (usize, usize, usize) {
        let l = self.cpp.fu * self.cpp.kappa;
        let t = 0;
        let g = t + self.r * l;
        let h = g + self.cpp.fg * npairs(self.r);
        (t, g, h)
    }
}

/// One round's proof message (C: proof).
#[derive(Debug, Clone)]
pub struct RoundProof {
    /// Input witness ranks and repartition (the verifier revalidates both).
    pub n_in: Vec<usize>,
    pub nu: Vec<usize>,
    pub cpp: ComParams,
    pub tail: bool,
    pub quadratic: bool,
    pub u1: PolyVec<64>,
    pub jlnonce: u64,
    pub p: Vec<i64>, // 256
    /// Lifted Z_q-constraint polynomials, constant term zeroed (LIFTS).
    pub bb: Vec<Rq>,
    pub u2: PolyVec<64>,
    /// Claimed output witness squared norm.
    pub normsq: u128,
}

/// Direct verification of a statement against an in-clear witness
/// (C: verify()) — the end of the recursion chain.
pub fn verify(key: &ComKey, st: &Statement, wt: &Witness) -> Result<(), String> {
    let cpp = &st.cpp;
    let (r, n, m) = (st.r, st.n, st.m);
    let l = cpp.fu * cpp.kappa;
    let (toff, goff, hoff) = st.aux_offsets();
    let pairs = npairs(r);

    /* norm bound */
    if wt.total_normsq() > st.betasq {
        return Err("total witness norm bigger than bound".into());
    }
    let expect_parts = cpp.f + usize::from(!st.tail);
    if wt.parts.len() != expect_parts {
        return Err("wrong witness arity".into());
    }
    for i in 0..cpp.f {
        if wt.parts[i].len() != n {
            return Err("wrong z-digit rank".into());
        }
    }

    /* materialize t / g / h */
    let (tvals, gvals, hvals): (Vec<Rq>, Vec<Rq>, Vec<Rq>) = if !st.tail {
        let v = &wt.parts[cpp.f];
        if v.len() != m {
            return Err("wrong aux rank".into());
        }
        /* outer commitment 1: sum_i B_i t-hat_i + C g-hat = u1 */
        let mut acc = PolyVec::zero(cpp.kappa1);
        for i in 0..r {
            let chunk = &v.polys[toff + i * l..toff + (i + 1) * l];
            acc = &acc + &key.mul(KeyDomain::B, cpp.kappa1, i * l, chunk);
        }
        if cpp.fg > 0 {
            let ghat = &v.polys[goff..goff + cpp.fg * pairs];
            acc = &acc + &key.mul(KeyDomain::C, cpp.kappa1, 0, ghat);
        }
        if acc != st.u1 {
            return Err("first outer commitment opening wrong".into());
        }

        /* outer commitment 2: D h-hat = u2 */
        let hhat = &v.polys[hoff..hoff + cpp.fu * pairs];
        if key.mul(KeyDomain::D, cpp.kappa1, 0, hhat) != st.u2 {
            return Err("second outer commitment opening wrong".into());
        }

        /* reconstruct */
        let recon = |base: usize, len: usize, parts: usize, logb: u32| -> Vec<Rq> {
            let mut out = v.polys[base..base + len].to_vec();
            for d in 1..parts {
                let scale = algebra::zq::Zq::new(1usize << (d as u32 * logb));
                for (k, ok) in out.iter_mut().enumerate() {
                    *ok = &*ok + &(&v.polys[base + d * len + k] * scale);
                }
            }
            out
        };
        // t: per chunk block, digit-major inside
        let mut tv = Vec::with_capacity(r * cpp.kappa);
        for i in 0..r {
            let base = toff + i * l;
            let mut tk = v.polys[base..base + cpp.kappa].to_vec();
            for d in 1..cpp.fu {
                let scale = algebra::zq::Zq::new(1usize << (d as u32 * cpp.bu));
                for (k, tkk) in tk.iter_mut().enumerate() {
                    *tkk = &*tkk + &(&v.polys[base + d * cpp.kappa + k] * scale);
                }
            }
            tv.extend(tk);
        }
        let gv = if cpp.fg > 0 {
            recon(goff, pairs, cpp.fg, cpp.bg)
        } else {
            Vec::new()
        };
        let hv = recon(hoff, pairs, cpp.fu, cpp.bu);
        (tv, gv, hv)
    } else {
        // tail: u1 = raw t (r*kappa) [ || raw g (pairs) ], u2 = raw h (2r-1)
        let tv = st.u1.polys[..r * cpp.kappa].to_vec();
        let gv = if st.quadratic {
            st.u1.polys[r * cpp.kappa..r * cpp.kappa + pairs].to_vec()
        } else {
            Vec::new()
        };
        let hv = st.u2.polys.clone();
        (tv, gv, hv)
    };

    /* reconstruct z */
    let mut z = wt.parts[0].clone();
    for d in 1..cpp.f {
        let scale = algebra::zq::Zq::new(1usize << (d as u32 * cpp.b));
        z = &z + &(&wt.parts[d] * scale);
    }

    /* A z = sum_i c_i t_i */
    let az = key.mul(KeyDomain::A, cpp.kappa, 0, &z.polys);
    let mut ct_sum = PolyVec::zero(cpp.kappa);
    for i in 0..r {
        let ti = PolyVec::from_polys(tvals[i * cpp.kappa..(i + 1) * cpp.kappa].to_vec());
        ct_sum = &ct_sum + &(&ti * &st.c[i]);
    }
    if az != ct_sum {
        return Err("amortized (inner commitment) opening wrong".into());
    }

    /* sum a_ij g_ij + sum h_ii = b */
    let mut acc = -(&st.b);
    for i in 0..r {
        let hii = if st.tail {
            if i == 0 {
                &hvals[0]
            } else {
                &hvals[2 * i]
            }
        } else {
            &hvals[pair_idx(i, i, r)]
        };
        acc = &acc + hii;
    }
    for (i, j, aij) in &st.a {
        acc = &acc + &(aij * &gvals[pair_idx(*i, *j, r)]);
    }
    if !acc.is_zero() {
        return Err("aggregated dot-product constraint doesn't hold".into());
    }

    /* <z,z> = sum_ij c_i c_j g_ij */
    if st.quadratic {
        let mut gfold = Rq::zero();
        for i in 0..r {
            for j in i..r {
                let w = &(&st.c[i] * &st.c[j]) * &gvals[pair_idx(i, j, r)];
                gfold = &gfold + &w;
                if i != j {
                    gfold = &gfold + &w;
                }
            }
        }
        let zz = z.inner_product(&z);
        if zz != gfold {
            return Err("quadratic garbage polynomials wrong".into());
        }
    }

    /* <phi, z> = folded h */
    let hfold = if !st.tail {
        let mut acc = Rq::zero();
        for i in 0..r {
            for j in i..r {
                acc = &acc + &(&(&st.c[i] * &st.c[j]) * &hvals[pair_idx(i, j, r)]);
            }
        }
        acc
    } else {
        // sequential tail folding: c0*(c0 h0) + sum_i c_i*(h_{2i-1} + c_i h_{2i})
        let mut acc = &st.c[0] * &(&st.c[0] * &hvals[0]);
        for i in 1..r {
            let inner = &hvals[2 * i - 1] + &(&st.c[i] * &hvals[2 * i]);
            acc = &acc + &(&st.c[i] * &inner);
        }
        acc
    };
    let phiz = st.phi.inner_product(&z);
    if phiz != hfold {
        return Err("linear garbage polynomials wrong".into());
    }

    /* commitment security under the claimed norm */
    let norm = (st.betasq as f64).sqrt();
    if !sis_secure(
        cpp.kappa,
        6.0 * OPNORM_T * SLACK * 2f64.powi(((cpp.f - 1) as u32 * cpp.b) as i32) * norm,
    ) {
        return Err("inner commitments not secure".into());
    }
    if !st.tail && !sis_secure(cpp.kappa1, 2.0 * SLACK * norm) {
        return Err("outer commitments not secure".into());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangular_index() {
        // r = 3: (0,0)=0 (0,1)=1 (0,2)=2 (1,1)=3 (1,2)=4 (2,2)=5
        assert_eq!(pair_idx(0, 0, 3), 0);
        assert_eq!(pair_idx(0, 2, 3), 2);
        assert_eq!(pair_idx(1, 1, 3), 3);
        assert_eq!(pair_idx(2, 2, 3), 5);
        assert_eq!(npairs(3), 6);
    }
}
