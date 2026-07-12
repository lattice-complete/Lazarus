//! # greyhound — lattice polynomial commitment scheme
//!
//! Rust port of the Greyhound PCS (Nguyen–Seiler, eprint 2024/1293) from
//! the reference implementation in lattice-dogs/labrador (`greyhound.c`),
//! built on the Lazarus `algebra` layer and compressed with the `labrador`
//! recursive proof system.
//!
//! A polynomial with N*len coefficients (packed as `len` ring elements) is
//! arranged into n chunks of m ring elements. Commit: per-chunk base-2^b
//! decomposition s-hat_i, inner Ajtai commitments t_i = A s-hat_i,
//! decomposed and committed in one outer commitment u1 = sum_i B_i t-hat_i.
//! Eval at x with claim y = f(x): the prover sends u2 = W w-hat (a
//! commitment to the partial evaluations w_i = s-hat_i evaluated on the
//! a-side), receives challenges c, and answers z = sum c_i s-hat_i. The
//! whole verification becomes 2*kappa1 + kappa + 2 sparse dot-product
//! constraints over (z-hat0, z-hat1, t-hat, w-hat) — a labrador principal
//! statement (C: polcom_reduce), proven recursively.

use algebra::ring::PolyVec;
use algebra::zq::Zq;
use labrador::comkey::{ComKey, KeyDomain};
use labrador::params::{sis_secure, LOGQ, OPNORM_T, SLACK};
use labrador::principal::{ct_poly, PrincipalStatement, SparseConstraint, SparseEntry};
use labrador::transcript::{self, tags};
use labrador::witness::Witness;
use labrador::Rq;

const DEGREE: usize = 64;
/// Challenge second moment (TAU1 + 4 TAU2).
const CHAL_VAR: f64 = algebra::sampling::challenge::CHALLENGE_NORMSQ as f64;

/// Commitment parameters (C: the comparams subset greyhound uses).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GreyhoundParams {
    pub len: usize,
    /// chunk length (ring elements) and chunk count: len <= m * n.
    pub m: usize,
    pub n: usize,
    /// witness decomposition: f parts, base 2^b.
    pub f: usize,
    pub b: u32,
    /// uniform decomposition: fu parts, base 2^bu.
    pub fu: usize,
    pub bu: u32,
    pub kappa: usize,
    pub kappa1: usize,
}

impl GreyhoundParams {
    pub fn that_len(&self) -> usize {
        self.n * self.kappa * self.fu
    }
    pub fn what_len(&self) -> usize {
        self.fu * self.n
    }
    pub fn zhat_len(&self) -> usize {
        self.m * self.f
    }
}

/// Parameter search (C: init_polcomctx).
pub fn select_params(len: usize) -> Result<GreyhoundParams, String> {
    assert!(len >= 1);
    let nd = DEGREE as f64;
    for f in 2..=8usize {
        let b = (LOGQ as usize + f / 2) / f;
        let mut chosen: Option<(usize, usize, usize, u32, usize, f64)> = None;
        for kappa in 1..=32usize {
            let m = ((len as f64 * (kappa + 1) as f64 / 2.0).sqrt())
                .round()
                .max(1.0) as usize;
            let n = len.div_ceil(m);
            let varz = 2f64.powi(2 * b as i32) / 12.0 * n as f64 * CHAL_VAR;
            let bu = (0.25 * (12.0 * varz).log2()).round().max(1.0) as u32;
            let fu = ((LOGQ as f64 / bu as f64).round() as usize).max(1);

            let mut normsq = (2f64.powi(2 * bu as i32) / 12.0 + varz / 2f64.powi(2 * bu as i32))
                * (m * f) as f64;
            normsq += (2f64.powi(2 * bu as i32) * (fu - 1) as f64
                + 2f64.powi(2 * (LOGQ as i32 - (fu as i32 - 1) * bu as i32)))
                / 12.0
                * ((kappa + 1) * n) as f64;
            normsq *= nd;

            if sis_secure(
                kappa,
                6.0 * OPNORM_T * SLACK * 2f64.powi(bu as i32) * normsq.sqrt(),
            ) {
                chosen = Some((kappa, m, n, bu, fu, normsq));
                break;
            }
        }
        let Some((kappa, m, n, bu, fu, normsq)) = chosen else {
            continue;
        };
        let mut kappa1 = 33usize;
        for k1 in 1..=32usize {
            if sis_secure(k1, 2.0 * SLACK * normsq.sqrt()) {
                kappa1 = k1;
                break;
            }
        }
        if kappa1 <= 32 {
            return Ok(GreyhoundParams {
                len,
                m,
                n,
                f,
                b: b as u32,
                fu,
                bu,
                kappa,
                kappa1,
            });
        }
    }
    Err("cannot make greyhound commitments secure".into())
}

/// Prover-side commitment context.
pub struct CommitCtx {
    pub pp: GreyhoundParams,
    /// Decomposed witness, chunk-major: sx[i] has m*f polys, digit-major
    /// (digit k of chunk element j at index k*m + j).
    pub sx: Vec<Vec<Rq>>,
    /// Decomposed inner commitments, chunk-major blocks of kappa*fu
    /// (digit-major within a block).
    pub that: Vec<Rq>,
    pub u1: PolyVec<64>,
    pub h: transcript::Hash,
}

/// Evaluation proof header (C: polcomprf). The heavy part is the labrador
/// witness / composite proof carried alongside.
#[derive(Debug, Clone)]
pub struct EvalProof {
    pub pp: GreyhoundParams,
    pub x: Zq,
    pub y: Zq,
    pub u1: PolyVec<64>,
    pub u2: PolyVec<64>,
    pub normsq: u128,
}

/// Commit to a polynomial given as `len` ring elements (N*len coefficients,
/// coefficient of X^(N*i + k) = coeffs[i].coeffs[k]).
pub fn commit(key: &ComKey, poly_vec: &[Rq]) -> Result<CommitCtx, String> {
    let pp = select_params(poly_vec.len())?;
    let (m, n) = (pp.m, pp.n);
    let l = pp.kappa * pp.fu;

    let mut sx = Vec::with_capacity(n);
    let mut that = vec![Rq::zero(); n * l];
    let mut u1 = PolyVec::zero(pp.kappa1);

    for i in 0..n {
        // chunk with zero padding
        let mut chunk: Vec<Rq> = Vec::with_capacity(m);
        for j in 0..m {
            chunk.push(poly_vec.get(i * m + j).cloned().unwrap_or_else(Rq::zero));
        }
        // digit-major decomposition: sx_i[k*m + j] = digit k of chunk[j]
        let digits = PolyVec::from_polys(chunk).decompose_topmost(1u64 << pp.b, pp.f);
        let mut sxi = vec![Rq::zero(); m * pp.f];
        for (k, dv) in digits.iter().enumerate() {
            for (j, p) in dv.polys.iter().enumerate() {
                sxi[k * m + j] = p.clone();
            }
        }
        // inner commitment + decomposition
        let ti = key.mul(KeyDomain::A, pp.kappa, 0, &sxi);
        let tdig = ti.decompose_topmost(1u64 << pp.bu, pp.fu);
        for (d, dv) in tdig.iter().enumerate() {
            for (k, p) in dv.polys.iter().enumerate() {
                that[i * l + d * pp.kappa + k] = p.clone();
            }
        }
        // outer commitment
        u1 = &u1 + &key.mul(KeyDomain::B, pp.kappa1, i * l, &that[i * l..(i + 1) * l]);
        sx.push(sxi);
    }

    let mut h = [0u8; 32];
    transcript::absorb_polyvec(&mut h, tags::U1, &u1);
    Ok(CommitCtx {
        pp,
        sx,
        that,
        u1,
        h,
    })
}

/// Evaluate the committed polynomial at scalar x (over Z_q; C: polzvec_eval).
pub fn eval_poly(poly_vec: &[Rq], x: Zq) -> Zq {
    let mut y = Zq::new(0);
    for p in poly_vec.iter().rev() {
        for k in (0..DEGREE).rev() {
            y = y * x + p.coeffs[k];
        }
    }
    y
}

fn zq_bytes(v: Zq) -> [u8; 8] {
    (v.value as u64).to_le_bytes()
}

/// Produce the evaluation proof and the labrador witness for it
/// (C: polcom_eval).
pub fn eval(key: &ComKey, ctx: &CommitCtx, x: Zq, y: Zq) -> (EvalProof, Witness) {
    let pp = &ctx.pp;
    let (m, n) = (pp.m, pp.n);

    // xx[j] = x^(N j)
    let xn = x.pow(DEGREE);
    let mut xx = Vec::with_capacity(m);
    let mut cur = Zq::new(1);
    for _ in 0..m {
        xx.push(cur);
        cur = cur * xn;
    }

    // w_i = sum_k 2^(k b) sum_j sx_i[k m + j] * x^(N j)
    let mut w = Vec::with_capacity(n);
    for sxi in ctx.sx.iter() {
        let mut acc = Rq::zero();
        for k in 0..pp.f {
            let mut part = Rq::zero();
            for (j, &xj) in xx.iter().enumerate() {
                part = &part + &(&sxi[k * m + j] * xj);
            }
            let scale = Zq::new(2).pow(k * pp.b as usize);
            acc = &acc + &(&part * scale);
        }
        w.push(acc);
    }
    // w-hat digit-major with stride n: what[d n + i] = digit d of w_i
    let wdig = PolyVec::from_polys(w).decompose_topmost(1u64 << pp.bu, pp.fu);
    let mut what = vec![Rq::zero(); pp.fu * n];
    for (d, dv) in wdig.iter().enumerate() {
        for (i, p) in dv.polys.iter().enumerate() {
            what[d * n + i] = p.clone();
        }
    }
    let u2 = key.mul(KeyDomain::W, pp.kappa1, 0, &what);

    // transcript: h_commit || x || y || u2 -> challenges
    let mut h = ctx.h;
    transcript::absorb(&mut h, tags::GREYHOUND_XY, &[&zq_bytes(x), &zq_bytes(y)]);
    transcript::absorb_polyvec(&mut h, tags::U2, &u2);
    let c = transcript::challenges(&h, 0, n);

    // z = sum_i c_i sx_i, decomposed into two base-2^bu digits
    let mut z = vec![Rq::zero(); pp.zhat_len()];
    for (i, sxi) in ctx.sx.iter().enumerate() {
        for (zl, s) in z.iter_mut().zip(sxi.iter()) {
            *zl = &*zl + &(s * &c[i]);
        }
    }
    let zdig = PolyVec::from_polys(z).decompose_topmost(1u64 << pp.bu, 2);

    let wt = Witness::new(vec![
        zdig[0].clone(),
        zdig[1].clone(),
        PolyVec::from_polys(ctx.that.clone()),
        PolyVec::from_polys(what),
    ]);
    let normsq = wt.total_normsq();

    (
        EvalProof {
            pp: *pp,
            x,
            y,
            u1: ctx.u1.clone(),
            u2,
            normsq,
        },
        wt,
    )
}

/// Rebuild the labrador principal statement a valid evaluation proof must
/// satisfy (C: polcom_reduce). Runs on the verifier side.
pub fn reduce(key: &ComKey, pi: &EvalProof) -> Result<PrincipalStatement, String> {
    let pp = &pi.pp;
    let (m, n) = (pp.m, pp.n);
    let l = pp.kappa * pp.fu;
    if pi.u1.len() != pp.kappa1 || pi.u2.len() != pp.kappa1 {
        return Err("wrong commitment rank".into());
    }

    /* security under the claimed norm */
    let norm = (pi.normsq as f64).sqrt();
    if !sis_secure(
        pp.kappa,
        6.0 * OPNORM_T * SLACK * 2f64.powi(pp.bu as i32) * norm,
    ) {
        return Err("inner commitments not secure".into());
    }
    if !sis_secure(pp.kappa1, 2.0 * SLACK * norm) {
        return Err("outer commitments not secure".into());
    }

    /* transcript replay */
    let mut h = [0u8; 32];
    transcript::absorb_polyvec(&mut h, tags::U1, &pi.u1);
    transcript::absorb(
        &mut h,
        tags::GREYHOUND_XY,
        &[&zq_bytes(pi.x), &zq_bytes(pi.y)],
    );
    transcript::absorb_polyvec(&mut h, tags::U2, &pi.u2);
    let c = transcript::challenges(&h, 0, n);

    let mut cnst: Vec<SparseConstraint> = Vec::new();

    /* outer commitment 1 openings: kappa1 ring constraints */
    for k in 0..pp.kappa1 {
        let entries = (0..n)
            .map(|i| SparseEntry {
                part: 2,
                off: i * l,
                phi: (0..l)
                    .map(|c_| key.entry(KeyDomain::B, k, i * l + c_))
                    .collect(),
            })
            .collect();
        cnst.push(SparseConstraint {
            ct_only: false,
            entries,
            a: vec![],
            b: pi.u1.polys[k].clone(),
        });
    }

    /* outer commitment 2 openings */
    for k in 0..pp.kappa1 {
        cnst.push(SparseConstraint {
            ct_only: false,
            entries: vec![SparseEntry {
                part: 3,
                off: 0,
                phi: (0..pp.what_len())
                    .map(|c_| key.entry(KeyDomain::W, k, c_))
                    .collect(),
            }],
            a: vec![],
            b: pi.u2.polys[k].clone(),
        });
    }

    /* <a_l, z> = sum_i c_i w_i */
    let xn = pi.x.pow(DEGREE);
    let mut phi0 = Vec::with_capacity(m * pp.f);
    for k in 0..pp.f {
        let mut s = Zq::new(2).pow(k * pp.b as usize);
        for _ in 0..m {
            let mut p = Rq::zero();
            p.coeffs[0] = s;
            phi0.push(p);
            s = s * xn;
        }
    }
    let bu_scale = Zq::new(2).pow(pp.bu as usize);
    let phi1: Vec<Rq> = phi0.iter().map(|p| p * bu_scale).collect();
    let mut phiw = vec![Rq::zero(); pp.what_len()];
    for d in 0..pp.fu {
        let s = Zq::new(2).pow(d * pp.bu as usize);
        for i in 0..n {
            phiw[d * n + i] = -&(&c[i] * s);
        }
    }
    cnst.push(SparseConstraint {
        ct_only: false,
        entries: vec![
            SparseEntry {
                part: 0,
                off: 0,
                phi: phi0,
            },
            SparseEntry {
                part: 1,
                off: 0,
                phi: phi1,
            },
            SparseEntry {
                part: 3,
                off: 0,
                phi: phiw,
            },
        ],
        a: vec![],
        b: Rq::zero(),
    });

    /* A z = sum_i c_i t_i : kappa ring constraints */
    for k in 0..pp.kappa {
        let arow: Vec<Rq> = (0..pp.zhat_len())
            .map(|c_| key.entry(KeyDomain::A, k, c_))
            .collect();
        let arow1: Vec<Rq> = arow.iter().map(|p| p * bu_scale).collect();
        let mut phit = vec![Rq::zero(); pp.that_len()];
        for i in 0..n {
            for d in 0..pp.fu {
                let s = Zq::new(2).pow(d * pp.bu as usize);
                phit[i * l + d * pp.kappa + k] = -&(&c[i] * s);
            }
        }
        cnst.push(SparseConstraint {
            ct_only: false,
            entries: vec![
                SparseEntry {
                    part: 0,
                    off: 0,
                    phi: arow,
                },
                SparseEntry {
                    part: 1,
                    off: 0,
                    phi: arow1,
                },
                SparseEntry {
                    part: 2,
                    off: 0,
                    phi: phit,
                },
            ],
            a: vec![],
            b: Rq::zero(),
        });
    }

    /* ct(<a_r, w-hat>) = y */
    // xvec: ct(xvec_poly * w) = w(x): xvec[0] = 1, xvec[i] = -x^(N-i)
    let mut xvec = [Zq::new(0); DEGREE];
    xvec[0] = Zq::new(1);
    xvec[DEGREE - 1] = -pi.x;
    for i in (1..DEGREE - 1).rev() {
        xvec[i] = xvec[i + 1] * pi.x;
    }
    let xnm = pi.x.pow(DEGREE * m); // x^(N m)
    let mut phiy = vec![Rq::zero(); pp.what_len()];
    for d in 0..pp.fu {
        let sd = Zq::new(2).pow(d * pp.bu as usize);
        let mut sj = Zq::new(1);
        for j in 0..n {
            let mut coeffs = [Zq::new(0); DEGREE];
            for (k, &xv) in xvec.iter().enumerate() {
                coeffs[k] = xv * sd * sj;
            }
            phiy[d * n + j] = Rq::from_coeffs(coeffs);
            sj = sj * xnm;
        }
    }
    cnst.push(SparseConstraint {
        ct_only: true,
        entries: vec![SparseEntry {
            part: 3,
            off: 0,
            phi: phiy,
        }],
        a: vec![],
        b: ct_poly(pi.y),
    });

    Ok(PrincipalStatement {
        n: vec![pp.zhat_len(), pp.zhat_len(), pp.that_len(), pp.what_len()],
        cnst,
        betasq: pi.normsq,
        h,
    })
}

/// Full pipeline: evaluation proof compressed with the labrador chain.
pub fn prove_eval_compressed(
    key: &ComKey,
    ctx: &CommitCtx,
    x: Zq,
    y: Zq,
) -> Result<(EvalProof, labrador::CompositeProof), String> {
    let (pi, wt) = eval(key, ctx, x, y);
    let pst = reduce(key, &pi)?;
    labrador::principal_verify(&pst, &wt).map_err(|e| format!("self-check failed: {e}"))?;
    let comp = labrador::composite_prove(key, &pst, &wt)?;
    Ok((pi, comp))
}

/// Verify a compressed evaluation proof.
pub fn verify_eval_compressed(
    key: &ComKey,
    pi: &EvalProof,
    comp: &labrador::CompositeProof,
) -> Result<(), String> {
    let pst = reduce(key, pi)?;
    labrador::composite_verify(key, &pst, comp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::sampling::uniform::poly_uniform;

    fn random_polyvec(len: usize, seed: u8) -> Vec<Rq> {
        (0..len)
            .map(|i| poly_uniform(&[seed; 32], i as u64))
            .collect()
    }

    #[test]
    fn eval_poly_matches_horner_on_monomial() {
        // f = X^(N + 2): eval at x -> x^66
        let mut v = vec![Rq::zero(); 2];
        v[1] = Rq::monomial(2);
        let x = Zq::new(43);
        assert_eq!(eval_poly(&v, x), x.pow(66));
    }

    #[test]
    fn commit_eval_uncompressed_roundtrip() {
        // len = 8 ring elements = degree < 512 polynomial
        let key = ComKey::new([77u8; 32]);
        let pv = random_polyvec(8, 1);
        let ctx = commit(&key, &pv).unwrap();
        let x = Zq::new(43);
        let y = eval_poly(&pv, x);

        let (pi, wt) = eval(&key, &ctx, x, y);
        let pst = reduce(&key, &pi).unwrap();
        labrador::principal_verify(&pst, &wt).unwrap();
    }

    #[test]
    fn wrong_value_rejected() {
        let key = ComKey::new([78u8; 32]);
        let pv = random_polyvec(8, 2);
        let ctx = commit(&key, &pv).unwrap();
        let x = Zq::new(7);
        let y = eval_poly(&pv, x);
        let wrong = y + Zq::new(1);

        // prover honestly answers for the wrong claim: the ct constraint
        // must fail
        let (pi, wt) = eval(&key, &ctx, x, wrong);
        let pst = reduce(&key, &pi).unwrap();
        assert!(labrador::principal_verify(&pst, &wt).is_err());
    }

    #[test]
    fn compressed_roundtrip() {
        let key = ComKey::new([79u8; 32]);
        let pv = random_polyvec(8, 3);
        let ctx = commit(&key, &pv).unwrap();
        let x = Zq::new(43);
        let y = eval_poly(&pv, x);

        let (pi, comp) = prove_eval_compressed(&key, &ctx, x, y).unwrap();
        verify_eval_compressed(&key, &pi, &comp).unwrap();
        assert!(!comp.rounds.is_empty());

        // tampered evaluation claim fails
        let mut bad = pi.clone();
        bad.y += Zq::new(1);
        assert!(verify_eval_compressed(&key, &bad, &comp).is_err());
    }
}
