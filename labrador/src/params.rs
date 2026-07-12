//! Protocol constants, security estimation and automatic parameter
//! selection (port of the reference labrador's `data.h` constants,
//! `sis_secure()` and the parameter search in `init_proof()`).

use crate::witness::Witness;
use algebra::sampling::challenge::{CHALLENGE_NORMSQ, OPNORM_BOUND, TAU1, TAU2};
use algebra::zq::Zq;

/// Ring degree (C: N).
pub const DEGREE: usize = 64;
pub const LOGQ: u32 = 32;
/// Parallel Z_q-constraint liftings for ~2^128 soundness: ceil(128/LOGQ).
pub const LIFTS: usize = 4;
/// Slack factor in extractability bounds (C: SLACK).
pub const SLACK: f64 = 2.0;
/// Challenge operator norm bound (C: T).
pub const OPNORM_T: f64 = OPNORM_BOUND;
/// log2 of the root Hermite factor targeted by the SIS estimator.
pub fn log_delta() -> f64 {
    1.00444f64.log2()
}
/// Second moment of a challenge polynomial (C: TAU1 + 4*TAU2 = 64).
pub const CHAL_VAR: f64 = CHALLENGE_NORMSQ as f64;
pub const _TAUS: (usize, usize) = (TAU1, TAU2);

/// Maximum witness norm the JL projection supports
/// (C jlproj.h: min(q/125, 2^28), squared).
pub fn jl_max_normsq() -> u128 {
    let a = (Zq::Q as u64) / 125;
    let b = 1u64 << 28;
    let m = a.min(b) as u128;
    m * m
}

/// MSIS hardness estimate (port of `sis_secure`): the instance with the
/// given module rank and l2 solution-norm bound is hard iff
/// log2(norm) < min(LOGQ, 2*sqrt(LOGQ * log2(delta) * d * rank)).
pub fn sis_secure(rank: usize, norm: f64) -> bool {
    let maxlog = 2.0 * (LOGQ as f64 * log_delta() * DEGREE as f64 * rank as f64).sqrt();
    let maxlog = maxlog.min(LOGQ as f64);
    norm.log2() < maxlog
}

/// Commitment / decomposition parameters for one reduction round
/// (C: comparams).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComParams {
    /// Amortized-opening decomposition parts / log2(base).
    pub f: usize,
    pub b: u32,
    /// Uniform (commitment & linear garbage) decomposition parts / log2(base).
    pub fu: usize,
    pub bu: u32,
    /// Quadratic-garbage decomposition parts / log2(base).
    pub fg: usize,
    pub bg: u32,
    /// Inner / outer commitment module ranks.
    pub kappa: usize,
    pub kappa1: usize,
    pub u1len: usize,
    pub u2len: usize,
}

/// Witness repartition plan for one round: the input witness vectors are
/// grouped (consecutive vectors between "boundary" markers) and each group
/// is split into `nu[boundary]` chunks of rank `n`, giving `r` amortization
/// chunks in total (port of the rank joining/splitting in `init_proof` /
/// `init_statement`).
#[derive(Debug, Clone)]
pub struct Plan {
    /// Input part ranks.
    pub n_in: Vec<usize>,
    /// Chunks per part; nonzero exactly at group boundaries.
    pub nu: Vec<usize>,
    /// Output multiplicity (sum of nu).
    pub r: usize,
    /// Output chunk rank.
    pub n: usize,
    /// Flat slot index (into the r*n layout) where each input part starts.
    pub part_base: Vec<usize>,
}

impl Plan {
    fn build(n_in: &[usize], markers: &[bool], k: usize) -> Plan {
        // group totals
        let mut nu = vec![0usize; n_in.len()];
        let mut group_max = 0usize;
        let mut t = 0usize;
        for i in 0..n_in.len() {
            t += n_in[i];
            if markers[i] {
                nu[i] = t;
                group_max = group_max.max(t);
                t = 0;
            }
        }
        let nn = group_max.div_ceil(k);
        for x in nu.iter_mut() {
            if *x > 0 {
                *x = x.div_ceil(nn);
            }
        }
        let r: usize = nu.iter().sum();

        // output chunk rank = max over groups of ceil(group_total / nu)
        // (init_statement's recomputation)
        let mut n = 0usize;
        let mut acc = 0usize;
        for i in 0..n_in.len() {
            acc += n_in[i];
            if nu[i] > 0 {
                n = n.max(acc.div_ceil(nu[i]));
                acc = 0;
            }
        }

        // flat base offset of each part
        let mut part_base = vec![0usize; n_in.len()];
        let mut chunk = 0usize; // group start chunk
        let mut off = 0usize; // offset within group
        for i in 0..n_in.len() {
            part_base[i] = chunk * n + off;
            off += n_in[i];
            if nu[i] > 0 {
                chunk += nu[i];
                off = 0;
            }
        }

        Plan {
            n_in: n_in.to_vec(),
            nu,
            r,
            n,
            part_base,
        }
    }

    /// Flat slot of (part, offset-within-part).
    pub fn slot(&self, part: usize, off: usize) -> usize {
        debug_assert!(off < self.n_in[part]);
        self.part_base[part] + off
    }

    /// True iff every input part occupies exactly one whole chunk (needed
    /// for index-preserving quadratic constraints).
    pub fn parts_unsplit(&self) -> bool {
        self.nu.iter().all(|&x| x == 1)
            && self
                .part_base
                .iter()
                .enumerate()
                .all(|(i, &b)| b == i * self.n)
    }
}

/// Everything `init_proof` decides for a round.
#[derive(Debug, Clone)]
pub struct RoundParams {
    pub plan: Plan,
    pub cpp: ComParams,
    pub tail: bool,
    pub quadratic: bool,
    /// Predicted output witness squared norm (drives security checks).
    pub predicted_normsq: f64,
}

/// Port of `init_proof`'s automatic parameter search (quadratic in {0,1};
/// the reference's quadratic==2 witness-doubling mode is out of scope).
pub fn select_params(wt: &Witness, quadratic: bool, tail: bool) -> Result<RoundParams, String> {
    let r_in = wt.parts.len();
    let n_in: Vec<usize> = wt.parts.iter().map(|p| p.len()).collect();
    let normsq: Vec<f64> = wt.norms().iter().map(|&x| x as f64).collect();
    let total_normsq: f64 = normsq.iter().sum();

    // group boundary markers: quadratic keeps every part separate,
    // linear joins everything into one group
    let mut markers = vec![quadratic; r_in];
    markers[r_in - 1] = true;

    let nd = DEGREE as f64;
    let mut best: Option<(Plan, ComParams, f64)> = None;

    for k in (1..=15usize).rev() {
        let plan = Plan::build(&n_in, &markers, k);
        let nn = plan.n as f64;
        let rr = plan.r;

        /* decomposition in width */
        let varz = total_normsq / (nn * nd) * CHAL_VAR;
        let mut decompose = !tail
            && !sis_secure(
                13,
                6.0 * OPNORM_T * SLACK * (2.0 * CHAL_VAR * varz * nn * nd).sqrt(),
            );
        decompose = decompose || 64.0 * varz > (1u64 << 28) as f64;
        let (f, b) = if decompose {
            (2usize, ((12f64.log2() + varz.log2()) / 4.0).round() as u32)
        } else {
            (1usize, ((12f64.log2() + varz.log2()) / 2.0).round() as u32)
        };
        let b = b.max(1);

        /* uniform decomposition */
        let (fu, bu) = if !tail {
            let fu = (LOGQ as usize + 2 * b as usize / 3) / b as usize;
            let bu = (LOGQ + fu as u32 / 2) / fu as u32;
            (fu, bu)
        } else {
            (1usize, LOGQ)
        };

        /* quadratic garbage decomposition */
        let pairs = (rr * rr + rr) / 2;
        let bg = b;
        let (fg, varg) = if !quadratic {
            (0usize, 0.0f64)
        } else if tail {
            (1usize, 0.0f64)
        } else {
            let mut varg = 0.0f64;
            let (mut t, mut u) = (0.0f64, 0usize);
            for i in 0..r_in {
                let vars = normsq[i] / (n_in[i] as f64 * nd);
                let mut j = n_in[i];
                while j >= plan.n - u {
                    j -= plan.n - u;
                    t += vars * vars * (plan.n - u) as f64;
                    varg = varg.max(t);
                    t = 0.0;
                    u = 0;
                }
                t += vars * vars * j as f64;
                u += j;
                if plan.nu[i] > 0 {
                    varg = varg.max(t);
                    t = 0.0;
                    u = 0;
                }
            }
            varg *= 2.0 * nd;
            let fg = (((12f64.log2() + varg.log2()) / (2.0 * b as f64)).ceil() as usize).max(1);
            (fg, varg)
        };

        /* commitment ranks */
        let mut kappa = 33usize;
        let mut predicted = 0.0f64;
        for kap in 1..=32usize {
            let mut p = (2f64.powi(2 * b as i32) / 12.0 * (f - 1) as f64
                + varz / 2f64.powi(2 * b as i32 * (f - 1) as i32))
                * nn;
            if !tail {
                p += (2f64.powi(2 * bu as i32) * (fu - 1) as f64
                    + 2f64.powi(2 * (LOGQ - (fu as u32 - 1) * bu) as i32))
                    / 12.0
                    * (rr * kap + pairs) as f64;
            }
            if !tail && quadratic {
                p += (2f64.powi(2 * bg as i32) / 12.0 * (fg - 1) as f64
                    + varg / 2f64.powi(2 * (fg as i32 - 1) * bg as i32))
                    * pairs as f64;
            }
            p *= nd;
            if sis_secure(
                kap,
                6.0 * OPNORM_T * SLACK * 2f64.powi(((f - 1) as u32 * b) as i32) * p.sqrt(),
            ) {
                kappa = kap;
                predicted = p;
                break;
            }
        }
        if kappa > 32 {
            continue;
        }

        let (kappa1, u1len, u2len, ok) = if !tail {
            let mut kappa1 = 33usize;
            for k1 in 1..=32usize {
                if sis_secure(k1, 2.0 * SLACK * predicted.sqrt()) {
                    kappa1 = k1;
                    break;
                }
            }
            let ok = kappa1 <= 32 && (fu * rr * kappa + (fu + fg) * pairs) as f64 <= 1.1 * nn;
            (kappa1, kappa1, kappa1, ok)
        } else {
            let u1len = rr * kappa + if quadratic { pairs } else { 0 };
            let u2len = 2 * rr - 1;
            let ok =
                ((u1len + u2len) as f64 * LOGQ as f64) <= 1.1 * nn * (varz.log2() / 2.0 + 2.05);
            (0usize, u1len, u2len, ok)
        };

        let cpp = ComParams {
            f,
            b,
            fu,
            bu,
            fg,
            bg,
            kappa,
            kappa1,
            u1len,
            u2len,
        };
        best = Some((plan, cpp, predicted));
        if ok && (tail || kappa1 <= 32) {
            break;
        }
    }

    let (plan, cpp, predicted_normsq) =
        best.ok_or_else(|| "cannot make inner commitments secure".to_string())?;
    if !tail && cpp.kappa1 > 32 {
        return Err("cannot make outer commitments secure".to_string());
    }
    Ok(RoundParams {
        plan,
        cpp,
        tail,
        quadratic,
        predicted_normsq,
    })
}

/// Rank of the auxiliary witness vector v = (t-hat, g-hat, h-hat)
/// (C: statement.m).
pub fn aux_rank(r: usize, cpp: &ComParams) -> usize {
    let pairs = (r * r + r) / 2;
    r * cpp.fu * cpp.kappa + (cpp.fu + cpp.fg) * pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::ring::PolyVec;

    #[test]
    fn sis_estimator_sanity() {
        // tiny norms are always secure at rank 1; astronomically large never
        assert!(sis_secure(1, 10.0));
        assert!(!sis_secure(32, 2f64.powi(40))); // above q: never secure
    }

    #[test]
    fn plan_join_and_split() {
        // linear: 4 parts of rank 8 join into one group of 32,
        // k = 2 -> chunks of 16
        let n_in = vec![8usize, 8, 8, 8];
        let markers = vec![false, false, false, true];
        let plan = Plan::build(&n_in, &markers, 2);
        assert_eq!(plan.n, 16);
        assert_eq!(plan.r, 2);
        assert_eq!(plan.part_base, vec![0, 8, 16, 24]);
        assert_eq!(plan.slot(2, 3), 19);

        // quadratic: each part own group, k = 1 -> one chunk per part
        let markers = vec![true; 4];
        let plan = Plan::build(&n_in, &markers, 1);
        assert_eq!((plan.r, plan.n), (4, 8));
        assert!(plan.parts_unsplit());
    }

    #[test]
    fn select_params_produces_secure_round() {
        use algebra::ring::Poly;
        // small ternary witness: 3 parts of rank 8
        let mut c = [0i64; 64];
        for (i, x) in c.iter_mut().enumerate() {
            *x = if i % 3 == 0 { 1 } else { -1 };
        }
        let p = Poly::<64>::from_centered(&c);
        let parts = vec![PolyVec::from_polys(vec![p; 8]); 3];
        let wt = Witness::new(parts);
        let rp = select_params(&wt, false, false).unwrap();
        assert!(rp.cpp.kappa >= 1 && rp.cpp.kappa <= 32);
        assert!(rp.cpp.f >= 1 && rp.cpp.fu >= 1);
        assert_eq!(rp.cpp.fg, 0);
        assert!(rp.plan.r >= 1);
    }
}
