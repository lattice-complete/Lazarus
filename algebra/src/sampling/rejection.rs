//! Rejection sampling (port of lazer's `rejection.c`).
//!
//! Used after committing to a masked response z = y + v (or z = y + (-1)^b v
//! in the bimodal variant) to make the distribution of z independent of the
//! secret-dependent shift v.
//!
//! The C code compares a 128-bit uniform integer u against MPFR 128-bit
//! transcendentals. This portable port performs the comparison in the log
//! domain with f64 (~2^-52 relative error), which perturbs the acceptance
//! probability by a statistically negligible amount but is NOT bit-identical
//! to the C sampler; see docs/MIGRATION.md.

use crate::rng::Xof;

fn dot(a: &[i64], b: &[i64]) -> i128 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| x as i128 * y as i128)
        .sum()
}

/// ln(u / 2^128) for a 128-bit uniform integer u (None for u = 0,
/// meaning "always accept").
fn ln_unit(u: u128) -> Option<f64> {
    if u == 0 {
        return None;
    }
    // f64 conversion of u128 keeps 53 significant bits — ample here.
    Some((u as f64).ln() - 128.0 * std::f64::consts::LN_2)
}

/// Standard rejection: accept z with probability
///   min(1, exp((-2<z,v> + <v,v>) / (2 sigma^2)) / M).
///
/// Returns true iff accepted. `m` is the repetition constant M > 1.
pub fn rej_standard<X: Xof>(xof: &mut X, z: &[i64], v: &[i64], m: f64, sigma: f64) -> bool {
    assert!(m >= 1.0 && sigma > 0.0);
    let u = xof.next_u128();
    let x = (dot(v, v) - 2 * dot(z, v)) as f64 / (2.0 * sigma * sigma);
    match ln_unit(u) {
        None => true,
        // accept iff u/2^128 <= exp(x)/M
        Some(ln_u) => ln_u + m.ln() <= x,
    }
}

/// Bimodal rejection (BLISS-style): accept with probability
///   1 / (M exp(-<v,v> / (2 sigma^2)) cosh(<z,v> / sigma^2)).
pub fn rej_bimodal<X: Xof>(xof: &mut X, z: &[i64], v: &[i64], m: f64, sigma: f64) -> bool {
    assert!(m >= 1.0 && sigma > 0.0);
    let u = xof.next_u128();
    let s2 = sigma * sigma;
    let vv = dot(v, v) as f64 / (2.0 * s2);
    let zv = dot(z, v) as f64 / s2;
    match ln_unit(u) {
        None => true,
        // accept iff u/2^128 <= exp(vv) / (M cosh(zv))
        Some(ln_u) => ln_u + m.ln() + ln_cosh(zv) <= vv,
    }
}

/// Overflow-safe ln(cosh(y)).
fn ln_cosh(y: f64) -> f64 {
    let a = y.abs();
    // cosh(y) = e^a (1 + e^-2a) / 2
    a + (-2.0 * a).exp().ln_1p() - std::f64::consts::LN_2
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Shake128Rng;

    #[test]
    fn zero_shift_acceptance_rate_is_one_over_m() {
        // With v = 0 the acceptance condition is u/2^128 <= 1/M exactly.
        let m = std::f64::consts::E; // rate should be ~ 1/e
        let z = [0i64; 4];
        let v = [0i64; 4];
        let mut xof = Shake128Rng::new(&[21u8; 32], 0);
        let trials = 20_000;
        let mut acc = 0;
        for _ in 0..trials {
            if rej_standard(&mut xof, &z, &v, m, 100.0) {
                acc += 1;
            }
        }
        let rate = acc as f64 / trials as f64;
        assert!(
            (rate - 1.0 / m).abs() < 0.02,
            "acceptance rate {rate} vs {}",
            1.0 / m
        );
    }

    #[test]
    fn aligned_z_accepts_more_than_opposed_z() {
        // exp((-2<z,v>+<v,v>)/(2s^2)): z aligned with v accepts less often
        // than z opposed to v — check the monotonicity both ways.
        let sigma = 50.0;
        let v: Vec<i64> = (0..8).map(|i| 10 + i).collect();
        let z_aligned: Vec<i64> = v.iter().map(|x| x * 3).collect();
        let z_opposed: Vec<i64> = v.iter().map(|x| -x * 3).collect();
        let m = 3.0;
        let (mut acc_a, mut acc_o) = (0, 0);
        let mut xof = Shake128Rng::new(&[22u8; 32], 0);
        for _ in 0..5000 {
            if rej_standard(&mut xof, &z_aligned, &v, m, sigma) {
                acc_a += 1;
            }
            if rej_standard(&mut xof, &z_opposed, &v, m, sigma) {
                acc_o += 1;
            }
        }
        assert!(
            acc_a < acc_o,
            "aligned {acc_a} should accept less than opposed {acc_o}"
        );
    }

    #[test]
    fn bimodal_symmetric_in_sign_of_zv() {
        // cosh is even: flipping the sign of z must not change acceptance.
        let sigma = 40.0;
        let v: Vec<i64> = (0..8).map(|i| 5 + i).collect();
        let z: Vec<i64> = v.iter().map(|x| x * 2).collect();
        let zneg: Vec<i64> = z.iter().map(|x| -x).collect();
        let m = 2.0;
        let mut a = Shake128Rng::new(&[23u8; 32], 7);
        let mut b = Shake128Rng::new(&[23u8; 32], 7);
        for _ in 0..2000 {
            assert_eq!(
                rej_bimodal(&mut a, &z, &v, m, sigma),
                rej_bimodal(&mut b, &zneg, &v, m, sigma)
            );
        }
    }

    #[test]
    fn ln_cosh_matches_naive() {
        for y in [-3.0f64, -0.5, 0.0, 0.1, 2.0, 20.0] {
            let want = y.cosh().ln();
            assert!((ln_cosh(y) - want).abs() < 1e-12);
        }
        // and stays finite where naive overflows
        assert!(ln_cosh(1000.0).is_finite());
    }
}
