//! LaBRADOR challenge-set sampler (port of `challenge()` in the reference
//! labrador `poly.c`).
//!
//! A challenge is a degree-64 polynomial whose coefficient multiset is
//! exactly {TAU1 x (+-1), TAU2 x (+-2), rest 0} with uniform signs and a
//! uniform placement (inside-out Fisher-Yates driven by rejection-filtered
//! XOF bytes), additionally rejected unless its operator norm is <= T.
//!
//! The operator norm of multiplication by c in R = Z[X]/(X^d + 1) equals
//! max_j |c(zeta_j)| over the primitive 2d-th roots of unity zeta_j (the
//! multiplication operator diagonalizes in the complex embedding); the
//! reference computes the same quantity with an FFT (`poly_opnorm`).

use crate::ring::Poly;
use crate::rng::Xof;

/// Number of +-1 coefficients.
pub const TAU1: usize = 32;
/// Number of +-2 coefficients.
pub const TAU2: usize = 8;
/// Operator norm rejection bound.
pub const OPNORM_BOUND: f64 = 14.0;
/// Second moment of a challenge: sum of squared coefficients.
pub const CHALLENGE_NORMSQ: u64 = (TAU1 + 4 * TAU2) as u64;

/// Operator norm of (multiplication by) a, computed on the centered
/// representative via the complex embedding: max |a(zeta_j)|,
/// zeta_j = exp(i pi (2j+1) / d), j = 0..d/2 (conjugate symmetry).
pub fn opnorm<const D: usize>(a: &Poly<D>) -> f64 {
    let c = a.centered();
    let mut max = 0.0f64;
    for j in 0..D / 2 {
        let theta = std::f64::consts::PI * (2 * j + 1) as f64 / D as f64;
        let (mut re, mut im) = (0.0f64, 0.0f64);
        for (k, &ck) in c.iter().enumerate() {
            let ang = theta * k as f64;
            re += ck as f64 * ang.cos();
            im += ck as f64 * ang.sin();
        }
        let norm = (re * re + im * im).sqrt();
        if norm > max {
            max = norm;
        }
    }
    max
}

/// Sample one challenge polynomial from the XOF (rejection-sampled on the
/// operator norm). Mirrors the reference algorithm exactly: 5 sign bytes,
/// then inside-out Fisher-Yates placement with 6-bit position bytes.
pub fn challenge_poly<X: Xof, const D: usize>(xof: &mut X) -> Poly<D> {
    assert_eq!(D, 64, "challenge multiset is calibrated for degree 64");
    loop {
        let mut signs: u64 = 0;
        let mut sb = [0u8; (TAU1 + TAU2).div_ceil(8)];
        xof.fill(&mut sb);
        for (i, &byte) in sb.iter().enumerate() {
            signs |= (byte as u64) << (8 * i);
        }

        let mut c = [0i64; D];
        let mut k = D - TAU1 - TAU2;
        while k < D {
            let mut byte = [0u8; 1];
            xof.fill(&mut byte);
            let b = (byte[0] & (D as u8 - 1)) as usize;
            if b <= k {
                c[k] = c[b];
                let mag: i64 = if k < D - TAU2 { 1 } else { 2 };
                c[b] = if signs & 1 == 1 { -mag } else { mag };
                signs >>= 1;
                k += 1;
            }
        }

        let poly = Poly::from_centered(&c);
        if opnorm(&poly) <= OPNORM_BOUND {
            return poly;
        }
    }
}

/// Sample `count` challenges.
pub fn challenge_polys<X: Xof, const D: usize>(xof: &mut X, count: usize) -> Vec<Poly<D>> {
    (0..count).map(|_| challenge_poly(xof)).collect()
}

/// "Quarternary" aggregation challenge: every coefficient uniform in
/// {-2, -1, 0, 1} (2 bits per coefficient, minus 2) — the cheap folding
/// challenges used by the reference `polyvec_quarternary`.
pub fn quarternary_poly<X: Xof, const D: usize>(xof: &mut X) -> Poly<D> {
    let mut bytes = vec![0u8; D / 4];
    xof.fill(&mut bytes);
    let mut c = [0i64; D];
    for (i, coeff) in c.iter_mut().enumerate() {
        let bits = (bytes[i / 4] >> (2 * (i % 4))) & 3;
        *coeff = bits as i64 - 2;
    }
    Poly::from_centered(&c)
}

pub fn quarternary_polys<X: Xof, const D: usize>(xof: &mut X, count: usize) -> Vec<Poly<D>> {
    (0..count).map(|_| quarternary_poly(xof)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Shake128Rng;
    use crate::zq::Zq;

    #[test]
    fn challenge_multiset_and_opnorm() {
        let mut xof = Shake128Rng::new(&[1u8; 32], 0);
        for _ in 0..8 {
            let c: Poly<64> = challenge_poly(&mut xof);
            let cent = c.centered();
            let ones = cent.iter().filter(|&&x| x.abs() == 1).count();
            let twos = cent.iter().filter(|&&x| x.abs() == 2).count();
            let zeros = cent.iter().filter(|&&x| x == 0).count();
            assert_eq!((ones, twos, zeros), (TAU1, TAU2, 64 - TAU1 - TAU2));
            assert_eq!(c.l2_norm_squared(), CHALLENGE_NORMSQ as u128);
            assert!(opnorm(&c) <= OPNORM_BOUND);
        }
    }

    #[test]
    fn opnorm_known_values() {
        // opnorm(constant k) = k; opnorm(X^j) = 1.
        let mut c = [0i64; 64];
        c[0] = 5;
        let p = Poly::<64>::from_centered(&c);
        assert!((opnorm(&p) - 5.0).abs() < 1e-9);
        let x7 = Poly::<64>::monomial(7);
        assert!((opnorm(&x7) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn opnorm_bounds_product_linf() {
        // For any a, b: linf(a*b) <= opnorm(a) * l2(b) (Cauchy-Schwarz in
        // the embedding); check on random-ish small polys.
        let mut xof = Shake128Rng::new(&[2u8; 32], 1);
        let a: Poly<64> = challenge_poly(&mut xof);
        let b: Poly<64> = quarternary_poly(&mut xof);
        let prod = &a * &b;
        let l2b = (b.l2_norm_squared() as f64).sqrt();
        assert!((prod.linf_norm() as f64) <= opnorm(&a) * l2b + 1e-6);
    }

    #[test]
    fn quarternary_range() {
        let mut xof = Shake128Rng::new(&[3u8; 32], 2);
        let p: Poly<64> = quarternary_poly(&mut xof);
        assert!(p.coeffs.iter().all(|&z| {
            let c = z.centered();
            (-2..=1).contains(&c)
        }));
        // determinism
        let mut xof2 = Shake128Rng::new(&[3u8; 32], 2);
        let p2: Poly<64> = quarternary_poly(&mut xof2);
        assert_eq!(p, p2);
        let _ = Zq::new(0);
    }
}
