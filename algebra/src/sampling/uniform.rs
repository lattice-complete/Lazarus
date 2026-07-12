//! Uniform samplers (port of lazer's `urandom.c`).

use crate::ring::{Poly, PolyVec};
use crate::rng::{DefaultRng, Xof};
use crate::zq::Zq;

/// Uniform element of Z_q by rejection on 32-bit words.
/// q = 2^32 - 99, so the rejection probability is ~2^-25.4 per draw.
pub fn zq_uniform<X: Xof>(xof: &mut X) -> Zq {
    loop {
        let x = xof.next_u32() as usize;
        if x < Zq::Q {
            return Zq { value: x };
        }
    }
}

/// Uniform integer in [0, m) by bit-truncated rejection
/// (the `_urandom_i64` strategy: draw ceil(log2 m) bits, reject >= m).
pub fn u64_uniform_mod<X: Xof>(xof: &mut X, m: u64) -> u64 {
    assert!(m > 0);
    if m == 1 {
        return 0;
    }
    let bits = 64 - (m - 1).leading_zeros();
    let mask = if bits == 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };
    loop {
        let x = xof.next_u64() & mask;
        if x < m {
            return x;
        }
    }
}

/// Uniform ring element from (seed, dom).
pub fn poly_uniform<const D: usize>(seed: &[u8; 32], dom: u64) -> Poly<D> {
    let mut xof = DefaultRng::new(seed, dom);
    let mut p = Poly::zero();
    for c in p.coeffs.iter_mut() {
        *c = zq_uniform(&mut xof);
    }
    p
}

/// Uniform vector of n ring elements from (seed, dom).
pub fn polyvec_uniform<const D: usize>(n: usize, seed: &[u8; 32], dom: u64) -> PolyVec<D> {
    let mut xof = DefaultRng::new(seed, dom);
    let mut v = PolyVec::zero(n);
    for p in v.polys.iter_mut() {
        for c in p.coeffs.iter_mut() {
            *c = zq_uniform(&mut xof);
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Shake128Rng;

    #[test]
    fn zq_uniform_in_range_and_deterministic() {
        let seed = [1u8; 32];
        let mut a = Shake128Rng::new(&seed, 0);
        let mut b = Shake128Rng::new(&seed, 0);
        for _ in 0..1000 {
            let x = zq_uniform(&mut a);
            assert!(x.value < Zq::Q);
            assert_eq!(x, zq_uniform(&mut b));
        }
    }

    #[test]
    fn u64_uniform_mod_bounds() {
        let mut xof = Shake128Rng::new(&[2u8; 32], 7);
        for m in [1u64, 2, 3, 5, 100, 1 << 33, u64::MAX] {
            for _ in 0..100 {
                assert!(u64_uniform_mod(&mut xof, m) < m);
            }
        }
    }

    #[test]
    fn poly_uniform_domain_separation() {
        let seed = [3u8; 32];
        let a: Poly<64> = poly_uniform(&seed, 0);
        let b: Poly<64> = poly_uniform(&seed, 0);
        let c: Poly<64> = poly_uniform(&seed, 1);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn uniform_mean_is_centered() {
        // Mean of centered representatives should be ~0 relative to q.
        let v: PolyVec<64> = polyvec_uniform(16, &[4u8; 32], 0);
        let n = (16 * 64) as f64;
        let mean: f64 = v
            .polys
            .iter()
            .flat_map(|p| p.centered())
            .map(|x| x as f64)
            .sum::<f64>()
            / n;
        // std of the mean is q/sqrt(12 n) ~ 4.2e7 for n=1024
        assert!(mean.abs() < 2e8, "suspicious uniform mean: {mean}");
    }
}
