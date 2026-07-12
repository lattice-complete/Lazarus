//! Centered binomial sampler (port of lazer's `brandom.c`).
//!
//! Produces values in [-k, k] distributed as the difference of two
//! Binomial(k, 1/2) draws — the standard cheap approximation of a narrow
//! discrete Gaussian (as in Kyber). The bit layout matches the C code
//! exactly: one contiguous 2k·n-bit stream, "plus" bits first
//! (bit i·k+j), then "minus" bits (bit k·n + i·k + j).

use crate::ring::Poly;
use crate::rng::{DefaultRng, Xof};

#[inline]
fn bit(stream: &[u8], idx: usize) -> i32 {
    ((stream[idx >> 3] >> (idx & 7)) & 1) as i32
}

/// n samples from CBD(k), deterministic in (seed, dom).
pub fn cbd_vec(n: usize, k: usize, seed: &[u8; 32], dom: u64) -> Vec<i32> {
    assert!((1..=127).contains(&k));
    let mut xof = DefaultRng::new(seed, dom);
    let outlen = (2 * k * n).div_ceil(8);
    let mut stream = vec![0u8; outlen];
    xof.fill(&mut stream);

    let mut out = vec![0i32; n];
    for (i, o) in out.iter_mut().enumerate() {
        for j in 0..k {
            *o += bit(&stream, i * k + j);
        }
    }
    for (i, o) in out.iter_mut().enumerate() {
        for j in 0..k {
            *o -= bit(&stream, k * n + i * k + j);
        }
    }
    out
}

/// Ring element with CBD(k) coefficients.
pub fn poly_cbd<const D: usize>(k: usize, seed: &[u8; 32], dom: u64) -> Poly<D> {
    let v = cbd_vec(D, k, seed, dom);
    let mut centered = [0i64; D];
    for (o, x) in centered.iter_mut().zip(v.iter()) {
        *o = *x as i64;
    }
    Poly::from_centered(&centered)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cbd_range_and_determinism() {
        let seed = [9u8; 32];
        for k in [1usize, 2, 3, 8] {
            let a = cbd_vec(1024, k, &seed, 3);
            let b = cbd_vec(1024, k, &seed, 3);
            assert_eq!(a, b);
            assert!(a.iter().all(|&x| x.unsigned_abs() as usize <= k));
        }
    }

    #[test]
    fn cbd_moments() {
        // CBD(k): mean 0, variance k/2.
        let k = 2usize;
        let n = 1 << 14;
        let v = cbd_vec(n, k, &[10u8; 32], 0);
        let mean: f64 = v.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
        let var: f64 = v.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.05, "CBD mean off: {mean}");
        let expect = k as f64 / 2.0;
        assert!(
            (var - expect).abs() / expect < 0.1,
            "CBD variance off: {var} vs {expect}"
        );
    }

    #[test]
    fn poly_cbd_matches_vec() {
        let seed = [11u8; 32];
        let p: Poly<64> = poly_cbd(3, &seed, 42);
        let v = cbd_vec(64, 3, &seed, 42);
        let c = p.centered();
        for i in 0..64 {
            assert_eq!(c[i], v[i] as i64);
        }
    }
}
