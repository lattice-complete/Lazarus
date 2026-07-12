//! Discrete Gaussian sampler (port of lazer's `grandom.c`).
//!
//! Structure (FALCON-style):
//! * a base half-Gaussian with sigma_0 = 1.55, sampled by 128-bit CDF
//!   table lookup (`CDF155`, constants taken verbatim from grandom.c),
//! * bimodal shift k <- b + (2b-1)k_0 for a random sign bit b,
//! * Bernoulli rejection `BerExp` to correct to the target center c,
//! * scaling: for sigma = 1.55 * 2^log2o, draw the low bits of the output
//!   uniformly (the fractional center) and the high part from the base
//!   sampler: z = (k << log2o) - urand with c = urand / 2^log2o.
//!
//! Difference from C: grandom.c caches sign bits in a function-static
//! variable that leaks state across calls; here the cache lives in the
//! sampler instance and is reset per call, which is thread-safe and
//! reproducible in isolation (documented in docs/MIGRATION.md).

use crate::ring::Poly;
use crate::rng::{DefaultRng, Xof};

/// CDF of the half-Gaussian with sigma = 1.55, in 128-bit fixed point,
/// split as (hi, lo) u64 pairs. Verbatim from lazer `src/grandom.c`.
const CDF155: [(u64, u64); 21] = [
    (10894764499197476522, 10804844707381617341),
    (4761708367981796450, 6209732027000382074),
    (1476784279527800432, 14108379346150813303),
    (316388870594767345, 17827298407763885637),
    (46043503515468600, 18385899657892021654),
    (4503729779335039, 3860889375818664979),
    (294122444862326, 13947176349836216550),
    (12769598070895, 14321894682751135119),
    (367552986472, 10286761328368440884),
    (7001273393, 2287787188970898528),
    (88153536, 17843977990435837663),
    (733119, 12174894787802692461),
    (4024, 18067426722645776197),
    (14, 10764017821655913055),
    (0, 643125733022530080),
    (0, 1014291014134832),
    (0, 1055215183460),
    (0, 724109373),
    (0, 327744),
    (0, 98),
    (0, 0),
];

/// exp(x) for |x| <= 0.5 ln 2, via the fdlibm-derived Remez kernel
/// used by grandom.c.
#[allow(clippy::excessive_precision)] // fdlibm literals kept verbatim for auditability
fn exp_small(x: f64) -> f64 {
    const C1: f64 = 1.66666666666666019037e-01;
    const C2: f64 = -2.77777777770155933842e-03;
    const C3: f64 = 6.61375632143793436117e-05;
    const C4: f64 = -1.65339022054652515390e-06;
    const C5: f64 = 4.13813679705723846039e-08;

    let t = x * x;
    let t = x - t * (C1 + t * (C2 + t * (C3 + t * (C4 + t * C5))));
    1.0 - ((x * t) / (t - 2.0) - x)
}

/// Base sampler: half-Gaussian(1.55) by lexicographic 128-bit CDF walk.
fn cdf_sampler<X: Xof>(xof: &mut X) -> u32 {
    let mut buf = [0u8; 16];
    xof.fill(&mut buf);
    let hi = u64::from_le_bytes(buf[..8].try_into().unwrap());
    let lo = u64::from_le_bytes(buf[8..].try_into().unwrap());

    let mut z = 0usize;
    while hi < CDF155[z].0 || (hi == CDF155[z].0 && lo < CDF155[z].1) {
        z += 1;
    }
    z as u32
}

/// Bernoulli(exp(-x)) for x >= 0, using u random bits + a 53-bit compare.
fn ber_exp<X: Xof>(xof: &mut X, x: f64) -> bool {
    debug_assert!(x >= 0.0);
    let mut buf = [0u8; 16];
    xof.fill(&mut buf);
    let w0 = u64::from_le_bytes(buf[..8].try_into().unwrap());
    let w1 = u64::from_le_bytes(buf[8..].try_into().unwrap());

    // Split x = u ln2 + r with r in [0, ln2).
    let mut u = (x * std::f64::consts::LOG2_E) as u64;
    let r = x - std::f64::consts::LN_2 * u as f64;
    if u > 63 {
        u = 63;
    }
    // b1 = 1 with probability 2^-u (iff the low u bits of w0 are all zero;
    // C: t ^= (t >> u) << u; b = 1 - ((t | -t) >> 63)).
    let low = if u == 0 { 0 } else { w0 & ((1u64 << u) - 1) };
    let b1 = low == 0;

    // b2 = 1 with probability exp(-r).
    let t53 = w1 & ((1u64 << 53) - 1);
    let thresh = (exp_small(-r) * (1u64 << 53) as f64) as u64;
    let b2 = t53 < thresh;

    b1 && b2
}

/// Streaming discrete Gaussian sampler over one XOF.
pub struct GaussianSampler<X: Xof> {
    xof: X,
    bits: u64,
    pos: u32,
}

impl<X: Xof> GaussianSampler<X> {
    pub fn new(xof: X) -> Self {
        GaussianSampler {
            xof,
            bits: 0,
            pos: 64,
        }
    }

    fn sign_bit(&mut self) -> u32 {
        if self.pos >= 64 {
            self.bits = self.xof.next_u64();
            self.pos = 0;
        }
        let b = (self.bits & 1) as u32;
        self.bits >>= 1;
        self.pos += 1;
        b
    }

    /// Sample z ~ D_{Z, 1.55, c} for c in [0, 1) (grandom.c `gaussian155`).
    pub fn sample_base(&mut self, c: f64) -> i32 {
        const DSS: f64 = 1.0 / (2.0 * 1.55 * 1.55);
        loop {
            let b = self.sign_bit() as i32;
            let k0 = cdf_sampler(&mut self.xof) as i32;
            // bimodal shift: b=1 -> k0+1, b=0 -> -k0
            let k = b + (2 * b - 1) * k0;
            let kf = k as f64;
            let bf = b as f64;
            let x = ((kf - c) * (kf - c) - (kf - bf) * (kf - bf)) * DSS;
            if ber_exp(&mut self.xof, x) {
                return k;
            }
        }
    }
}

/// n samples of the centered discrete Gaussian with sigma = 1.55 * 2^log2o,
/// deterministic in (seed, dom). Port of `_grandom_sample_i32`, widened to
/// i64 (safe for log2o up to ~45).
pub fn gaussian_vec(n: usize, log2o: u32, seed: &[u8; 32], dom: u64) -> Vec<i64> {
    assert!(log2o <= 45, "sigma too large for i64 output");
    let mut xof = DefaultRng::new(seed, dom);

    // The C code draws all fractional-center bits first, then samples.
    let outlen = (log2o as usize * n).div_ceil(8);
    let mut stream = vec![0u8; outlen];
    if outlen > 0 {
        xof.fill(&mut stream);
    }

    let mut sampler = GaussianSampler::new(xof);
    let mut out = vec![0i64; n];
    for (i, o) in out.iter_mut().enumerate() {
        let mut urand: u64 = 0;
        for j in 0..log2o as usize {
            let idx = i * log2o as usize + j;
            let bit = ((stream[idx >> 3] >> (idx & 7)) & 1) as u64;
            urand |= bit << j;
        }
        let c = urand as f64 / (1u64 << log2o) as f64;
        let k = sampler.sample_base(c) as i64;
        *o = (k << log2o) - urand as i64;
    }
    out
}

/// Ring element with discrete Gaussian coefficients, sigma = 1.55 * 2^log2o.
pub fn poly_gaussian<const D: usize>(log2o: u32, seed: &[u8; 32], dom: u64) -> Poly<D> {
    let v = gaussian_vec(D, log2o, seed, dom);
    let mut centered = [0i64; D];
    centered.copy_from_slice(&v);
    Poly::from_centered(&centered)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exp_small_accuracy() {
        for x in [-0.34, -0.2, -0.05, 0.0, 0.05, 0.2, 0.34] {
            let got = exp_small(x);
            let want = x.exp();
            assert!(
                (got - want).abs() < 1e-14,
                "exp_small({x}) = {got} vs {want}"
            );
        }
    }

    #[test]
    fn gaussian_deterministic() {
        let seed = [13u8; 32];
        assert_eq!(
            gaussian_vec(256, 3, &seed, 1),
            gaussian_vec(256, 3, &seed, 1)
        );
        assert_ne!(
            gaussian_vec(256, 3, &seed, 1),
            gaussian_vec(256, 3, &seed, 2)
        );
    }

    #[test]
    fn gaussian_moments() {
        // sigma = 1.55 * 2^3 = 12.4; check sample mean/variance on a
        // fixed-seed draw (deterministic, so no flakiness).
        let n = 1 << 14;
        let v = gaussian_vec(n, 3, &[14u8; 32], 0);
        let sigma = 1.55 * 8.0;
        let mean = v.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
        let var = v.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.5, "gaussian mean off: {mean}");
        let expect = sigma * sigma;
        assert!(
            (var - expect).abs() / expect < 0.1,
            "gaussian variance off: {var} vs {expect}"
        );
        // tail sanity: nothing should stray beyond 8 sigma at n=2^14
        assert!(v.iter().all(|&x| (x as f64).abs() < 8.0 * sigma));
    }

    #[test]
    fn gaussian_sigma_one() {
        // log2o = 0: pure base sampler, sigma = 1.55, integer support.
        let v = gaussian_vec(4096, 0, &[15u8; 32], 0);
        let var = v.iter().map(|&x| (x * x) as f64).sum::<f64>() / v.len() as f64;
        assert!(
            (var - 1.55 * 1.55).abs() / (1.55 * 1.55) < 0.15,
            "var {var}"
        );
    }
}
