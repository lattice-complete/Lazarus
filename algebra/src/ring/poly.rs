//! Fixed-degree polynomials over R_q = Z_q[X]/(X^D + 1).
//!
//! This is the modern replacement for the loosely-typed `PolynomialRing`
//! (variable-length coefficient vector): the ring degree lives in the type,
//! storage is a flat array, and multiplication runs through the CRT/NTT
//! engine. Ported from lazer's `poly.c` / labrador's `poly.c`+`polz.c`,
//! minus the AVX-512 layout.

use crate::ntt;
use crate::polynomial_ring::PolynomialRing;
use crate::zq::Zq;
use std::ops::{Add, Mul, Neg, Sub};

/// Dense polynomial with exactly D coefficients (degree bound D, X^D = -1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly<const D: usize> {
    pub coeffs: [Zq; D],
}

/// The LaBRADOR ring: degree 64, q = 2^32 - 99.
pub type Rq64 = Poly<64>;

impl<const D: usize> Default for Poly<D> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<const D: usize> Poly<D> {
    pub const DEGREE: usize = D;

    pub fn zero() -> Self {
        Poly {
            coeffs: [Zq { value: 0 }; D],
        }
    }

    pub fn one() -> Self {
        let mut p = Self::zero();
        p.coeffs[0] = Zq::new(1);
        p
    }

    /// X^i (i < D).
    pub fn monomial(i: usize) -> Self {
        let mut p = Self::zero();
        p.coeffs[i] = Zq::new(1);
        p
    }

    pub fn from_coeffs(coeffs: [Zq; D]) -> Self {
        Poly { coeffs }
    }

    /// Build from a slice of length <= D, zero-padding the tail.
    pub fn from_slice(c: &[Zq]) -> Self {
        assert!(c.len() <= D, "coefficient slice longer than ring degree");
        let mut p = Self::zero();
        p.coeffs[..c.len()].copy_from_slice(c);
        p
    }

    /// Build from centered representatives.
    pub fn from_centered(c: &[i64; D]) -> Self {
        let mut p = Self::zero();
        for (o, &v) in p.coeffs.iter_mut().zip(c.iter()) {
            *o = Zq::from_i64(v);
        }
        p
    }

    /// Centered representatives in (-q/2, q/2].
    pub fn centered(&self) -> [i64; D] {
        let mut out = [0i64; D];
        for (o, c) in out.iter_mut().zip(self.coeffs.iter()) {
            *o = c.centered();
        }
        out
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.value == 0)
    }

    /// Product via the CRT/NTT engine (exact integer negacyclic product,
    /// then reduction mod q). Requires D to be a power of two.
    pub fn mul_ntt(&self, other: &Self) -> Self {
        let a = self.centered();
        let b = other.centered();
        let prod = ntt::negacyclic_mul_centered(&a, &b);
        let mut out = Self::zero();
        let q = Zq::Q as i128;
        for (o, &v) in out.coeffs.iter_mut().zip(prod.iter()) {
            let mut r = v % q;
            if r < 0 {
                r += q;
            }
            *o = Zq { value: r as usize };
        }
        out
    }

    /// Reference schoolbook product; kept public for cross-testing.
    pub fn mul_schoolbook(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for i in 0..D {
            let a = self.coeffs[i];
            if a.value == 0 {
                continue;
            }
            for j in 0..D {
                let prod = a * other.coeffs[j];
                if i + j < D {
                    out.coeffs[i + j] += prod;
                } else {
                    out.coeffs[i + j - D] += -prod;
                }
            }
        }
        out
    }

    /// Conjugation automorphism sigma(a)(X) = a(X^{-1}).
    /// Since X^{-i} = -X^{D-i} in Z_q[X]/(X^D+1):
    /// result[0] = a[0], result[D-i] = -a[i] for i >= 1.
    pub fn conjugate(&self) -> Self {
        let mut out = Self::zero();
        out.coeffs[0] = self.coeffs[0];
        for i in 1..D {
            out.coeffs[D - i] = -self.coeffs[i];
        }
        out
    }

    /// Infinity norm of the centered representative.
    pub fn linf_norm(&self) -> u64 {
        self.coeffs
            .iter()
            .map(|c| c.centered().unsigned_abs())
            .max()
            .unwrap_or(0)
    }

    /// Squared l2 norm of the centered representative.
    pub fn l2_norm_squared(&self) -> u128 {
        self.coeffs
            .iter()
            .map(|c| {
                let v = c.centered() as i128;
                (v * v) as u128
            })
            .sum()
    }

    /// Balanced base-b decomposition into `parts` polynomials with
    /// coefficients in [-b/2, b/2]: self = sum_i out[i] * b^i.
    /// (The gadget decomposition used throughout LaBRADOR/lazer.)
    ///
    /// Panics if `parts` is insufficient to absorb every coefficient.
    pub fn decompose(&self, b: u64, parts: usize) -> Vec<Self> {
        assert!(b >= 2);
        let bi = b as i64;
        let mut rem = self.centered();
        let mut out = Vec::with_capacity(parts);
        for _ in 0..parts {
            let mut digits = [0i64; D];
            for (dig, r) in digits.iter_mut().zip(rem.iter_mut()) {
                let mut d = r.rem_euclid(bi);
                // Balanced digit in [-b/2, b/2]; on the b/2 tie (even b) take
                // the sign of the remainder so the quotient always shrinks.
                if d > bi / 2 || (2 * d == bi && *r < 0) {
                    d -= bi;
                }
                *dig = d;
                *r = (*r - d) / bi;
            }
            out.push(Self::from_centered(&digits));
        }
        assert!(
            rem.iter().all(|&r| r == 0),
            "decompose: {parts} parts too few for base {b}"
        );
        out
    }

    /// Like `decompose`, but the topmost part absorbs the remaining
    /// quotient instead of asserting it is zero (the labrador C library's
    /// `polxvec_decompose` convention): parts 0..t-1 are balanced digits
    /// in [-b/2, b/2], part t-1 is the exact remainder and may be larger.
    /// Always satisfies self = sum_i out[i] * b^i.
    pub fn decompose_topmost(&self, b: u64, parts: usize) -> Vec<Self> {
        assert!(b >= 2 && parts >= 1);
        let bi = b as i64;
        let mut rem = self.centered();
        let mut out = Vec::with_capacity(parts);
        for _ in 0..parts - 1 {
            let mut digits = [0i64; D];
            for (dig, r) in digits.iter_mut().zip(rem.iter_mut()) {
                let mut d = r.rem_euclid(bi);
                if d > bi / 2 || (2 * d == bi && *r < 0) {
                    d -= bi;
                }
                *dig = d;
                *r = (*r - d) / bi;
            }
            out.push(Self::from_centered(&digits));
        }
        out.push(Self::from_centered(&rem));
        out
    }

    /// Recompose from a balanced decomposition (inverse of `decompose`).
    pub fn recompose(parts: &[Self], b: u64) -> Self {
        let mut acc = Self::zero();
        let base = Zq::new(b as usize);
        let mut pow = Zq::new(1);
        for p in parts {
            acc = &acc + &(p * pow);
            pow = pow * base;
        }
        acc
    }

    /// Serialize as D little-endian u32 words (q < 2^32).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(4 * D);
        for c in &self.coeffs {
            out.extend_from_slice(&(c.value as u32).to_le_bytes());
        }
        out
    }

    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 4 * D {
            return None;
        }
        let mut p = Self::zero();
        for (o, ch) in p.coeffs.iter_mut().zip(bytes.chunks_exact(4)) {
            let v = u32::from_le_bytes(ch.try_into().ok()?) as usize;
            if v >= Zq::Q {
                return None;
            }
            *o = Zq { value: v };
        }
        Some(p)
    }
}

// ---- operator impls (by reference; owned variants delegate) ----

impl<const D: usize> Add for &Poly<D> {
    type Output = Poly<D>;
    fn add(self, rhs: &Poly<D>) -> Poly<D> {
        let mut out = self.clone();
        for (o, r) in out.coeffs.iter_mut().zip(rhs.coeffs.iter()) {
            *o += *r;
        }
        out
    }
}

impl<const D: usize> Sub for &Poly<D> {
    type Output = Poly<D>;
    fn sub(self, rhs: &Poly<D>) -> Poly<D> {
        let mut out = self.clone();
        for (o, r) in out.coeffs.iter_mut().zip(rhs.coeffs.iter()) {
            *o = *o - *r;
        }
        out
    }
}

impl<const D: usize> Neg for &Poly<D> {
    type Output = Poly<D>;
    fn neg(self) -> Poly<D> {
        let mut out = self.clone();
        for o in out.coeffs.iter_mut() {
            *o = -*o;
        }
        out
    }
}

impl<const D: usize> Mul for &Poly<D> {
    type Output = Poly<D>;
    fn mul(self, rhs: &Poly<D>) -> Poly<D> {
        if D.is_power_of_two() && D >= 2 {
            self.mul_ntt(rhs)
        } else {
            self.mul_schoolbook(rhs)
        }
    }
}

impl<const D: usize> Mul<Zq> for &Poly<D> {
    type Output = Poly<D>;
    fn mul(self, rhs: Zq) -> Poly<D> {
        let mut out = self.clone();
        for o in out.coeffs.iter_mut() {
            *o = *o * rhs;
        }
        out
    }
}

impl<const D: usize> Add for Poly<D> {
    type Output = Poly<D>;
    fn add(self, rhs: Poly<D>) -> Poly<D> {
        &self + &rhs
    }
}

impl<const D: usize> Sub for Poly<D> {
    type Output = Poly<D>;
    fn sub(self, rhs: Poly<D>) -> Poly<D> {
        &self - &rhs
    }
}

impl<const D: usize> Mul for Poly<D> {
    type Output = Poly<D>;
    fn mul(self, rhs: Poly<D>) -> Poly<D> {
        &self * &rhs
    }
}

// ---- bridges to the legacy PolynomialRing type ----

impl<const D: usize> From<&PolynomialRing> for Poly<D> {
    fn from(p: &PolynomialRing) -> Self {
        Poly::from_slice(&p.coefficients)
    }
}

impl<const D: usize> From<&Poly<D>> for PolynomialRing {
    fn from(p: &Poly<D>) -> Self {
        PolynomialRing {
            coefficients: p.coeffs.to_vec(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    fn random_poly(rng: &mut impl Rng) -> Rq64 {
        let mut p = Rq64::zero();
        for c in p.coeffs.iter_mut() {
            *c = Zq::new(rng.gen_range(0..Zq::Q));
        }
        p
    }

    #[test]
    fn ntt_mul_matches_schoolbook() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        for _ in 0..10 {
            let a = random_poly(&mut rng);
            let b = random_poly(&mut rng);
            assert_eq!(a.mul_ntt(&b), a.mul_schoolbook(&b));
        }
    }

    #[test]
    fn mul_matches_legacy_polynomial_ring() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(2);
        let a = random_poly(&mut rng);
        let b = random_poly(&mut rng);
        let legacy_a = PolynomialRing::from(&a);
        let legacy_b = PolynomialRing::from(&b);
        let modern = &a * &b;
        let legacy = &legacy_a * &legacy_b;
        assert_eq!(PolynomialRing::from(&modern), legacy);
    }

    #[test]
    fn ring_identities() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(3);
        let a = random_poly(&mut rng);
        // a * 1 = a
        assert_eq!(&a * &Rq64::one(), a);
        // a * X^64 = -a  <=>  a * X * X^63 = -a
        let x = Rq64::monomial(1);
        let x63 = Rq64::monomial(63);
        let ax64 = &(&a * &x) * &x63;
        assert_eq!(ax64, -(&a));
        // conjugation is an involution
        assert_eq!(a.conjugate().conjugate(), a);
        // conjugation is multiplicative: sigma(ab) = sigma(a) sigma(b)
        let b = random_poly(&mut rng);
        assert_eq!((&a * &b).conjugate(), &a.conjugate() * &b.conjugate());
    }

    #[test]
    fn conjugate_constant_term_of_product_is_inner_product() {
        // <a, b> = ct(sigma(a) * b): the identity LaBRADOR uses to encode
        // inner products of coefficient vectors as ring products.
        let mut rng = rand::rngs::StdRng::seed_from_u64(4);
        let a = random_poly(&mut rng);
        let b = random_poly(&mut rng);
        let expected: Zq = a
            .coeffs
            .iter()
            .zip(b.coeffs.iter())
            .map(|(x, y)| *x * *y)
            .sum();
        let got = (&a.conjugate() * &b).coeffs[0];
        assert_eq!(got, expected);
    }

    #[test]
    fn decompose_recompose_roundtrip() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        let a = random_poly(&mut rng);
        for b in [2u64, 16, 1 << 12] {
            let parts = 40; // plenty for q < 2^32 even at base 2
            let dec = a.decompose(b, parts);
            assert!(dec.iter().all(|p| p.linf_norm() <= b / 2));
            assert_eq!(Rq64::recompose(&dec, b), a, "base {b}");
        }
    }

    #[test]
    fn decompose_topmost_roundtrip() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(51);
        let a = random_poly(&mut rng);
        for (b, parts) in [(1u64 << 8, 2usize), (1 << 12, 3), (1 << 30, 1)] {
            let dec = a.decompose_topmost(b, parts);
            assert_eq!(dec.len(), parts);
            for d in &dec[..parts - 1] {
                assert!(d.linf_norm() <= b / 2);
            }
            assert_eq!(Rq64::recompose(&dec, b), a, "base {b} parts {parts}");
        }
        // parts = 1 is the identity
        assert_eq!(a.decompose_topmost(1 << 10, 1)[0], a);
    }

    #[test]
    fn serialization_roundtrip() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(6);
        let a = random_poly(&mut rng);
        let bytes = a.to_bytes();
        assert_eq!(bytes.len(), 256);
        assert_eq!(Rq64::from_bytes(&bytes).unwrap(), a);
        assert!(Rq64::from_bytes(&bytes[1..]).is_none());
    }

    #[test]
    fn norms() {
        let mut c = [0i64; 64];
        c[0] = -3;
        c[5] = 2;
        let p = Rq64::from_centered(&c);
        assert_eq!(p.linf_norm(), 3);
        assert_eq!(p.l2_norm_squared(), 13);
    }
}
