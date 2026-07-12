//! Vectors over R_q (port of lazer's `polyvec.c`).

use super::poly::Poly;
use crate::zq::Zq;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PolyVec<const D: usize> {
    pub polys: Vec<Poly<D>>,
}

impl<const D: usize> PolyVec<D> {
    pub fn zero(n: usize) -> Self {
        PolyVec {
            polys: vec![Poly::zero(); n],
        }
    }

    pub fn from_polys(polys: Vec<Poly<D>>) -> Self {
        PolyVec { polys }
    }

    pub fn len(&self) -> usize {
        self.polys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.polys.is_empty()
    }

    /// <a, b> = sum_i a_i * b_i in R_q.
    pub fn inner_product(&self, other: &Self) -> Poly<D> {
        assert_eq!(self.len(), other.len(), "inner_product: length mismatch");
        self.polys
            .iter()
            .zip(other.polys.iter())
            .fold(Poly::zero(), |acc, (a, b)| &acc + &(a * b))
    }

    /// Squared l2 norm of the full centered coefficient vector.
    pub fn l2_norm_squared(&self) -> u128 {
        self.polys.iter().map(|p| p.l2_norm_squared()).sum()
    }

    /// Infinity norm over all centered coefficients.
    pub fn linf_norm(&self) -> u64 {
        self.polys.iter().map(|p| p.linf_norm()).max().unwrap_or(0)
    }

    /// Coefficient-wise balanced base-b decomposition; returns `parts`
    /// vectors so that self = sum_i out[i] * b^i.
    pub fn decompose(&self, b: u64, parts: usize) -> Vec<Self> {
        let per_poly: Vec<Vec<Poly<D>>> =
            self.polys.iter().map(|p| p.decompose(b, parts)).collect();
        (0..parts)
            .map(|i| PolyVec {
                polys: per_poly.iter().map(|d| d[i].clone()).collect(),
            })
            .collect()
    }

    /// Like `decompose`, but the topmost part absorbs the remainder
    /// (see `Poly::decompose_topmost`).
    pub fn decompose_topmost(&self, b: u64, parts: usize) -> Vec<Self> {
        let per_poly: Vec<Vec<Poly<D>>> = self
            .polys
            .iter()
            .map(|p| p.decompose_topmost(b, parts))
            .collect();
        (0..parts)
            .map(|i| PolyVec {
                polys: per_poly.iter().map(|d| d[i].clone()).collect(),
            })
            .collect()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.len() * 4 * D);
        for p in &self.polys {
            out.extend_from_slice(&p.to_bytes());
        }
        out
    }
}

impl<const D: usize> Add for &PolyVec<D> {
    type Output = PolyVec<D>;
    fn add(self, rhs: &PolyVec<D>) -> PolyVec<D> {
        assert_eq!(self.len(), rhs.len());
        PolyVec {
            polys: self
                .polys
                .iter()
                .zip(rhs.polys.iter())
                .map(|(a, b)| a + b)
                .collect(),
        }
    }
}

impl<const D: usize> Sub for &PolyVec<D> {
    type Output = PolyVec<D>;
    fn sub(self, rhs: &PolyVec<D>) -> PolyVec<D> {
        assert_eq!(self.len(), rhs.len());
        PolyVec {
            polys: self
                .polys
                .iter()
                .zip(rhs.polys.iter())
                .map(|(a, b)| a - b)
                .collect(),
        }
    }
}

/// Scale by a ring element.
impl<const D: usize> Mul<&Poly<D>> for &PolyVec<D> {
    type Output = PolyVec<D>;
    fn mul(self, rhs: &Poly<D>) -> PolyVec<D> {
        PolyVec {
            polys: self.polys.iter().map(|a| a * rhs).collect(),
        }
    }
}

/// Scale by a scalar.
impl<const D: usize> Mul<Zq> for &PolyVec<D> {
    type Output = PolyVec<D>;
    fn mul(self, rhs: Zq) -> PolyVec<D> {
        PolyVec {
            polys: self.polys.iter().map(|a| a * rhs).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    fn random_vec(n: usize, rng: &mut impl Rng) -> PolyVec<64> {
        PolyVec::from_polys(
            (0..n)
                .map(|_| {
                    let mut p = Poly::<64>::zero();
                    for c in p.coeffs.iter_mut() {
                        *c = Zq::new(rng.gen_range(0..Zq::Q));
                    }
                    p
                })
                .collect(),
        )
    }

    #[test]
    fn inner_product_bilinear() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let a = random_vec(4, &mut rng);
        let b = random_vec(4, &mut rng);
        let c = random_vec(4, &mut rng);
        let left = (&a + &b).inner_product(&c);
        let right = &a.inner_product(&c) + &b.inner_product(&c);
        assert_eq!(left, right);
        // symmetry
        assert_eq!(a.inner_product(&b), b.inner_product(&a));
    }

    #[test]
    fn decompose_recompose_vector() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(12);
        let a = random_vec(3, &mut rng);
        let b = 1u64 << 16;
        let parts = a.decompose(b, 3);
        let mut acc = PolyVec::zero(3);
        let mut pow = Zq::new(1);
        for part in &parts {
            acc = &acc + &(part * pow);
            pow = pow * Zq::new(b as usize);
        }
        assert_eq!(acc, a);
    }
}
