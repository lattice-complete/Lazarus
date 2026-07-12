//! Matrices over R_q (port of lazer's `polymat.c`).

use super::poly::Poly;
use super::polyvec::PolyVec;

/// Row-major dense matrix over R_q.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyMat<const D: usize> {
    pub nrows: usize,
    pub ncols: usize,
    entries: Vec<Poly<D>>,
}

impl<const D: usize> PolyMat<D> {
    pub fn zero(nrows: usize, ncols: usize) -> Self {
        PolyMat {
            nrows,
            ncols,
            entries: vec![Poly::zero(); nrows * ncols],
        }
    }

    pub fn from_rows(rows: Vec<Vec<Poly<D>>>) -> Self {
        let nrows = rows.len();
        let ncols = rows.first().map_or(0, |r| r.len());
        assert!(rows.iter().all(|r| r.len() == ncols), "ragged matrix");
        PolyMat {
            nrows,
            ncols,
            entries: rows.into_iter().flatten().collect(),
        }
    }

    pub fn get(&self, i: usize, j: usize) -> &Poly<D> {
        &self.entries[i * self.ncols + j]
    }

    pub fn set(&mut self, i: usize, j: usize, v: Poly<D>) {
        self.entries[i * self.ncols + j] = v;
    }

    pub fn row(&self, i: usize) -> PolyVec<D> {
        PolyVec::from_polys(self.entries[i * self.ncols..(i + 1) * self.ncols].to_vec())
    }

    /// Matrix-vector product A * v.
    pub fn mul_vec(&self, v: &PolyVec<D>) -> PolyVec<D> {
        assert_eq!(self.ncols, v.len(), "mul_vec: dimension mismatch");
        PolyVec::from_polys(
            (0..self.nrows)
                .map(|i| self.row(i).inner_product(v))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zq::Zq;

    #[test]
    fn identity_action() {
        let mut eye = PolyMat::<64>::zero(3, 3);
        for i in 0..3 {
            eye.set(i, i, Poly::one());
        }
        let v = PolyVec::from_polys(vec![
            Poly::monomial(1),
            Poly::monomial(2) * Poly::from_slice(&[Zq::new(5)]),
            Poly::one(),
        ]);
        assert_eq!(eye.mul_vec(&v), v);
    }

    #[test]
    fn linearity() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(21);
        let mut a = PolyMat::<64>::zero(2, 3);
        for i in 0..2 {
            for j in 0..3 {
                let mut p = Poly::zero();
                for c in p.coeffs.iter_mut() {
                    *c = Zq::new(rng.gen_range(0..Zq::Q));
                }
                a.set(i, j, p);
            }
        }
        let mut v = PolyVec::<64>::zero(3);
        let mut w = PolyVec::<64>::zero(3);
        for k in 0..3 {
            let mut p = Poly::zero();
            let mut q = Poly::zero();
            for (c, e) in p.coeffs.iter_mut().zip(q.coeffs.iter_mut()) {
                *c = Zq::new(rng.gen_range(0..Zq::Q));
                *e = Zq::new(rng.gen_range(0..Zq::Q));
            }
            v.polys[k] = p;
            w.polys[k] = q;
        }
        let left = a.mul_vec(&(&v + &w));
        let right = &a.mul_vec(&v) + &a.mul_vec(&w);
        assert_eq!(left, right);
    }
}
