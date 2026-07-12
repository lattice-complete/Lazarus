//! Witness: a tuple of R_q vectors with tracked norms (C: witness).

use crate::Rq;
use algebra::ring::PolyVec;

#[derive(Debug, Clone)]
pub struct Witness {
    pub parts: Vec<PolyVec<64>>,
}

impl Witness {
    pub fn new(parts: Vec<PolyVec<64>>) -> Self {
        Witness { parts }
    }

    /// Squared l2 norm of each part (centered representatives).
    pub fn norms(&self) -> Vec<u128> {
        self.parts.iter().map(|p| p.l2_norm_squared()).collect()
    }

    pub fn total_normsq(&self) -> u128 {
        self.norms().iter().sum()
    }

    /// Flatten the witness into the plan's r*n chunk layout, zero-padding
    /// group tails (C: the sx assembly in commit()).
    pub fn assemble(&self, plan: &crate::params::Plan) -> Vec<Rq> {
        let mut sx = vec![Rq::zero(); plan.r * plan.n];
        for (i, part) in self.parts.iter().enumerate() {
            let base = plan.part_base[i];
            for (j, p) in part.polys.iter().enumerate() {
                sx[base + j] = p.clone();
            }
        }
        sx
    }
}
