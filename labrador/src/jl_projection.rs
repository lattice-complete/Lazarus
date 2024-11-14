use lark_algebra::{
    ConfigRingGoldilocks256, ConfigZZpX, ConfigZZpXGoldilocks256, Goldilocks, Polynomial,
    RingGoldilock256,
};
use rand::distributions::{Distribution, WeightedIndex};
use rand::Rng;
use std::ops::Neg;

struct PiDistribution;

impl Distribution<i8> for PiDistribution {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> i8 {
        let values = [-1, 0, 1];
        let weights = [1, 2, 1];
        let dist = WeightedIndex::new(&weights).unwrap();
        values[dist.sample(rng)]
    }
}

/// Projects a matrix of polynomial ring into a lower-dimensional space using a random projection.
///
/// Function Overview:
/// - Input: A matrix `D` over a polynomial ring with dimensions `(rows, original_domain)`.
/// - Output: A projected matrix `P` over the same polynomial ring with dimensions `(rows, target_domain)`.
/// - Process:
///   - Generate a random projection matrix `Π` of dimensions `(target_domain, original_domain)`.
///     - Each entry in `Π` is a polynomial ring with coefficients randomly selected from the set `{ -1, 0, 1 }` with probabilities `{ 1/4, 1/2, 1/4 }` respectively.
///   - Compute the projected matrix `P` by fot product `D` with `Π`
/// - Property:
///   - The L2 norms (Euclidean norms) of `D` and `P` are approximately preserved due to the Johnson-Lindenstrauss lemma.
///
/// Notes:
/// - This technique is useful for dimensionality reduction while maintaining the structural integrity of the data.
/// - The random projection is efficient and suitable for large-scale data processing in polynomial rings.

pub fn jl_projection(
    data: Vec<Vec<RingGoldilock256>>,
    target_dim: usize,
) -> Vec<Vec<RingGoldilock256>> {
    let original_dim = ConfigRingGoldilocks256::DIM;
    let mut projection_matrix: Vec<Vec<RingGoldilock256>> = Vec::with_capacity(target_dim);
    for _ in 0..target_dim {
        let mut row: Vec<RingGoldilock256> = Vec::with_capacity(original_dim);
        for _ in 0..original_dim {
            // generate random polynomial ring with coefficients from the set {-1, 0, 1}
            let cell = match PiDistribution.sample(&mut rand::thread_rng()) {
                -1 => Goldilocks::from(1u64).neg(),
                0 => Goldilocks::from(0u64),
                1 => Goldilocks::from(1u64),
                _ => unreachable!(),
            };
            let poly = vec![cell; ConfigZZpXGoldilocks256::DIM];
            row.push(RingGoldilock256::from_coefficients_vec_unchecked(poly));
        }
        projection_matrix.push(row);
    }

    let mut projected_data = Vec::with_capacity(data.len());
    for data_vec in &data {
        let mut projected_vector: Vec<RingGoldilock256> = Vec::with_capacity(target_dim);
        for projection_row in &projection_matrix {
            let projection_value = RingGoldilock256::dot_product(&data_vec, &projection_row);
            projected_vector.push(projection_value);
        }
        projected_data.push(projected_vector);
    }

    projected_data
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_jl_projection() {
        let mut coeffs_s1_0 = vec![
            Goldilocks::from(1u64), // Coefficient for x^0
            Goldilocks::from(2u64), // Coefficient for x^1
        ];
        coeffs_s1_0.resize(ConfigZZpXGoldilocks256::DIM, Goldilocks::from(0u64));
        let s1_0 = RingGoldilock256::from_coefficients_vec_unchecked(coeffs_s1_0);

        let mut coeffs_s1_1 = vec![
            Goldilocks::from(3u64), // Coefficient for x^0
            Goldilocks::from(4u64), // Coefficient for x^1
        ];
        coeffs_s1_1.resize(ConfigZZpXGoldilocks256::DIM, Goldilocks::from(0u64));
        let s1_1 = RingGoldilock256::from_coefficients_vec_unchecked(coeffs_s1_1);

        let s1 = vec![s1_0, s1_1];
        let target_dim = 128;

        let projection = jl_projection(vec![s1.clone(), s1.clone()], target_dim);
        assert!(
            projection.len() == 2,
            "projection row count should be the same as the original data"
        );
        assert!(
            projection[0].len() == target_dim,
            "projection column len should be the same as the target dimension"
        );
    }
}
