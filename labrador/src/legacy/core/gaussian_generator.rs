use algebra::Zq;
use rand::Rng;

pub fn generate_gaussian_distribution(nd: Zq) -> Vec<Vec<Zq>> {
    let nd_usize: usize = nd.value();
    let modulus: usize = Zq::modulus();
    let mut rng = rand::thread_rng();
    let mut matrix = vec![vec![Zq::from(0); nd_usize]; 256]; // Initialize a 256 x nd matrix

    for row in matrix.iter_mut() {
        for cell in row.iter_mut() {
            let random_value: f32 = rng.gen(); // Generate a random float between 0 and 1
            *cell = if random_value < 0.25 {
                // todo: should we use symmetric distribution from -q/2 to q/2?
                Zq::from(modulus - 1) // 1/4 probability
            } else if random_value < 0.75 {
                Zq::from(0) // 1/2 probability
            } else {
                Zq::from(1) // 1/4 probability
            };
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_gaussian_distribution() {
        let nd = Zq::from(10);
        let matrix = generate_gaussian_distribution(nd);
        assert_eq!(matrix.len(), 256);
        assert_eq!(matrix[0].len(), nd.value());
        assert_eq!(matrix[1].len(), nd.value());
        assert!(matrix.iter().all(|row| row
            .iter()
            .all(|&val| val.value == Zq::Q - 1 || val.value == 0 || val.value == 1)));
    }
}
