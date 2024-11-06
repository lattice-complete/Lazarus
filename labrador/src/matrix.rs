use rand::Rng;

pub fn generate_symmetric_tridiagonal_matrix(n: usize) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();

    // Initialize a matrix with all zeros
    let mut matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        // Set the main diagonal elements
        matrix[i][i] = rng.gen_range(1.0..10.0); // Random non-zero value

        // Set the elements just above and below the main diagonal
        if i > 0 {
            let value = rng.gen_range(1.0..10.0);
            matrix[i][i - 1] = value; // Below the main diagonal
            matrix[i - 1][i] = value; // Symmetric element above the main diagonal
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_symmatric_tridiagonal_matrix() {
        let n = 3;
        let matrix = generate_symmetric_tridiagonal_matrix(n);
        println!("{:?}", matrix);
    }
}
