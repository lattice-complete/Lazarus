use algebra::{polynomial_ring::PolynomialRing, zq::Zq};

// Calculate the sum of squared norms for a single PolynomialRing instance.
pub fn poly_norm_squared(poly: &PolynomialRing) -> Zq {
    poly.coefficients
        .iter()
        .fold(Zq::new(0), |acc, coeff| acc + coeff.pow(2))
}

// Calculate the sum of squared norms for a vector of PolynomialRing instances.
pub fn poly_vec_norm_squared(polys: &[PolynomialRing]) -> Zq {
    polys
        .iter()
        .fold(Zq::new(0), |acc, poly| acc + poly_norm_squared(poly))
}

// Calculate the sum of squared norms for a matrix of PolynomialRing instances.
pub fn poly_matrix_norm_squared(poly_matrix: &[Vec<PolynomialRing>]) -> Zq {
    poly_matrix.iter().fold(Zq::new(0), |acc, vector| {
        acc + poly_vec_norm_squared(vector)
    })
}

// Calculate the sum of squared norms for a 3D vector of PolynomialRing instances.
pub fn poly_3d_norm_squared(polymat3d: &[Vec<Vec<PolynomialRing>>]) -> Zq {
    polymat3d.iter().fold(Zq::new(0), |acc, poly_matrix| {
        acc + poly_matrix_norm_squared(poly_matrix)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::vec;

    #[test]
    fn test_poly_norm_squared() {
        let poly = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let expected = Zq::from(14); // 1^2 + 2^2 + 3^2 = 14
        let result = poly_norm_squared(&poly);
        assert_eq!(
            result, expected,
            "poly_norm_squared should return the sum of squared coefficients"
        );
    }

    #[test]
    fn test_poly_matrix_norm_squared() {
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let poly_matrix = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly2.clone(), poly1.clone()],
        ];
        // poly_norm_squared(poly1) = 1 + 4 + 9 = 14
        // poly_norm_squared(poly2) = 16 + 25 + 36 = 77
        // Total sum: 14 + 77 + 77 + 14 = 182
        let expected = Zq::from(182);
        let result = poly_matrix_norm_squared(&poly_matrix);
        assert_eq!(result, expected, "poly_matrix_norm_squared should return the sum of squared norms of all polynomials in the matrix");
    }

    #[test]
    fn test_poly_vec_norm_squared() {
        // Arrange
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let vec_polys = vec![poly1.clone(), poly2.clone()];

        // Act
        let result = poly_vec_norm_squared(&vec_polys);

        // Assert
        // poly1 norm: 1^2 + 2^2 + 3^2 = 14
        // poly2 norm: 4^2 + 5^2 + 6^2 = 77
        let expected = Zq::from(14) + Zq::from(77);
        assert_eq!(
            result, expected,
            "poly_vec_norm_squared did not return the correct sum of squared norms"
        );
    }

    #[test]
    fn test_poly_3d_norm_squared() {
        // Arrange
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let poly_matrix = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly1.clone(), poly2.clone()],
        ];
        let polymat3d = vec![poly_matrix.clone(), poly_matrix.clone()];

        // Act
        let result = poly_3d_norm_squared(&polymat3d);

        // Assert
        // Each poly_matrix contains two vectors of polynomials, each vector has 2 polynomials with norms 14 and 77
        // Each matrix: 2 vectors * (14 + 77) = 2 * 91 = 182
        // Total: 2 matrices * 182 = 364
        let expected = Zq::from(364);
        assert_eq!(
            result, expected,
            "poly_3d_norm_squared did not return the correct sum of squared norms"
        );
    }
}
