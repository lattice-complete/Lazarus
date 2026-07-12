use algebra::{polynomial_ring::PolynomialRing, zq::Zq};
// convert number to basis
// 42 = 0 * 2^7 + 1 * 2^6 + 0 * 2^5 + 1 * 2^4 + 0 * 2^3 + 1 * 2^2 + 0 * 2^1 + 0 * 2^0
// first digit: 42 / 2 = 21, result_i = 0
// second digit: 21 / 2 = 10, result_i = 1
// third digit: 10 / 2 = 5, result_i = 0
// forth digit: 5 / 2 = 2, result_i = 1
// fifth digit: 2 / 2 = 1, result_i = 0
// sixth digit: 1 / 2 = 0, result_i = 1

pub fn num_to_basis(num: Zq, basis: Zq, digits: Zq) -> Vec<Zq> {
    let mut result = Vec::new();
    let mut remainder = num;

    let zero = Zq::from(0);
    let base = basis;

    for _ in 0..digits.value() {
        let digit = remainder % base;
        result.push(digit);
        remainder = Zq::from(remainder.value() / base.value());
    }

    while result.len() < digits.value() {
        // push 0 to the highest position
        result.push(zero);
    }

    result
}

// convert ring polynomial to basis
pub fn ring_polynomial_to_basis(poly: &PolynomialRing, basis: Zq, digits: Zq) -> Vec<Vec<Zq>> {
    poly.coefficients
        .iter()
        .map(|&coeff| num_to_basis(coeff, basis, digits))
        .collect()
}

pub fn poly_vec_decompose_to_basis(
    poly: &[PolynomialRing],
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<Zq>>> {
    poly.iter()
        .map(|poly_i| ring_polynomial_to_basis(poly_i, basis, digits))
        .collect::<Vec<Vec<Vec<Zq>>>>()
}

pub fn poly_matrix_decompose_to_basis(
    poly: &[Vec<PolynomialRing>],
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<Vec<Zq>>>> {
    poly.iter()
        .map(|poly_i| poly_vec_decompose_to_basis(poly_i, basis, digits))
        .collect::<Vec<Vec<Vec<Vec<Zq>>>>>()
}

// TODO(junochiu): check where best to put this function / any renaming needed
// TODO(junochiu): this function is giving error when applying changes to fix ptr_arg, debug needed
// aggregate basis form of a vector of PolynomialRing
#[allow(clippy::ptr_arg)]
pub fn aggregate_poly_vec_basis_form(
    poly_basis_form: &Vec<Vec<Vec<Zq>>>,
) -> Vec<Vec<PolynomialRing>> {
    poly_basis_form
        .iter()
        .map(|poly_i_j_basis_form| {
            let num_loop_needed = poly_i_j_basis_form.first().map_or(0, |v| v.len());
            (0..num_loop_needed)
                .map(|k| {
                    let coefficients = poly_i_j_basis_form
                        .iter()
                        .map(|basis_part| basis_part.get(k).cloned().unwrap_or(Zq::from(0)))
                        .collect();
                    PolynomialRing { coefficients }
                })
                .collect()
        })
        .collect()
}

// TODO(junochiu): check where best to put this function / any renaming needed
pub fn poly_matrix_decompose_and_aggregate(
    poly: &[Vec<PolynomialRing>],
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<PolynomialRing>>> {
    // Decompose h_ij into basis t_1 parts
    let poly_basis_form = poly_matrix_decompose_to_basis(poly, basis, digits);

    // Pick elements at each position across all inner vectors and aggregate them
    poly_basis_form
        .iter()
        .map(aggregate_poly_vec_basis_form)
        .collect()
}

// TODO(junochiu): check where best to put this function / any renaming needed
pub fn poly_vec_decompose_and_aggregate(
    poly: &[PolynomialRing],
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<PolynomialRing>> {
    // Decompose each PolynomialRing into basis parts
    let poly_basis_form = poly_vec_decompose_to_basis(poly, basis, digits);
    aggregate_poly_vec_basis_form(&poly_basis_form)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_num_to_basis() {
        let num = Zq::from(42);
        let basis = Zq::from(2);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(3);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(1),
                Zq::from(0),
                Zq::from(2),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(6);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(4),
                Zq::from(4),
                Zq::from(2),
                Zq::from(0),
                Zq::from(0),
                Zq::from(0)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(10);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(0)
            ]
        );
    }

    #[test]
    fn test_basis_to_num_vector() {
        let basis = Zq::from(10);
        let digits = Zq::from(3);
        let vec1 = [8, 46, 61, 71, 33, 33, 18];
        let vec2 = [20, 54, 94, 93, 70, 33, 14];
        let vec3 = [24, 40, 100, 85, 121, 57, 56];
        let vec4 = [14, 37, 91, 118, 159, 109, 72];
        for vec in [vec1, vec2, vec3, vec4] {
            let mut temp_vec: Vec<Vec<Zq>> = Vec::new();
            for i in vec {
                let num = num_to_basis(Zq::from(i), basis, digits);
                temp_vec.push(num);
            }
        }
    }

    #[test]
    fn test_ring_polynomial_to_basis() {
        let poly = PolynomialRing {
            coefficients: vec![Zq::from(42), Zq::from(100), Zq::from(100)],
        };
        let basis = Zq::from(2);
        let digits = Zq::from(8);
        let expected_result = vec![
            vec![
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
            ],
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(1),
                Zq::from(0),
            ],
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(1),
                Zq::from(0),
            ],
        ];
        let result = ring_polynomial_to_basis(&poly, basis, digits);
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_decompose_poly_to_basis_form() {
        // Arrange: Create sample input polynomial rings
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(123), Zq::from(456), Zq::from(789)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(12), Zq::from(45), Zq::from(78)],
        };
        let poly_input = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly1.clone(), poly2.clone()],
        ];
        let basis = Zq::from(10);
        let digits = Zq::from(3);

        // Act: Call the function to decompose the polynomial
        let result = poly_matrix_decompose_and_aggregate(&poly_input, basis, digits);

        let expected_row1 = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(6), Zq::from(9)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(0), Zq::from(0), Zq::from(0)],
                },
            ],
        ];

        let expected = vec![expected_row1.clone(), expected_row1.clone()];
        assert_eq!(
            result, expected,
            "The decomposition did not match the expected output."
        );
    }

    #[test]
    fn test_decompose_poly_ring_vector_to_basis_form() {
        // Arrange: Create sample input vector of PolynomialRing
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(123), Zq::from(456), Zq::from(789)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(12), Zq::from(45), Zq::from(78)],
        };
        let poly_vec = vec![poly1.clone(), poly2.clone()];
        let basis = Zq::from(10);
        let digits = Zq::from(3);

        // Act: Call the decomposition function
        let result = poly_vec_decompose_and_aggregate(&poly_vec, basis, digits);

        // Define expected output
        let expected = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(6), Zq::from(9)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(0), Zq::from(0), Zq::from(0)],
                },
            ],
        ];

        // Assert
        assert_eq!(
            result, expected,
            "The decomposition did not match the expected output."
        );
    }
}
