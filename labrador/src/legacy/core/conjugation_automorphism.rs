use algebra::{polynomial_ring::PolynomialRing, zq::Zq};

// Conjugation Automorphism σ_{-1}
// for polynomial ring a = 1+2x+3x^2, since x^64 = -1, apply this method to a, will get 1+2*(Zq.modulus()-1) * x^(64-1) +3*(Zq.modulus()-1) * x^(64-2)
//todo: Aut(Rq) ∼= Z×2d what is this???
pub fn conjugation_automorphism(poly: &PolynomialRing) -> PolynomialRing {
    let modulus_minus_one = Zq::from(Zq::modulus() - 1);
    let transformed_coeffs: Vec<Zq> = (0..PolynomialRing::DEGREE_BOUND)
        .map(|i| {
            if i < poly.coefficients.len() {
                if i == 0 {
                    poly.coefficients[i]
                } else {
                    poly.coefficients[i] * modulus_minus_one
                }
            } else {
                Zq::from(0)
            }
        })
        .collect();
    // reverse the coefficients except constant term
    let reversed_coefficients = transformed_coeffs
        .iter()
        .take(1)
        .cloned()
        .chain(transformed_coeffs.iter().skip(1).rev().cloned())
        .collect::<Vec<Zq>>();
    PolynomialRing {
        coefficients: reversed_coefficients,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::utils::inner_product_zq_vector;
    #[test]
    fn test_conjugation_automorphism() {
        // Create example PolynomialRings a and b
        let a = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let b = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };

        // Compute <a, b>
        let inner_ab = inner_product_zq_vector(&a.coefficients, &b.coefficients);
        assert_eq!(inner_ab.value(), 32);
        // Compute σ_{-1}(a)
        let sigma_inv_a = conjugation_automorphism(&a);
        // Compute <σ_{-1}(a), b>
        let inner_sigma_inv_a_b = &sigma_inv_a * &b;
        // = -12x^{65}-28x^{64}-23x^{63}-12x^{62}+6x^{2}+5x+4
        // = 23x^{63}-12x^{62}+6x^{2}+17x+32 (since x^64 = -1)
        assert_eq!(inner_sigma_inv_a_b.coefficients.len(), 64);
        assert_eq!(inner_sigma_inv_a_b.coefficients[0], Zq::from(32));
        assert_eq!(inner_sigma_inv_a_b.coefficients[1], Zq::from(17));
        assert_eq!(inner_sigma_inv_a_b.coefficients[2], Zq::from(6));
        assert_eq!(inner_sigma_inv_a_b.coefficients[62], Zq::from(Zq::Q - 12));
        assert_eq!(inner_sigma_inv_a_b.coefficients[63], Zq::from(Zq::Q - 23));
        for i in 3..62 {
            assert_eq!(inner_sigma_inv_a_b.coefficients[i], Zq::from(0));
        }

        // Get the constant term of <σ_{-1}(a), b>
        let ct_inner_sigma_inv_a_b = inner_sigma_inv_a_b.coefficients[0];

        // Assert that <a, b> == ct <σ_{-1}(a), b>
        assert_eq!(
            inner_ab, ct_inner_sigma_inv_a_b,
            "<a, b> should equal the constant term of <σ-1(a), b>"
        );
    }
}
