use crate::zq::Zq;
use std::cmp::PartialEq;
use std::ops::Add;
use std::ops::Mul;

#[derive(Debug, Clone, Eq)]

pub struct PolynomialRing {
    pub coefficients: Vec<Zq>,
}

impl PolynomialRing {
    pub const DEGREE_BOUND: usize = 64;

    pub fn new(coefficients: Vec<Zq>) -> Self {
        // Ensure coefficients are in Zq and degree is less than 64
        assert!(
            coefficients.len() <= Self::DEGREE_BOUND,
            "Polynomial degree must be less than 64"
        );
        PolynomialRing {
            coefficients: coefficients.clone(),
        }
    }

    // todo: use NTT to speed up
    // Multiply two polynomials in R = Zq[X]/(X^64 + 1)
    fn multiply_by_polynomial_ring(&self, other: &PolynomialRing) -> PolynomialRing {
        // Initialize a vector to hold the intermediate multiplication result
        let mut result_coefficients =
            vec![Zq::new(0); self.coefficients.len() + other.coefficients.len() - 1];
        for (i, &coeff1) in self.coefficients.iter().enumerate() {
            for (j, &coeff2) in other.coefficients.iter().enumerate() {
                result_coefficients[i + j] += coeff1 * coeff2;
            }
        }

        // Reduce modulo X^64 + 1
        if result_coefficients.len() > Self::DEGREE_BOUND {
            let modulus_minus_one = Zq::from(Zq::modulus() - 1);
            let (front, back) = result_coefficients.split_at_mut(Self::DEGREE_BOUND);
            for (i, &overflow) in back.iter().enumerate() {
                front[i] += overflow * modulus_minus_one;
            }
            result_coefficients.truncate(Self::DEGREE_BOUND);
        }
        PolynomialRing {
            coefficients: result_coefficients,
        }
    }

    fn add_polynomial_ring(&self, other: &PolynomialRing) -> PolynomialRing {
        let len_a = self.coefficients.len();
        let len_b = other.coefficients.len();
        let max_len = std::cmp::max(len_a, len_b);
        let mut result_coefficients = Vec::with_capacity(max_len);

        result_coefficients.extend(
            self.coefficients
                .iter()
                .zip(other.coefficients.iter())
                .map(|(&a, &b)| a + b),
        );

        match len_a.cmp(&len_b) {
            std::cmp::Ordering::Greater => {
                result_coefficients.extend_from_slice(&self.coefficients[len_b..]);
            }
            std::cmp::Ordering::Less => {
                result_coefficients.extend_from_slice(&other.coefficients[len_a..]);
            }
            std::cmp::Ordering::Equal => {
                // Do nothing
            }
        }
        PolynomialRing {
            coefficients: result_coefficients,
        }
    }
}

impl Mul for PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: PolynomialRing) -> PolynomialRing {
        self.multiply_by_polynomial_ring(&other)
    }
}

impl Mul<&PolynomialRing> for PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: &PolynomialRing) -> PolynomialRing {
        self.multiply_by_polynomial_ring(other)
    }
}

impl Mul<PolynomialRing> for &PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: PolynomialRing) -> PolynomialRing {
        self.multiply_by_polynomial_ring(&other)
    }
}

impl Mul<&PolynomialRing> for &PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: &PolynomialRing) -> PolynomialRing {
        self.multiply_by_polynomial_ring(other)
    }
}

impl Mul<Zq> for PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: Zq) -> PolynomialRing {
        let new_coefficients = self.coefficients.iter().map(|c| *c * other).collect();
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl Mul<Zq> for &PolynomialRing {
    type Output = PolynomialRing;

    fn mul(self, other: Zq) -> PolynomialRing {
        let new_coefficients = self.coefficients.iter().map(|c| *c * other).collect();
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl Add for PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: PolynomialRing) -> PolynomialRing {
        self.add_polynomial_ring(&other)
    }
}

impl Add<&PolynomialRing> for PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: &PolynomialRing) -> PolynomialRing {
        self.add_polynomial_ring(other)
    }
}

impl Add<PolynomialRing> for &PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: PolynomialRing) -> PolynomialRing {
        self.add_polynomial_ring(&other)
    }
}

impl Add<&PolynomialRing> for &PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: &PolynomialRing) -> PolynomialRing {
        self.add_polynomial_ring(other)
    }
}

impl Add<Zq> for PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: Zq) -> PolynomialRing {
        let mut new_coefficients = self.coefficients.clone();
        if let Some(first) = new_coefficients.get_mut(0) {
            *first += other;
        } else {
            new_coefficients.push(other);
        }
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl Add<&Zq> for PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: &Zq) -> PolynomialRing {
        let mut new_coefficients = self.coefficients.clone();
        if let Some(first) = new_coefficients.get_mut(0) {
            *first += *other;
        } else {
            new_coefficients.push(*other);
        }
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl Add<Zq> for &PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: Zq) -> PolynomialRing {
        let mut new_coefficients = self.coefficients.clone();
        if let Some(first) = new_coefficients.get_mut(0) {
            *first += other;
        } else {
            new_coefficients.push(other);
        }
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl Add<&Zq> for &PolynomialRing {
    type Output = PolynomialRing;

    fn add(self, other: &Zq) -> PolynomialRing {
        let mut new_coefficients = self.coefficients.clone();
        if let Some(first) = new_coefficients.get_mut(0) {
            *first += *other;
        } else {
            new_coefficients.push(*other);
        }
        PolynomialRing {
            coefficients: new_coefficients,
        }
    }
}

impl PartialEq for PolynomialRing {
    fn eq(&self, other: &Self) -> bool {
        // Compare coefficients, ignoring trailing zeros
        let self_coeffs = self
            .coefficients
            .iter()
            .rev()
            .skip_while(|&&x| x == Zq::from(0))
            .collect::<Vec<_>>();
        let other_coeffs = other
            .coefficients
            .iter()
            .rev()
            .skip_while(|&&x| x == Zq::from(0))
            .collect::<Vec<_>>();

        self_coeffs == other_coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply_by_polynomial_ring() {
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(3), Zq::from(4)],
        };
        let result = poly1 * poly2;
        assert_eq!(
            result.coefficients,
            vec![Zq::from(3), Zq::from(10), Zq::from(8)]
        ); // 1*3, 1*4 + 2*3, 2*4
    }

    #[test]
    fn test_polynomial_ring_mul_overflow() {
        // Create two polynomials that will cause overflow when multiplied
        // For example, (X^63 + 1) * (X^63 + x) = X^126 + X^64 + X^63 + X
        // Modulo X^64 + 1, X^64 = -1, so X^126 = X^(2*64 -2) = X^-2 = X^62, X*X^63 = -1
        // Thus, X^126 + X^64 + X^63 + X mod X^64+1 = (-1)*X^62 + (-1) + X + X^63 = - 1 + X - X^62 + X^63

        // Initialize poly1 as X^63 + 1
        let mut poly1_coeffs = vec![Zq::from(0); 64];
        poly1_coeffs[0] = Zq::from(1); // Constant term
        poly1_coeffs[63] = Zq::from(1); // X^63 term
        let poly1 = PolynomialRing {
            coefficients: poly1_coeffs,
        };

        // Initialize poly1 as X^63 + X
        let mut poly2_coeffs = vec![Zq::from(0); 64];
        poly2_coeffs[1] = Zq::from(1); // X term
        poly2_coeffs[63] = Zq::from(1); // X^63 term
        let poly2 = PolynomialRing {
            coefficients: poly2_coeffs,
        };

        // Multiply poly1 by poly2
        let product = poly1.clone() * poly2.clone();

        // Expected coefficients after reduction modulo X^64 + 1:
        // coefficients[0] = 1
        // coefficients[62] = Zq::modulus() - 1  (since -1 mod q)
        // coefficients[63] = 2
        // All other coefficients should be 0
        let mut expected_coeffs = vec![Zq::from(0); 64];
        expected_coeffs[0] = Zq::from(Zq::modulus() - 1); // Constant term
        expected_coeffs[1] = Zq::from(1); // X term
        expected_coeffs[62] = Zq::from(Zq::modulus() - 1); // X^62 term
        expected_coeffs[63] = Zq::from(1); // X^63 term

        // Assert that the product has the correct degree bound
        assert_eq!(
            product.coefficients.len(),
            64,
            "Product should be truncated to DEGREE_BOUND"
        );

        // Assert that the coefficients match the expected values
        assert_eq!(
            product.coefficients, expected_coeffs,
            "Overflow handling in multiplication is incorrect"
        );
    }

    /// Regression test for the LaBRADOR ring-bug write-up (against commit d752488).
    ///
    /// With the buggy modulus q = 2^32, choosing u with every coefficient = 2^31
    /// makes u·α = 0 in R = Zq[X]/(X^64+1) with probability 1/2 over uniform α
    /// (2^31 is a zero divisor), reducing protocol soundness to 1 bit. With q
    /// prime, Zq is a field and α·u = 0 forces α = 0 — vanishingly rare. This
    /// guards against any future regression that swaps q back to a power of two.
    #[test]
    fn ring_bug_regression_zero_divisor() {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let u = PolynomialRing {
            coefficients: vec![Zq::from(1usize << 31); PolynomialRing::DEGREE_BOUND],
        };

        const TRIALS: usize = 1_000;
        let mut zero_products = 0;
        for _ in 0..TRIALS {
            let alpha = PolynomialRing {
                coefficients: (0..PolynomialRing::DEGREE_BOUND)
                    .map(|_| Zq::from(rng.gen_range(0..Zq::Q)))
                    .collect(),
            };
            if (&alpha * &u).coefficients.iter().all(|c| c.value() == 0) {
                zero_products += 1;
            }
        }

        assert_eq!(
            zero_products, 0,
            "ring bug regression: u = (2^31, …, 2^31) collided with random α \
             {zero_products}/{TRIALS} times — q must be prime, not a power of 2"
        );
    }

    // add test for polynomial addition and multiplication with overload
    #[test]
    fn test_polynomial_addition_and_multiplication() {
        let a = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let b = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let c = &a + &b;
        assert_eq!(c.coefficients, vec![Zq::from(5), Zq::from(7), Zq::from(9)]);
        let d = &a * &b;
        assert_eq!(
            d.coefficients,
            vec![
                Zq::from(4),
                Zq::from(13),
                Zq::from(28),
                Zq::from(27),
                Zq::from(18)
            ]
        );
    }
}
