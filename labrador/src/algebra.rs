use std::cmp::PartialEq;
use std::fmt::Display;
use std::iter::Sum;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Rem;
use std::ops::Sub;

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

// Let q be a modulus, and let Zq be the ring of integers mod q
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Zq {
    pub value: usize,
}

impl Zq {
    // todo: use symmetric from -Q/2 to Q/2
    pub const Q: usize = 2usize.pow(32);

    pub fn modulus() -> usize {
        Self::Q
    }
    pub fn new(value: usize) -> Self {
        Zq {
            value: value % Self::Q,
        }
    }

    pub fn value(&self) -> usize {
        self.value
    }

    pub fn pow(&self, other: usize) -> Self {
        Zq::new(self.value.pow(other as u32))
    }
}

impl PartialOrd for Zq {
    fn partial_cmp(&self, other: &Zq) -> Option<std::cmp::Ordering> {
        Some(self.value.cmp(&other.value))
    }
}

impl Display for Zq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl Rem for Zq {
    type Output = Zq;

    fn rem(self, other: Zq) -> Zq {
        Zq::new(self.value % other.value)
    }
}

impl Add for Zq {
    type Output = Zq;

    fn add(self, other: Zq) -> Zq {
        // assert!(self.value + other.value < Self::Q, "Addition result exceeds modulus");
        Zq::new(self.value + other.value)
    }
}

impl AddAssign for Zq {
    fn add_assign(&mut self, other: Zq) {
        self.value = (self.value + other.value) % Self::Q;
    }
}

impl Sub for Zq {
    type Output = Zq;

    fn sub(self, other: Zq) -> Zq {
        Zq::new((self.value + Self::Q) - other.value)
    }
}

impl Mul for Zq {
    type Output = Zq;

    fn mul(self, other: Zq) -> Zq {
        // assert!(self.value * other.value < Self::Q, "Multiplication result exceeds modulus");
        Zq::new(self.value * other.value)
    }
}

impl From<usize> for Zq {
    fn from(value: usize) -> Self {
        Zq::new(value)
    }
}

impl Sum for Zq {
    fn sum<I: Iterator<Item = Zq>>(iter: I) -> Self {
        iter.fold(Zq::new(0), |acc, x| acc + x)
    }
}

#[derive(Debug)]
pub struct RqMatrix {
    pub values: Vec<Vec<PolynomialRing>>, // matrix of PolynomialRing values
}

impl RqMatrix {
    pub fn new(kappa: Zq, size_n: Zq) -> Self {
        let size_kappa_usize: usize = kappa.value();
        let size_n_usize: usize = size_n.value();
        let mut rng = rand::thread_rng();
        let values: Vec<Vec<PolynomialRing>> = (0..size_kappa_usize)
            .map(|_| {
                (0..size_n_usize)
                    .map(|_| PolynomialRing {
                        coefficients: (0..size_n_usize)
                            .map(|_| Zq::from(2)) // we just use 2 as random number to facilitate test
                            .collect(),
                    })
                    .collect()
            })
            .collect();
        assert_eq!(
            values.len(),
            size_kappa_usize,
            "values must have the same length as size_kappa"
        );
        assert_eq!(
            values[0].len(),
            size_n_usize,
            "values[0] must have the same length as size_n"
        );
        RqMatrix { values }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zq_addition() {
        let a = Zq::new(10);
        let b = Zq::new(20);
        let result = a + b;
        assert_eq!(result.value, 30);
    }

    #[test]
    fn test_zq_subtraction() {
        let a = Zq::new(10);
        let b = Zq::new(5);
        let result = a - b;
        assert_eq!(result.value, 5);
    }

    #[test]
    fn test_zq_multiplication() {
        let a = Zq::new(6);
        let b = Zq::new(7);
        let result = a * b;
        assert_eq!(result.value, 42);
    }

    #[test]
    fn test_zq_multiplication_overflow() {
        let a = Zq::new(Zq::Q - 1);
        let b = Zq::new(2);
        let result = a * b;
        let expected = Zq::new(Zq::Q - 2); // -2
        assert_eq!(result, expected);
    }

    #[test]
    fn test_zq_overflow() {
        let a = Zq::new(Zq::Q - 1);
        let b = Zq::new(2);
        let result = a + b;
        assert_eq!(result.value, 1); // (2^32 - 1) + 2 mod 2^32 = 1
    }

    #[test]
    fn test_zq_new() {
        let value = 4294967297; // Q + 1
        let zq = Zq::new(value);
        assert_eq!(zq.value, 1);
    }

    #[test]
    fn test_zq_remainder() {
        let a = Zq::new(10);
        let b = Zq::new(3);
        let result = a % b;
        assert_eq!(result.value, 1);
    }

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
}
