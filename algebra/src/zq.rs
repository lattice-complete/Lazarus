use std::cmp::PartialEq;
use std::fmt::Display;
use std::iter::Sum;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::Rem;
use std::ops::Sub;

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

    pub fn value(self) -> usize {
        self.value
    }

    pub fn pow(self, other: usize) -> Self {
        Zq::new(self.value.pow(other as u32))
    }

    // bypass clippy for now as we'll implement Zq in isize/i64 (no div_ceil function supported) later
    #[allow(clippy::manual_div_ceil)]
    pub fn div_ceil(self, other: Zq) -> Zq {
        Zq::new((self.value + other.value - 1) / other.value)
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
    fn test_zq_ceil_div() {
        let mut a;
        let mut b;
        a = Zq::new(10);
        b = Zq::new(3);
        assert_eq!(a.div_ceil(b).value(), 4);

        a = Zq::new(9);
        b = Zq::new(3);
        assert_eq!(a.div_ceil(b).value(), 3);
    }
}
