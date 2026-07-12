use std::cmp::PartialEq;
use std::fmt::Display;
use std::iter::Sum;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::Rem;
use std::ops::Sub;

// The arithmetic below relies on usize being 64 bits wide so that products of
// two values < q < 2^32 never overflow ((q-1)^2 < 2^64).
const _: () = assert!(usize::BITS >= 64, "algebra::zq requires a 64-bit target");

// Let q be a modulus, and let Zq be the ring of integers mod q
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Zq {
    pub value: usize,
}

impl Zq {
    // todo: use symmetric from -Q/2 to Q/2
    // q = 2^32 - 99 = 4_294_967_197, the prime modulus used by the LaBRADOR reference
    // implementation (lattice-dogs/labrador, data32.c: LOGQ=32, QOFF=99).
    pub const Q: usize = (1usize << 32) - 99;

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

    /// Modular exponentiation by square-and-multiply.
    ///
    /// (The previous implementation computed `self.value.pow(e)` in machine
    /// integers, which overflows for virtually every non-trivial input.)
    pub fn pow(self, mut exp: usize) -> Self {
        let mut base = self.value as u128;
        let q = Self::Q as u128;
        let mut acc: u128 = 1;
        base %= q;
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc * base % q;
            }
            base = base * base % q;
            exp >>= 1;
        }
        Zq {
            value: acc as usize,
        }
    }

    /// Multiplicative inverse. q is prime, so a^(q-2) = a^{-1} (Fermat).
    ///
    /// Panics on zero.
    pub fn inv(self) -> Self {
        assert!(self.value != 0, "Zq::inv: zero has no inverse");
        self.pow(Self::Q - 2)
    }

    /// Centered (symmetric) representative in (-q/2, q/2].
    pub fn centered(self) -> i64 {
        let v = self.value as i64;
        let q = Self::Q as i64;
        if v > q / 2 {
            v - q
        } else {
            v
        }
    }

    /// Build an element from a (possibly negative) integer.
    pub fn from_i64(v: i64) -> Self {
        let q = Self::Q as i64;
        let mut r = v % q;
        if r < 0 {
            r += q;
        }
        Zq { value: r as usize }
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

impl std::ops::Neg for Zq {
    type Output = Zq;

    /// Additive inverse.
    fn neg(self) -> Zq {
        if self.value == 0 {
            self
        } else {
            Zq {
                value: Self::Q - self.value,
            }
        }
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
        let value = Zq::Q + 1;
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
    fn test_zq_pow_large_exponent_no_overflow() {
        // (q-1)^2 mod q = 1, and generally pow must not overflow machine ints.
        let a = Zq::new(Zq::Q - 1); // -1 mod q
        assert_eq!(a.pow(2).value(), 1);
        assert_eq!(a.pow(3), a);
        // Fermat: a^(q-1) = 1 for a != 0
        let b = Zq::new(123456789);
        assert_eq!(b.pow(Zq::Q - 1).value(), 1);
    }

    #[test]
    fn test_zq_inv() {
        for v in [1usize, 2, 99, 1 << 31, Zq::Q - 1] {
            let a = Zq::new(v);
            assert_eq!((a * a.inv()).value(), 1, "inv failed for {v}");
        }
    }

    #[test]
    fn test_zq_centered_roundtrip() {
        for v in [0i64, 1, -1, (Zq::Q as i64 - 1) / 2, -(Zq::Q as i64) / 2] {
            assert_eq!(Zq::from_i64(v).centered(), v, "centered failed for {v}");
        }
        assert_eq!(Zq::new(Zq::Q - 1).centered(), -1);
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
