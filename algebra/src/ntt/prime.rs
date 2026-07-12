//! NTT-friendly primes with self-verifying construction.
//!
//! Instead of trusting hard-coded twiddle constants (a classic source of
//! silent breakage in the C world — lazer generates them with a SageMath
//! script), primes are checked with a deterministic Miller–Rabin test and
//! roots of unity are derived and order-verified at runtime. Table
//! construction cost is microseconds and paid once per process.

/// A prime p = m * 2^a + 1 suitable for radix-2 NTTs up to size 2^a.
#[derive(Debug, Clone, Copy)]
pub struct NttPrime {
    pub p: u64,
    /// Largest a with 2^a | p - 1.
    pub two_adicity: u32,
    /// An element of multiplicative order exactly 2^two_adicity.
    root: u64,
}

impl NttPrime {
    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }

    pub fn pow(&self, mut base: u64, mut exp: u64) -> u64 {
        let mut acc: u64 = 1;
        base %= self.p;
        while exp > 0 {
            if exp & 1 == 1 {
                acc = self.mul(acc, base);
            }
            base = self.mul(base, base);
            exp >>= 1;
        }
        acc
    }

    /// A primitive n-th root of unity mod p (n must be a power of two
    /// dividing 2^two_adicity).
    pub fn root_of_order(&self, n: u64) -> u64 {
        assert!(n.is_power_of_two());
        let k = n.trailing_zeros();
        assert!(
            k <= self.two_adicity,
            "no root of order 2^{k} mod {} (2-adicity {})",
            self.p,
            self.two_adicity
        );
        let w = self.pow(self.root, 1u64 << (self.two_adicity - k));
        debug_assert_eq!(self.pow(w, n), 1);
        debug_assert_ne!(self.pow(w, n / 2), 1);
        w
    }

    /// Construct from a claimed (prime, two_adicity) pair, verifying both
    /// claims and deriving an order-2^a root.
    pub fn from_verified(p: u64, two_adicity: u32) -> Self {
        assert!(is_prime_u64(p), "{p} is not prime");
        assert_eq!(
            (p - 1).trailing_zeros(),
            two_adicity,
            "wrong 2-adicity for {p}"
        );
        let root = find_root_2adic(p, two_adicity);
        NttPrime {
            p,
            two_adicity,
            root,
        }
    }

    /// Deterministic downward search from 2^62 over p = m * 2^20 + 1.
    /// Returns the largest such prime below 2^62. Yields the same prime on
    /// every run/platform — a named constant found by construction.
    pub fn search_below_2_62() -> Self {
        let mut m = (1u64 << 42) - 1; // (m << 20) + 1 < 2^62
        loop {
            let p = (m << 20) + 1;
            if is_prime_u64(p) {
                // two_adicity may exceed 20 if m is even; recompute exactly.
                let a = (p - 1).trailing_zeros();
                return NttPrime::from_verified(p, a);
            }
            m -= 1;
        }
    }
}

/// Find an element of multiplicative order exactly 2^a mod p, where
/// 2^a || p - 1. For any x, y = x^((p-1)/2^a) has order dividing 2^a;
/// the order is exactly 2^a iff y^(2^(a-1)) = -1. No factorization of
/// p - 1 is needed.
fn find_root_2adic(p: u64, a: u32) -> u64 {
    let pr = NttPrime {
        p,
        two_adicity: a,
        root: 0,
    };
    let odd = (p - 1) >> a;
    for x in 2u64.. {
        let y = pr.pow(x, odd);
        if pr.pow(y, 1u64 << (a - 1)) == p - 1 {
            return y;
        }
    }
    unreachable!("half of all elements are non-residues");
}

/// Deterministic Miller–Rabin, complete for all u64 with this base set.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &sp in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == sp {
            return true;
        }
        if n.is_multiple_of(sp) {
            return false;
        }
    }
    let d = n - 1;
    let s = d.trailing_zeros();
    let d = d >> s;
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut b: u64, mut e: u64| {
        let mut acc = 1u64;
        b %= n;
        while e > 0 {
            if e & 1 == 1 {
                acc = mulmod(acc, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        acc
    };
    'witness: for &w in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(w, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn miller_rabin_known_values() {
        assert!(is_prime_u64(2));
        assert!(is_prime_u64(0xFFFF_FFFF_0000_0001)); // Goldilocks
        assert!(is_prime_u64((1u64 << 32) - 99)); // the LaBRADOR q
        assert!(!is_prime_u64(1u64 << 32));
        assert!(!is_prime_u64(3215031751)); // strong pseudoprime to bases 2,3,5,7
    }

    #[test]
    fn goldilocks_roots() {
        let p = NttPrime::from_verified(0xFFFF_FFFF_0000_0001, 32);
        let w = p.root_of_order(128);
        assert_eq!(p.pow(w, 128), 1);
        assert_ne!(p.pow(w, 64), 1);
        assert_eq!(p.pow(w, 64), p.p - 1); // psi^d = -1: negacyclic wrap
    }

    #[test]
    fn searched_prime_is_stable() {
        let a = NttPrime::search_below_2_62();
        let b = NttPrime::search_below_2_62();
        assert_eq!(a.p, b.p);
        assert!(a.two_adicity >= 20);
        assert!(is_prime_u64(a.p));
    }
}
