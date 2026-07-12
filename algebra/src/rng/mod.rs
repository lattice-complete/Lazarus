//! Deterministic, domain-separated RNG layer (port of lazer's `rng.c`).
//!
//! lazer expands all randomness from a 32-byte seed plus a 64-bit domain
//! separator, through either SHAKE128 (absorb seed || dom_le, then squeeze)
//! or AES256-CTR (key = seed, 16-byte IV = dom_le || 0). Both backends are
//! provided here on top of the RustCrypto crates; SHAKE128 is the default
//! and is bit-compatible with any standard Keccak implementation.
//!
//! Note on AES-CTR compatibility: the C code lays the IV out as the two
//! little-endian words { dom, 0 } and increments the trailing counter word;
//! `Ctr64LE` reproduces that layout. Cross-validation against C output is
//! tracked in docs/MIGRATION.md.

use aes::cipher::{KeyIvInit, StreamCipher};
use sha3::digest::{ExtendableOutput, Update, XofReader};

/// A deterministic byte stream seeded by (seed, dom).
pub trait Xof {
    fn fill(&mut self, out: &mut [u8]);

    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.fill(&mut b);
        u32::from_le_bytes(b)
    }

    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.fill(&mut b);
        u64::from_le_bytes(b)
    }

    fn next_u128(&mut self) -> u128 {
        let mut b = [0u8; 16];
        self.fill(&mut b);
        u128::from_le_bytes(b)
    }
}

/// SHAKE128-based XOF: absorb seed || dom (little endian), squeeze forever.
pub struct Shake128Rng {
    reader: sha3::Shake128Reader,
}

impl Shake128Rng {
    pub fn new(seed: &[u8; 32], dom: u64) -> Self {
        let mut hasher = sha3::Shake128::default();
        hasher.update(seed);
        hasher.update(&dom.to_le_bytes());
        Shake128Rng {
            reader: hasher.finalize_xof(),
        }
    }
}

impl Xof for Shake128Rng {
    fn fill(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }
}

/// AES256-CTR keystream: key = seed, IV = dom_le || counter (LE, from 0).
pub struct Aes256CtrRng {
    cipher: ctr::Ctr64LE<aes::Aes256>,
}

impl Aes256CtrRng {
    pub fn new(seed: &[u8; 32], dom: u64) -> Self {
        let mut iv = [0u8; 16];
        iv[..8].copy_from_slice(&dom.to_le_bytes());
        Aes256CtrRng {
            cipher: ctr::Ctr64LE::<aes::Aes256>::new(seed.into(), &iv.into()),
        }
    }
}

impl Xof for Aes256CtrRng {
    fn fill(&mut self, out: &mut [u8]) {
        out.fill(0);
        self.cipher.apply_keystream(out);
    }
}

/// The default RNG backend, matching lazer's reproducible-build default.
pub type DefaultRng = Shake128Rng;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shake128_nist_empty_message_vector() {
        // SHAKE128(""), first 32 bytes — standard KAT value.
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"");
        let mut reader = hasher.finalize_xof();
        let mut out = [0u8; 32];
        reader.read(&mut out);
        let expected = [
            0x7f, 0x9c, 0x2b, 0xa4, 0xe8, 0x8f, 0x82, 0x7d, 0x61, 0x60, 0x45, 0x50, 0x76, 0x05,
            0x85, 0x3e, 0xd7, 0x3b, 0x80, 0x93, 0xf6, 0xef, 0xbc, 0x88, 0xeb, 0x1a, 0x6e, 0xac,
            0xfa, 0x66, 0xef, 0x26,
        ];
        assert_eq!(out, expected);
    }

    #[test]
    fn deterministic_and_domain_separated() {
        let seed = [7u8; 32];
        let mut a = Shake128Rng::new(&seed, 0);
        let mut b = Shake128Rng::new(&seed, 0);
        let mut c = Shake128Rng::new(&seed, 1);
        let (mut xa, mut xb, mut xc) = ([0u8; 64], [0u8; 64], [0u8; 64]);
        a.fill(&mut xa);
        b.fill(&mut xb);
        c.fill(&mut xc);
        assert_eq!(xa, xb);
        assert_ne!(xa, xc);

        let mut d = Aes256CtrRng::new(&seed, 0);
        let mut e = Aes256CtrRng::new(&seed, 0);
        let mut f = Aes256CtrRng::new(&seed, 9);
        let (mut xd, mut xe, mut xf) = ([0u8; 64], [0u8; 64], [0u8; 64]);
        d.fill(&mut xd);
        e.fill(&mut xe);
        f.fill(&mut xf);
        assert_eq!(xd, xe);
        assert_ne!(xd, xf);
        assert_ne!(xa, xd);
    }

    #[test]
    fn streaming_equals_one_shot() {
        let seed = [42u8; 32];
        let mut one = Shake128Rng::new(&seed, 5);
        let mut big = [0u8; 48];
        one.fill(&mut big);

        let mut two = Shake128Rng::new(&seed, 5);
        let mut parts = [0u8; 48];
        two.fill(&mut parts[..7]);
        two.fill(&mut parts[7..29]);
        two.fill(&mut parts[29..]);
        assert_eq!(big, parts);
    }
}
