//! # labrador — recursive lattice proof system
//!
//! Rust implementation of the LaBRADOR proof system (Beullens–Seiler,
//! CRYPTO'23) over R_q = Z_q[X]/(X^64 + 1), q = 2^32 - 99, ported from the
//! reference C implementation (lattice-dogs/labrador) onto the Lazarus
//! `algebra` foundation layer.
//!
//! Layers, mirroring the C architecture:
//!
//! * [`principal`] — user-facing sparse dot-product constraint systems
//!   (C: chihuahua's `prncplstmnt`);
//! * [`round`] — one reduction: Ajtai commit, modular JL projection,
//!   Z_q-constraint lifting, constraint aggregation, amortization
//!   (C: `prove`/`reduce`);
//! * [`composite`] — the recursion driver: fold until the witness is small,
//!   finish with a tail round and an in-clear witness (C: pack.c);
//! * [`statement`] — reduced statements and their direct verification;
//! * [`params`] — MSIS security estimation and automatic per-round
//!   parameter selection (C: `sis_secure`, `init_proof`).
//!
//! The legacy one-round prototype lives in [`legacy`].

pub mod legacy;

pub mod comkey;
pub mod composite;
pub mod jl;
pub mod params;
pub mod principal;
pub mod round;
pub mod statement;
pub mod transcript;
pub mod witness;

/// The working ring: degree 64, q = 2^32 - 99.
pub type Rq = algebra::ring::Poly<64>;

pub use comkey::{ComKey, KeyDomain};
pub use composite::{composite_prove, composite_verify, CompositeProof};
pub use principal::{principal_verify, PrincipalStatement, SparseConstraint, SparseEntry};
pub use round::{prove_round, reduce_round, Input};
pub use statement::{verify, RoundProof, Statement};
pub use witness::Witness;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::principal::ct_poly;
    use algebra::ring::{Poly, PolyVec};
    use algebra::zq::Zq;

    /// Build a small linear principal statement with a known witness:
    /// two parts of rank 4 with narrow entries, one ring constraint
    /// <phi, s_0> = b, one ct constraint ct(<psi, s_1>) = y.
    fn toy_instance(seed: u64) -> (PrincipalStatement, Witness) {
        use algebra::rng::{Shake128Rng, Xof};
        use algebra::sampling::binomial::poly_cbd;

        let nrank = 4usize;
        let parts: Vec<PolyVec<64>> = (0..2)
            .map(|i| {
                PolyVec::from_polys(
                    (0..nrank)
                        .map(|j| poly_cbd(2, &[seed as u8; 32], (i * nrank + j) as u64))
                        .collect(),
                )
            })
            .collect();
        let wt = Witness::new(parts);

        // random-ish public phi / psi
        let mut xof = Shake128Rng::new(&[9u8; 32], seed);
        let phi: Vec<Rq> = (0..nrank)
            .map(|_| {
                let mut c = [0i64; 64];
                for x in c.iter_mut() {
                    let mut b = [0u8; 1];
                    xof.fill(&mut b);
                    *x = (b[0] % 5) as i64 - 2;
                }
                Poly::from_centered(&c)
            })
            .collect();
        let psi: Vec<Rq> = (0..nrank).map(|i| Poly::monomial(i * 3 % 64)).collect();

        let b = phi
            .iter()
            .zip(wt.parts[0].polys.iter())
            .fold(Rq::zero(), |acc, (f, s)| &acc + &(f * s));
        let yfull = psi
            .iter()
            .zip(wt.parts[1].polys.iter())
            .fold(Rq::zero(), |acc, (f, s)| &acc + &(f * s));

        let cnst_ring = SparseConstraint {
            ct_only: false,
            entries: vec![SparseEntry {
                part: 0,
                off: 0,
                phi,
            }],
            a: vec![],
            b,
        };
        let cnst_ct = SparseConstraint {
            ct_only: true,
            entries: vec![SparseEntry {
                part: 1,
                off: 0,
                phi: psi,
            }],
            a: vec![],
            b: ct_poly(yfull.coeffs[0]),
        };

        let pst = PrincipalStatement {
            n: vec![nrank, nrank],
            cnst: vec![cnst_ring, cnst_ct],
            betasq: wt.total_normsq(),
            h: [seed as u8; 32],
        };
        principal_verify(&pst, &wt).unwrap();
        (pst, wt)
    }

    #[test]
    fn one_round_reduces_and_verifies() {
        let (pst, wt) = toy_instance(1);
        let key = ComKey::new([42u8; 32]);
        let (st, owt, pf) = prove_round(&key, &Input::Principal(&pst), &wt, false).unwrap();

        // verifier reduces to the same statement
        let vst = reduce_round(&key, &Input::Principal(&pst), &pf).unwrap();
        assert_eq!(vst.h, st.h, "transcripts diverged");
        assert_eq!(vst.b, st.b, "aggregated constant diverged");
        assert_eq!(vst.phi, st.phi, "folded phi diverged");
        assert_eq!(vst.c, st.c);
        assert_eq!(vst.betasq, st.betasq);

        // and the output witness satisfies it
        verify(&key, &vst, &owt).unwrap();
    }

    #[test]
    fn tail_round_verifies() {
        let (pst, wt) = toy_instance(2);
        let key = ComKey::new([43u8; 32]);
        let (st, owt, pf) = prove_round(&key, &Input::Principal(&pst), &wt, true).unwrap();
        let vst = reduce_round(&key, &Input::Principal(&pst), &pf).unwrap();
        assert_eq!(vst.h, st.h);
        verify(&key, &vst, &owt).unwrap();
    }

    #[test]
    fn two_rounds_chain() {
        let (pst, wt) = toy_instance(3);
        let key = ComKey::new([44u8; 32]);
        let (st1, wt1, pf1) = prove_round(&key, &Input::Principal(&pst), &wt, false).unwrap();
        let (st2, wt2, pf2) = prove_round(&key, &Input::Core(&st1), &wt1, false).unwrap();

        let v1 = reduce_round(&key, &Input::Principal(&pst), &pf1).unwrap();
        let v2 = reduce_round(&key, &Input::Core(&v1), &pf2).unwrap();
        assert_eq!(v2.h, st2.h, "second-round transcripts diverged");
        verify(&key, &v2, &wt2).unwrap();
    }

    #[test]
    fn composite_roundtrip_and_tamper_rejection() {
        let (pst, wt) = toy_instance(4);
        let key = ComKey::new([45u8; 32]);
        let proof = composite_prove(&key, &pst, &wt).unwrap();
        composite_verify(&key, &pst, &proof).unwrap();

        // tamper: flip one coefficient of the final witness
        let mut bad = proof.clone();
        let old = bad.final_witness.parts[0].polys[0].coeffs[0];
        bad.final_witness.parts[0].polys[0].coeffs[0] = old + Zq::new(1);
        assert!(composite_verify(&key, &pst, &bad).is_err());

        // tamper: modify a lifted constraint polynomial (transcript forks)
        let mut bad = proof.clone();
        bad.rounds[0].bb[0].coeffs[1] += Zq::new(1);
        assert!(composite_verify(&key, &pst, &bad).is_err());

        // tamper: different statement (wrong b)
        let mut bad_pst = pst.clone();
        bad_pst.cnst[0].b = &bad_pst.cnst[0].b + &Rq::one();
        assert!(composite_verify(&key, &bad_pst, &proof).is_err());
    }

    #[test]
    fn quadratic_constraint_roundtrip() {
        // witness with a quadratic relation: <s_0, s_0> = b
        use algebra::sampling::binomial::poly_cbd;
        let nrank = 3usize;
        let parts = vec![PolyVec::from_polys(
            (0..nrank)
                .map(|j| poly_cbd::<64>(2, &[7u8; 32], j as u64))
                .collect(),
        )];
        let wt = Witness::new(parts);
        let b = wt.parts[0].inner_product(&wt.parts[0]);
        let pst = PrincipalStatement {
            n: vec![nrank],
            cnst: vec![SparseConstraint {
                ct_only: false,
                entries: vec![],
                a: vec![(0, 0, Rq::one())],
                b,
            }],
            betasq: wt.total_normsq(),
            h: [5u8; 32],
        };
        principal_verify(&pst, &wt).unwrap();

        let key = ComKey::new([46u8; 32]);
        let (st, owt, pf) = prove_round(&key, &Input::Principal(&pst), &wt, false).unwrap();
        assert!(st.quadratic);
        let vst = reduce_round(&key, &Input::Principal(&pst), &pf).unwrap();
        assert_eq!(vst.h, st.h);
        verify(&key, &vst, &owt).unwrap();

        // and chain one more quadratic round
        let (st2, owt2, pf2) = prove_round(&key, &Input::Core(&st), &owt, false).unwrap();
        let v2 = reduce_round(&key, &Input::Core(&vst), &pf2).unwrap();
        assert_eq!(v2.h, st2.h);
        verify(&key, &v2, &owt2).unwrap();
    }
}
