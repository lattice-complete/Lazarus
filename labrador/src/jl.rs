//! Modular Johnson-Lindenstrauss projection (port of `jlproj.c` /
//! `project()`): 256 random +-1 rows over the whole coefficient vector,
//! used to prove smallness, then collapsed into Z_q constant-term
//! constraints via the conjugation automorphism.
//!
//! The +-1 matrix is never materialized: each witness poly slot gets a
//! 256 x 64 bit block, streamed from a transcript-derived XOF (the C uses
//! AES128; deviation documented in docs/MIGRATION.md). Zero-padded layout
//! slots contribute zero blocks, mirroring the reference's memset.

use crate::transcript::{self, Hash};
use crate::Rq;
use algebra::rng::{Shake128Rng, Xof};
use algebra::zq::Zq;

pub const JL_ROWS: usize = 256;
/// Bytes of matrix bits per polynomial slot: 256 rows x 64 coeffs / 8.
pub const BLOCK_BYTES: usize = JL_ROWS * 64 / 8;

/// Per-part matrix bits, one block per polynomial of the part.
pub type JlMats = Vec<Vec<u8>>;

fn sign(bits: &[u8], row: usize, coeff: usize) -> i64 {
    let idx = row * 64 + coeff;
    if (bits[idx >> 3] >> (idx & 7)) & 1 == 1 {
        -1
    } else {
        1
    }
}

/// Draw the matrices for all real witness slots (part-major, poly-major).
pub fn draw_mats(h: &Hash, nonce: u64, part_lens: &[usize]) -> JlMats {
    let mut xof = Shake128Rng::new(h, nonce);
    part_lens
        .iter()
        .map(|&len| {
            let mut bits = vec![0u8; len * BLOCK_BYTES];
            xof.fill(&mut bits);
            bits
        })
        .collect()
}

/// p += Pi * s for one witness part (blocks laid poly-major).
pub fn project_add(p: &mut [i64; JL_ROWS], part: &algebra::ring::PolyVec<64>, bits: &[u8]) {
    for (j, poly) in part.polys.iter().enumerate() {
        let block = &bits[j * BLOCK_BYTES..(j + 1) * BLOCK_BYTES];
        let cent = poly.centered();
        for (row, pr) in p.iter_mut().enumerate() {
            let mut acc = 0i64;
            for (c, &v) in cent.iter().enumerate() {
                if v != 0 {
                    acc += sign(block, row, c) * v;
                }
            }
            *pr += acc;
        }
    }
}

pub fn proj_normsq(p: &[i64; JL_ROWS]) -> u128 {
    p.iter().map(|&x| (x as i128 * x as i128) as u128).sum()
}

fn next2power(a: u64) -> u64 {
    a.next_power_of_two()
}

/// Prover-side projection with the reference's rejection loop: retry with
/// a fresh nonce until every |p_i| < next2power(4 sqrt(normsq)) and
/// ||p||^2 <= 256 * normsq. Returns (nonce, p, mats).
pub fn project(
    h: &Hash,
    wt: &crate::witness::Witness,
) -> Result<(u64, [i64; JL_ROWS], JlMats), String> {
    let normsq = wt.total_normsq();
    if normsq > crate::params::jl_max_normsq() {
        return Err("total witness norm too big for JL projection".to_string());
    }
    let part_lens: Vec<usize> = wt.parts.iter().map(|p| p.len()).collect();
    let test = next2power((4.0 * (normsq as f64).sqrt()) as u64 + 1);
    let bound = 256u128 * normsq;

    let mut nonce = 0u64;
    loop {
        nonce += 1;
        let mats = draw_mats(h, nonce, &part_lens);
        let mut p = [0i64; JL_ROWS];
        for (part, bits) in wt.parts.iter().zip(mats.iter()) {
            project_add(&mut p, part, bits);
        }
        let coord_ok = p.iter().all(|&x| (x.unsigned_abs()) < test);
        if coord_ok && proj_normsq(&p) <= bound {
            return Ok((nonce, p, mats));
        }
    }
}

/// Collapse the JL relation into one Z_q constant-term constraint under
/// 256 uniform scalars alpha (C: collaps_jlproj_raw):
///   phi[slot] = sigma_{-1}( sum_k alpha_k * pi_{k,slot} ),
///   b = sum_k alpha_k p_k,
/// so that ct(<phi, s>) = b whenever p = Pi s.
pub struct JlCollapse {
    /// phi per flat layout slot (padding slots zero).
    pub phi: Vec<Rq>,
    /// claimed constant term.
    pub b: Zq,
}

pub fn collapse(
    h: &mut Hash,
    plan: &crate::params::Plan,
    mats: &JlMats,
    p: &[i64; JL_ROWS],
) -> JlCollapse {
    transcript::absorb(h, transcript::tags::JLCOLLAPS, &[]);
    let alphas = transcript::zq_scalars(h, 0, JL_ROWS);

    let mut phi = vec![Rq::zero(); plan.r * plan.n];
    for (part, bits) in mats.iter().enumerate() {
        let len = plan.n_in[part];
        for j in 0..len {
            let block = &bits[j * BLOCK_BYTES..(j + 1) * BLOCK_BYTES];
            let mut coeffs = [Zq::new(0); 64];
            for (c, coeff) in coeffs.iter_mut().enumerate() {
                let mut acc = Zq::new(0);
                for (row, &a) in alphas.iter().enumerate() {
                    if sign(block, row, c) == 1 {
                        acc += a;
                    } else {
                        acc = acc - a;
                    }
                }
                *coeff = acc;
            }
            phi[plan.slot(part, j)] = Rq::from_coeffs(coeffs).conjugate();
        }
    }

    let mut b = Zq::new(0);
    for (row, &a) in alphas.iter().enumerate() {
        b += a * Zq::from_i64(p[row]);
    }
    JlCollapse { phi, b }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Plan;
    use crate::witness::Witness;
    use algebra::ring::{Poly, PolyVec};

    fn small_witness() -> Witness {
        let mut c = [0i64; 64];
        c[0] = 2;
        c[13] = -1;
        c[40] = 1;
        let p = Poly::<64>::from_centered(&c);
        Witness::new(vec![
            PolyVec::from_polys(vec![p.clone(), Poly::monomial(3)]),
            PolyVec::from_polys(vec![p]),
        ])
    }

    #[test]
    fn projection_bounds_hold() {
        let wt = small_witness();
        let h = [9u8; 32];
        let (nonce, p, _) = project(&h, &wt).unwrap();
        assert!(nonce >= 1);
        assert!(proj_normsq(&p) <= 256 * wt.total_normsq());
        // deterministic re-draw gives same projection
        let mats = draw_mats(&h, nonce, &[2, 1]);
        let mut p2 = [0i64; JL_ROWS];
        for (part, bits) in wt.parts.iter().zip(mats.iter()) {
            project_add(&mut p2, part, bits);
        }
        assert_eq!(p, p2);
    }

    /// The collapse soundness identity: ct(<phi, s>) == b for honest p.
    #[test]
    fn collapse_constant_term_identity() {
        let wt = small_witness();
        let h0 = [11u8; 32];
        let (_, p, mats) = project(&h0, &wt).unwrap();

        // plan: join both parts into one group of rank 3, chunks of 3
        let plan = Plan {
            n_in: vec![2, 1],
            nu: vec![0, 1],
            r: 1,
            n: 3,
            part_base: vec![0, 2],
        };
        let sx = wt.assemble(&plan);

        let mut h = h0;
        let col = collapse(&mut h, &plan, &mats, &p);
        let acc = col
            .phi
            .iter()
            .zip(sx.iter())
            .fold(Rq::zero(), |a, (f, s)| &a + &(f * s));
        assert_eq!(acc.coeffs[0], col.b, "ct(<phi,s>) != <alpha,p>");
    }
}
