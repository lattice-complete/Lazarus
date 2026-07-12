//! Recursive proof composition (port of `pack.c`): keep self-folding while
//! it shrinks the message, then try one tail round, and ship the final
//! witness in the clear.

use crate::comkey::ComKey;
use crate::principal::{principal_verify, PrincipalStatement};
use crate::round::{prove_round, reduce_round, Input};
use crate::statement::{verify, RoundProof, Statement};
use crate::witness::Witness;

pub const MAX_ROUNDS: usize = 16;

/// A chained proof: one round per entry plus the final in-clear witness.
#[derive(Debug, Clone)]
pub struct CompositeProof {
    pub rounds: Vec<RoundProof>,
    pub final_witness: Witness,
}

/// Rough message sizes in ring elements, used only for the greedy
/// stop rule (the C uses entropy-based KB estimates; the greedy shape is
/// the same: stop when a round stops shrinking the remaining message).
fn proof_size(pf: &RoundProof) -> usize {
    pf.u1.len() + pf.u2.len() + pf.bb.len() + 4 /* p, scalars */
}

fn witness_size(wt: &Witness) -> usize {
    wt.parts.iter().map(|p| p.len()).sum()
}

/// Prove a principal statement recursively (C: composite_prove_principle).
pub fn composite_prove(
    key: &ComKey,
    pst: &PrincipalStatement,
    wt: &Witness,
) -> Result<CompositeProof, String> {
    let mut rounds = Vec::new();

    // first round: principal -> core
    let (mut st, mut cur_wt, pf) = prove_round(key, &Input::Principal(pst), wt, false)?;
    rounds.push(pf);
    let mut cur_size = witness_size(&cur_wt);

    // keep folding while it pays off
    while rounds.len() < MAX_ROUNDS {
        let (nst, nwt, pf) = prove_round(key, &Input::Core(&st), &cur_wt, false)?;
        let nsize = witness_size(&nwt) + proof_size(&pf);
        if nsize >= cur_size {
            break;
        }
        st = nst;
        cur_wt = nwt;
        cur_size = witness_size(&cur_wt);
        rounds.push(pf);
    }

    // try one tail round
    if rounds.len() < MAX_ROUNDS {
        let (nst, nwt, pf) = prove_round(key, &Input::Core(&st), &cur_wt, true)?;
        let nsize = witness_size(&nwt) + proof_size(&pf);
        if nsize < cur_size {
            st = nst;
            cur_wt = nwt;
            rounds.push(pf);
        }
    }

    let _ = st;
    Ok(CompositeProof {
        rounds,
        final_witness: cur_wt,
    })
}

/// Verify a composite proof against the principal statement
/// (C: composite_verify_principle).
pub fn composite_verify(
    key: &ComKey,
    pst: &PrincipalStatement,
    proof: &CompositeProof,
) -> Result<(), String> {
    if proof.rounds.is_empty() {
        // degenerate: witness in the clear against the principal statement
        return principal_verify(pst, &proof.final_witness);
    }
    let mut st: Statement = reduce_round(key, &Input::Principal(pst), &proof.rounds[0])?;
    for pf in &proof.rounds[1..] {
        st = reduce_round(key, &Input::Core(&st), pf)?;
    }
    verify(key, &st, &proof.final_witness)
}
