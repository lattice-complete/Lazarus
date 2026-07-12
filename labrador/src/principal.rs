//! Principal statements: user-facing sparse dot-product constraint systems
//! (port of chihuahua's `sparsecnst` / `prncplstmnt`, flattened — a
//! vector-valued constraint of extension degree e in the C API becomes e
//! scalar ring constraints here).

use crate::transcript::Hash;
use crate::witness::Witness;
use crate::Rq;
use algebra::zq::Zq;

/// One block of a sparse linear form: applies `phi` to
/// witness part `part` starting at offset `off`.
#[derive(Debug, Clone)]
pub struct SparseEntry {
    pub part: usize,
    pub off: usize,
    pub phi: Vec<Rq>,
}

/// One sparse constraint over the witness parts:
///   sum_e <e.phi, s_[e.part][e.off..]> + sum a_ij <s_i, s_j> = b,
/// as full ring elements, or as constant terms only if `ct_only`.
#[derive(Debug, Clone)]
pub struct SparseConstraint {
    pub ct_only: bool,
    pub entries: Vec<SparseEntry>,
    /// Quadratic part over whole witness parts (i <= j).
    pub a: Vec<(usize, usize, Rq)>,
    pub b: Rq,
}

impl SparseConstraint {
    pub fn eval(&self, wt: &Witness) -> Rq {
        let mut acc = Rq::zero();
        for e in &self.entries {
            let part = &wt.parts[e.part];
            assert!(e.off + e.phi.len() <= part.len(), "entry out of range");
            for (l, f) in e.phi.iter().enumerate() {
                if !f.is_zero() {
                    acc = &acc + &(f * &part.polys[e.off + l]);
                }
            }
        }
        for (i, j, aij) in &self.a {
            let ip = wt.parts[*i].inner_product(&wt.parts[*j]);
            acc = &acc + &(aij * &ip);
        }
        acc
    }

    pub fn check(&self, wt: &Witness) -> bool {
        let v = self.eval(wt);
        if self.ct_only {
            v.coeffs[0] == self.b.coeffs[0]
        } else {
            v == self.b
        }
    }
}

/// A collection of sparse constraints plus a norm bound (C: prncplstmnt).
#[derive(Debug, Clone)]
pub struct PrincipalStatement {
    /// Witness part ranks.
    pub n: Vec<usize>,
    pub cnst: Vec<SparseConstraint>,
    pub betasq: u128,
    pub h: Hash,
}

impl PrincipalStatement {
    pub fn quadratic(&self) -> bool {
        self.cnst.iter().any(|c| !c.a.is_empty())
    }

    pub fn deg0(&self) -> impl Iterator<Item = &SparseConstraint> {
        self.cnst.iter().filter(|c| c.ct_only)
    }

    pub fn ring_cnsts(&self) -> impl Iterator<Item = &SparseConstraint> {
        self.cnst.iter().filter(|c| !c.ct_only)
    }
}

/// Direct verification (C: principle_verify).
pub fn principal_verify(pst: &PrincipalStatement, wt: &Witness) -> Result<(), String> {
    if wt.parts.len() != pst.n.len() {
        return Err("wrong witness arity".into());
    }
    for (i, part) in wt.parts.iter().enumerate() {
        if part.len() != pst.n[i] {
            return Err(format!("wrong rank of witness part {i}"));
        }
    }
    if wt.total_normsq() > pst.betasq {
        return Err("witness norm bigger than bound".into());
    }
    for (k, c) in pst.cnst.iter().enumerate() {
        if !c.check(wt) {
            return Err(format!("constraint {k} does not hold"));
        }
    }
    Ok(())
}

/// Convenience: build the constant-term evaluation target for a ct-only
/// constraint from a scalar.
pub fn ct_poly(v: Zq) -> Rq {
    let mut b = Rq::zero();
    b.coeffs[0] = v;
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::ring::{Poly, PolyVec};

    #[test]
    fn sparse_eval_and_check() {
        // witness: one part [s0, s1], constraint <(1, X), s> = s0 + X*s1
        let mut c = [0i64; 64];
        c[2] = 3;
        let s0 = Poly::<64>::from_centered(&c);
        let s1 = Poly::monomial(5);
        let wt = Witness::new(vec![PolyVec::from_polys(vec![s0.clone(), s1.clone()])]);

        let b = &s0 + &(&Poly::monomial(1) * &s1);
        let cnst = SparseConstraint {
            ct_only: false,
            entries: vec![SparseEntry {
                part: 0,
                off: 0,
                phi: vec![Poly::one(), Poly::monomial(1)],
            }],
            a: vec![],
            b,
        };
        assert!(cnst.check(&wt));

        // quadratic: <s, s> = a-term
        let ip = wt.parts[0].inner_product(&wt.parts[0]);
        let q = SparseConstraint {
            ct_only: false,
            entries: vec![],
            a: vec![(0, 0, Poly::one())],
            b: ip,
        };
        assert!(q.check(&wt));

        let pst = PrincipalStatement {
            n: vec![2],
            cnst: vec![cnst, q],
            betasq: wt.total_normsq(),
            h: [0u8; 32],
        };
        principal_verify(&pst, &wt).unwrap();
    }
}
