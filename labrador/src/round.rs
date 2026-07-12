//! One LaBRADOR reduction round: commit, JL-project, lift & aggregate,
//! amortize (C: prove()/reduce() with the principal front end of
//! chihuahua's principle_prove/principle_reduce).
//!
//! Prover and verifier walk the same transcript; the verifier rebuilds the
//! aggregated constraint from public data and the proof message only.

use crate::comkey::{ComKey, KeyDomain};
use crate::jl;
use crate::params::{self, aux_rank, select_params, ComParams, Plan, LIFTS, OPNORM_T, SLACK};
use crate::principal::PrincipalStatement;
use crate::statement::{npairs, pair_idx, RoundProof, Statement};
use crate::transcript::{self, tags, Hash};
use crate::witness::Witness;
use crate::Rq;
use algebra::ring::PolyVec;
use algebra::zq::Zq;
use std::collections::BTreeMap;

/// What a round reduces: either a user-facing principal statement (first
/// round) or a previously reduced core statement.
pub enum Input<'a> {
    Principal(&'a PrincipalStatement),
    Core(&'a Statement),
}

impl Input<'_> {
    pub fn hash(&self) -> Hash {
        match self {
            Input::Principal(p) => p.h,
            Input::Core(s) => s.h,
        }
    }
    pub fn quadratic(&self) -> bool {
        match self {
            Input::Principal(p) => p.quadratic(),
            Input::Core(s) => s.quadratic,
        }
    }
    pub fn betasq(&self) -> u128 {
        match self {
            Input::Principal(p) => p.betasq,
            Input::Core(s) => s.betasq,
        }
    }
    /// Ranks the input witness must have.
    pub fn expected_ranks(&self) -> Vec<usize> {
        match self {
            Input::Principal(p) => p.n.clone(),
            Input::Core(s) => {
                let mut v = vec![s.n; s.cpp.f];
                if !s.tail {
                    v.push(s.m);
                }
                v
            }
        }
    }
}

fn ip(a: &[Rq], b: &[Rq]) -> Rq {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .fold(Rq::zero(), |acc, (x, y)| &acc + &(x * y))
}

fn scale_slice(v: &[Rq], s: Zq) -> Vec<Rq> {
    v.iter().map(|p| p * s).collect()
}

/// Accumulator for the output aggregated constraint over the r x n layout.
struct Agg {
    phi: Vec<Rq>,
    a: BTreeMap<(usize, usize), Rq>,
    b: Rq,
}

impl Agg {
    fn new(slots: usize) -> Self {
        Agg {
            phi: vec![Rq::zero(); slots],
            a: BTreeMap::new(),
            b: Rq::zero(),
        }
    }
    fn add_phi(&mut self, slot: usize, v: &Rq) {
        self.phi[slot] = &self.phi[slot] + v;
    }
    fn add_a(&mut self, i: usize, j: usize, v: &Rq) {
        let e = self.a.entry((i, j)).or_insert_with(Rq::zero);
        *e = &*e + v;
    }
    fn a_vec(&self) -> Vec<(usize, usize, Rq)> {
        self.a
            .iter()
            .filter(|(_, v)| !v.is_zero())
            .map(|(&(i, j), v)| (i, j, v.clone()))
            .collect()
    }
}

/// Zq-constant-term constraint being lifted (JL collapse + deg-0 sparse
/// constraints of a principal input).
struct LiftCnst {
    phi: Vec<Rq>,
    a: Vec<(usize, usize, Rq)>,
    b_ct: Zq,
}

fn eval_full(l: &LiftCnst, sx: &[Rq], plan: &Plan) -> Rq {
    let n = plan.n;
    let mut acc = Rq::zero();
    for (slot, f) in l.phi.iter().enumerate() {
        if !f.is_zero() {
            acc = &acc + &(f * &sx[slot]);
        }
    }
    for (i, j, aij) in &l.a {
        let gij = ip(&sx[i * n..(i + 1) * n], &sx[j * n..(j + 1) * n]);
        acc = &acc + &(aij * &gij);
    }
    acc
}

/// Fold the deg-0 constraints of a principal statement into the lift
/// constraint under fresh uniform Z_q scalars (C: collaps_sparsecnst).
fn collapse_deg0(h: &mut Hash, pst: &PrincipalStatement, plan: &Plan, lift: &mut LiftCnst) {
    let deg0: Vec<_> = pst.deg0().collect();
    if deg0.is_empty() {
        return;
    }
    transcript::absorb(h, tags::SPARSE, &[]);
    let scalars = transcript::zq_scalars(h, 0, deg0.len());
    for (cnst, &s) in deg0.iter().zip(scalars.iter()) {
        lift.b_ct += s * cnst.b.coeffs[0];
        for e in &cnst.entries {
            for (lidx, f) in e.phi.iter().enumerate() {
                let slot = plan.slot(e.part, e.off + lidx);
                lift.phi[slot] = &lift.phi[slot] + &(f * s);
            }
        }
        for (i, j, aij) in &cnst.a {
            assert!(
                plan.parts_unsplit(),
                "quadratic sparse constraints need unsplit parts"
            );
            lift.a.push((*i, *j, aij * s));
        }
    }
}

/// Fold the full-ring constraints of a principal statement into the
/// aggregated constraint (C: aggregate_sparsecnst, flattened).
fn aggregate_principal(h: &mut Hash, pst: &PrincipalStatement, plan: &Plan, agg: &mut Agg) {
    transcript::absorb(h, tags::AGG, &[]);
    for (idx, cnst) in pst.ring_cnsts().enumerate() {
        let alpha = transcript::quarternary(h, idx as u64, 1).pop().unwrap();
        agg.b = &agg.b + &(&alpha * &cnst.b);
        for e in &cnst.entries {
            for (lidx, f) in e.phi.iter().enumerate() {
                let slot = plan.slot(e.part, e.off + lidx);
                agg.add_phi(slot, &(&alpha * f));
            }
        }
        for (i, j, aij) in &cnst.a {
            assert!(
                plan.parts_unsplit(),
                "quadratic sparse constraints need unsplit parts"
            );
            agg.add_a(*i, *j, &(&alpha * aij));
        }
    }
}

/// Fold a previously reduced core statement's verification equations into
/// the aggregated constraint over the current layout (C: aggregate()).
fn aggregate_core(h: &mut Hash, ist: &Statement, plan: &Plan, key: &ComKey, agg: &mut Agg) {
    let cpp = &ist.cpp;
    let r_in = ist.r;
    let l_in = cpp.fu * cpp.kappa;
    let pairs_in = npairs(r_in);
    let (toff, goff, hoff) = ist.aux_offsets();
    let vpart = cpp.f; // index of the aux part in the current input witness

    transcript::absorb(h, tags::AGG, &[]);
    let chal = transcript::quarternary(h, 0, 2 * cpp.kappa1 + cpp.kappa + 2);
    let alpha = &chal[..cpp.kappa1];
    let beta = &chal[cpp.kappa1..2 * cpp.kappa1];
    let gamma = &chal[2 * cpp.kappa1..2 * cpp.kappa1 + cpp.kappa];
    let delta0 = &chal[2 * cpp.kappa1 + cpp.kappa];
    let delta1 = &chal[2 * cpp.kappa1 + cpp.kappa + 1];

    /* B v_t (+ C v_g) = u1 */
    for i in 0..r_in {
        let cols = key.collapse(KeyDomain::B, alpha, i * l_in, l_in);
        for (l, cval) in cols.iter().enumerate() {
            agg.add_phi(plan.slot(vpart, toff + i * l_in + l), cval);
        }
    }
    if cpp.fg > 0 {
        let cols = key.collapse(KeyDomain::C, alpha, 0, cpp.fg * pairs_in);
        for (l, cval) in cols.iter().enumerate() {
            agg.add_phi(plan.slot(vpart, goff + l), cval);
        }
    }
    agg.b = &agg.b + &ip(alpha, &ist.u1.polys);

    /* D v_h = u2 */
    let cols = key.collapse(KeyDomain::D, beta, 0, cpp.fu * pairs_in);
    for (l, cval) in cols.iter().enumerate() {
        agg.add_phi(plan.slot(vpart, hoff + l), cval);
    }
    agg.b = &agg.b + &ip(beta, &ist.u2.polys);

    /* z digits: phi[z_i] -= 2^(i b) (phi_in + gamma^T A) */
    let mut tmp: Vec<Rq> = key.collapse(KeyDomain::A, gamma, 0, ist.n);
    for (l, t) in tmp.iter_mut().enumerate() {
        *t = &*t + &ist.phi.polys[l];
    }
    for digit in 0..cpp.f {
        if digit > 0 {
            let s = Zq::new(2).pow(cpp.b as usize);
            tmp = scale_slice(&tmp, s);
        }
        for (l, t) in tmp.iter().enumerate() {
            let neg = -t;
            agg.add_phi(plan.slot(digit, l), &neg);
        }
    }

    /* + sum_i c_i t_i : phi[v_t chunk i digit j] += 2^(j bu) c_i gamma */
    for i in 0..r_in {
        let mut cvec: Vec<Rq> = gamma.iter().map(|g| g * &ist.c[i]).collect();
        for digit in 0..cpp.fu {
            if digit > 0 {
                let s = Zq::new(2).pow(cpp.bu as usize);
                cvec = scale_slice(&cvec, s);
            }
            for (k, cv) in cvec.iter().enumerate() {
                agg.add_phi(
                    plan.slot(vpart, toff + i * l_in + digit * cpp.kappa + k),
                    cv,
                );
            }
        }
    }

    /* pair weights w_k = c_i c_j */
    let mut w: Vec<Rq> = Vec::with_capacity(pairs_in);
    for i in 0..r_in {
        for j in i..r_in {
            w.push(&ist.c[i] * &ist.c[j]);
        }
    }

    /* quadratic garbage block (only if the input had quadratic terms) */
    if ist.quadratic {
        let mut gw = vec![Rq::zero(); pairs_in];
        for i in 0..r_in {
            for j in i..r_in {
                let k = pair_idx(i, j, r_in);
                let mut v = delta0 * &w[k];
                if i != j {
                    v = &v + &v; // 2 delta0 c_i c_j off-diagonal
                }
                gw[k] = v;
            }
        }
        for (i, j, aij) in &ist.a {
            let k = pair_idx(*i, *j, r_in);
            gw[k] = &gw[k] + &(delta1 * aij);
        }
        for digit in 0..cpp.fg {
            let s = Zq::new(2).pow((digit as u32 * cpp.bg) as usize);
            for (k, g) in gw.iter().enumerate() {
                agg.add_phi(plan.slot(vpart, goff + digit * pairs_in + k), &(g * s));
            }
        }

        /* output quadratic form: -delta0 <z,z> over the z-digit chunks */
        let nu0 = plan.nu[0];
        assert!(
            plan.nu[..cpp.f].iter().all(|&x| x == nu0),
            "z digits must split uniformly"
        );
        for di in 0..cpp.f {
            for dk in di..cpp.f {
                let scale = -2i64 + i64::from(di == dk); // -1 diag, -2 off
                let pow = Zq::new(2).pow(((di + dk) as u32 * cpp.b) as usize);
                let coeff = &(delta0 * pow) * Zq::from_i64(scale);
                for j in 0..nu0 {
                    agg.add_a(di * nu0 + j, dk * nu0 + j, &coeff);
                }
            }
        }
    }

    /* linear garbage block: weights c_i c_j + delta1 on the diagonal */
    for i in 0..r_in {
        let k = pair_idx(i, i, r_in);
        w[k] = &w[k] + delta1;
    }
    for digit in 0..cpp.fu {
        let s = Zq::new(2).pow((digit as u32 * cpp.bu) as usize);
        for (k, wv) in w.iter().enumerate() {
            agg.add_phi(plan.slot(vpart, hoff + digit * pairs_in + k), &(wv * s));
        }
    }

    /* b */
    agg.b = &agg.b + &(delta1 * &ist.b);
}

/// Everything commit() produces (in tail mode u1 carries the raw inner
/// commitments and quadratic garbage directly).
struct Committed {
    u1: PolyVec<64>,
    vt: Vec<Rq>,
    vg: Vec<Rq>,
}

fn commit(
    key: &ComKey,
    plan: &Plan,
    cpp: &ComParams,
    tail: bool,
    quadratic: bool,
    sx: &[Rq],
) -> Committed {
    let (r, n) = (plan.r, plan.n);
    let l = cpp.fu * cpp.kappa;
    let pairs = npairs(r);
    let mut vt = vec![Rq::zero(); if tail { 0 } else { r * l }];
    let mut u1acc = PolyVec::zero(cpp.kappa1);
    let mut t_raw = Vec::new();

    for i in 0..r {
        let chunk = &sx[i * n..(i + 1) * n];
        let ti = key.mul(KeyDomain::A, cpp.kappa, 0, chunk);
        if tail {
            t_raw.extend(ti.polys);
        } else {
            let digits = ti.decompose_topmost(1u64 << cpp.bu, cpp.fu);
            for (d, dv) in digits.iter().enumerate() {
                for (k, p) in dv.polys.iter().enumerate() {
                    vt[i * l + d * cpp.kappa + k] = p.clone();
                }
            }
            u1acc = &u1acc + &key.mul(KeyDomain::B, cpp.kappa1, i * l, &vt[i * l..(i + 1) * l]);
        }
    }

    let mut vg = Vec::new();
    let mut g_raw = Vec::new();
    if quadratic {
        let mut g = Vec::with_capacity(pairs);
        for i in 0..r {
            for j in i..r {
                g.push(ip(&sx[i * n..(i + 1) * n], &sx[j * n..(j + 1) * n]));
            }
        }
        if tail {
            g_raw = g;
        } else {
            let gv = PolyVec::from_polys(g).decompose_topmost(1u64 << cpp.bg, cpp.fg.max(1));
            vg = gv.into_iter().flat_map(|d| d.polys).collect();
            u1acc = &u1acc + &key.mul(KeyDomain::C, cpp.kappa1, 0, &vg);
        }
    }

    let u1 = if tail {
        let mut v = t_raw;
        v.extend(g_raw);
        PolyVec::from_polys(v)
    } else {
        u1acc
    };
    Committed { u1, vt, vg }
}

fn p_bytes(p: &[i64]) -> Vec<u8> {
    p.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Prover side of one reduction (C: prove / principle_prove).
pub fn prove_round(
    key: &ComKey,
    input: &Input,
    iwt: &Witness,
    tail: bool,
) -> Result<(Statement, Witness, RoundProof), String> {
    if let Input::Core(ist) = input {
        if ist.tail {
            return Err("cannot reduce a tail statement further".into());
        }
    }
    let quadratic = input.quadratic();
    let rp = select_params(iwt, quadratic, tail)?;
    let (plan, cpp) = (&rp.plan, &rp.cpp);
    let (r, n) = (plan.r, plan.n);
    let m = aux_rank(r, cpp);
    let pairs = npairs(r);
    let sx = iwt.assemble(plan);
    let mut h = input.hash();

    /* commit */
    let com = commit(key, plan, cpp, tail, quadratic, &sx);
    transcript::absorb_polyvec(&mut h, tags::U1, &com.u1);

    /* project */
    let (jlnonce, p, mats) = jl::project(&h, iwt)?;
    transcript::absorb(&mut h, tags::JLPROJ, &[&p_bytes(&p)]);

    /* lift & aggregate Zq constraints */
    let mut agg = Agg::new(r * n);
    let mut bb = Vec::with_capacity(LIFTS);
    for li in 0..LIFTS {
        let col = jl::collapse(&mut h, plan, &mats, &p);
        let mut lift = LiftCnst {
            phi: col.phi,
            a: Vec::new(),
            b_ct: col.b,
        };
        if let Input::Principal(pst) = input {
            collapse_deg0(&mut h, pst, plan, &mut lift);
        }
        let b_full = eval_full(&lift, &sx, plan);
        debug_assert_eq!(b_full.coeffs[0], lift.b_ct, "lift ct identity broken");
        let mut bbi = b_full.clone();
        bbi.coeffs[0] = Zq::new(0);
        transcript::absorb_poly(&mut h, tags::LIFT, &bbi);
        let alpha = transcript::quarternary(&h, li as u64, 1).pop().unwrap();
        for (slot, f) in lift.phi.iter().enumerate() {
            if !f.is_zero() {
                agg.add_phi(slot, &(&alpha * f));
            }
        }
        for (i, j, aij) in &lift.a {
            agg.add_a(*i, *j, &(&alpha * aij));
        }
        agg.b = &agg.b + &(&alpha * &b_full);
        bb.push(bbi);
    }

    /* aggregate full-ring constraints */
    match input {
        Input::Principal(pst) => aggregate_principal(&mut h, pst, plan, &mut agg),
        Input::Core(ist) => aggregate_core(&mut h, ist, plan, key, &mut agg),
    }

    /* amortize */
    let phi_rows: Vec<&[Rq]> = (0..r).map(|i| &agg.phi[i * n..(i + 1) * n]).collect();
    let chunk = |i: usize| &sx[i * n..(i + 1) * n];

    let (u2, c, vh) = if !tail {
        let mut hvals = Vec::with_capacity(pairs);
        for i in 0..r {
            for j in i..r {
                if i == j {
                    hvals.push(ip(phi_rows[i], chunk(i)));
                } else {
                    hvals.push(&ip(phi_rows[i], chunk(j)) + &ip(phi_rows[j], chunk(i)));
                }
            }
        }
        let hv = PolyVec::from_polys(hvals).decompose_topmost(1u64 << cpp.bu, cpp.fu);
        let vh: Vec<Rq> = hv.into_iter().flat_map(|d| d.polys).collect();
        let u2 = key.mul(KeyDomain::D, cpp.kappa1, 0, &vh);
        transcript::absorb_polyvec(&mut h, tags::U2, &u2);
        let c = transcript::challenges(&h, 0, r);
        (u2, c, vh)
    } else {
        // sequential tail amortization
        let mut u2vals = Vec::with_capacity(2 * r - 1);
        let h00 = ip(phi_rows[0], chunk(0));
        transcript::absorb_poly(&mut h, tags::TAIL_H, &h00);
        u2vals.push(h00);
        let mut c = vec![transcript::challenges(&h, 0, 1).pop().unwrap()];
        let mut zacc: Vec<Rq> = chunk(0).iter().map(|p| p * &c[0]).collect();
        let mut phiacc: Vec<Rq> = phi_rows[0].iter().map(|p| p * &c[0]).collect();
        #[allow(clippy::needless_range_loop)] // parallel index over chunk(i) + phi_rows[i]
        for i in 1..r {
            let hodd = &ip(phi_rows[i], &zacc) + &ip(&phiacc, chunk(i));
            let heven = ip(phi_rows[i], chunk(i));
            transcript::absorb(&mut h, tags::TAIL_H, &[&hodd.to_bytes(), &heven.to_bytes()]);
            u2vals.push(hodd);
            u2vals.push(heven);
            let ci = transcript::challenges(&h, 0, 1).pop().unwrap();
            for (za, s) in zacc.iter_mut().zip(chunk(i).iter()) {
                *za = &*za + &(s * &ci);
            }
            for (pa, f) in phiacc.iter_mut().zip(phi_rows[i].iter()) {
                *pa = &*pa + &(f * &ci);
            }
            c.push(ci);
        }
        (PolyVec::from_polys(u2vals), c, Vec::new())
    };

    /* z and the output witness */
    let mut z = vec![Rq::zero(); n];
    #[allow(clippy::needless_range_loop)] // parallel index over chunk(i) + c[i]
    for i in 0..r {
        for (zl, s) in z.iter_mut().zip(chunk(i).iter()) {
            *zl = &*zl + &(s * &c[i]);
        }
    }
    let zdigits = PolyVec::from_polys(z).decompose_topmost(1u64 << cpp.b, cpp.f);
    let mut parts: Vec<PolyVec<64>> = zdigits;
    if !tail {
        let mut v = com.vt.clone();
        v.extend(com.vg.iter().cloned());
        v.extend(vh.iter().cloned());
        debug_assert_eq!(v.len(), m);
        parts.push(PolyVec::from_polys(v));
    }
    let owt = Witness::new(parts);
    let betasq = owt.total_normsq();

    /* folded phi for the output statement */
    let mut phi_folded = vec![Rq::zero(); n];
    for i in 0..r {
        for (pf, f) in phi_folded.iter_mut().zip(phi_rows[i].iter()) {
            *pf = &*pf + &(f * &c[i]);
        }
    }

    let ost = Statement {
        r,
        n,
        m,
        tail,
        quadratic,
        cpp: *cpp,
        u1: com.u1.clone(),
        u2: u2.clone(),
        c,
        a: agg.a_vec(),
        phi: PolyVec::from_polys(phi_folded),
        b: agg.b.clone(),
        betasq,
        h,
    };
    let proof = RoundProof {
        n_in: plan.n_in.clone(),
        nu: plan.nu.clone(),
        cpp: *cpp,
        tail,
        quadratic,
        u1: com.u1,
        jlnonce,
        p: p.to_vec(),
        bb,
        u2,
        normsq: betasq,
    };
    Ok((ost, owt, proof))
}

/// Verifier side of one reduction (C: reduce / principle_reduce).
pub fn reduce_round(key: &ComKey, input: &Input, pf: &RoundProof) -> Result<Statement, String> {
    /* validate shape */
    if let Input::Core(ist) = input {
        if ist.tail {
            return Err("cannot reduce a tail statement further".into());
        }
    }
    let expected = input.expected_ranks();
    if pf.n_in != expected {
        return Err("proof witness ranks do not match statement".into());
    }
    if pf.quadratic != input.quadratic() {
        return Err("quadratic flag mismatch".into());
    }
    if pf.p.len() != 256 || pf.bb.len() != LIFTS {
        return Err("malformed proof".into());
    }
    let cpp = &pf.cpp;
    if cpp.f == 0 || cpp.fu == 0 || cpp.b == 0 || cpp.bu == 0 || cpp.kappa == 0 {
        return Err("degenerate parameters".into());
    }
    if pf.quadratic && !pf.tail && cpp.fg == 0 {
        return Err("missing quadratic garbage parameters".into());
    }

    /* rebuild the plan and validate the marker pattern */
    let r_in = pf.n_in.len();
    let markers: Vec<bool> = (0..r_in).map(|i| pf.quadratic || i == r_in - 1).collect();
    for (i, &marker) in markers.iter().enumerate() {
        if marker != (pf.nu[i] > 0) {
            return Err("invalid repartition".into());
        }
    }
    let plan = rebuild_plan(&pf.n_in, &pf.nu);
    let (r, n) = (plan.r, plan.n);
    let total_in: usize = pf.n_in.iter().sum();
    if r == 0 || n == 0 || r * n > 2 * total_in + 64 * r {
        return Err("implausible repartition".into());
    }
    let m = aux_rank(r, cpp);
    let pairs = npairs(r);
    if !pf.tail {
        if pf.u1.len() != cpp.kappa1 || pf.u2.len() != cpp.kappa1 {
            return Err("wrong commitment rank".into());
        }
    } else {
        let want_u1 = r * cpp.kappa + if pf.quadratic { pairs } else { 0 };
        if pf.u1.len() != want_u1 || pf.u2.len() != 2 * r - 1 {
            return Err("wrong tail message length".into());
        }
    }

    let mut h = input.hash();
    transcript::absorb_polyvec(&mut h, tags::U1, &pf.u1);

    /* projection */
    let parr: [i64; 256] =
        pf.p.clone()
            .try_into()
            .map_err(|_| "malformed projection".to_string())?;
    let bound = 256u128 * params::jl_max_normsq().min(input.betasq());
    if jl::proj_normsq(&parr) > bound {
        return Err("witness projection longer than bound".into());
    }
    let mats = jl::draw_mats(&h, pf.jlnonce, &pf.n_in);
    transcript::absorb(&mut h, tags::JLPROJ, &[&p_bytes(&pf.p)]);

    /* lifts */
    let mut agg = Agg::new(r * n);
    for li in 0..LIFTS {
        let col = jl::collapse(&mut h, &plan, &mats, &parr);
        let mut lift = LiftCnst {
            phi: col.phi,
            a: Vec::new(),
            b_ct: col.b,
        };
        if let Input::Principal(pst) = input {
            collapse_deg0(&mut h, pst, &plan, &mut lift);
        }
        // b believed = prover's bb with the claimed constant term
        let mut b_ver = pf.bb[li].clone();
        if b_ver.coeffs[0].value != 0 {
            return Err("lift polynomial has nonzero constant term".into());
        }
        b_ver.coeffs[0] = lift.b_ct;
        transcript::absorb_poly(&mut h, tags::LIFT, &pf.bb[li]);
        let alpha = transcript::quarternary(&h, li as u64, 1).pop().unwrap();
        for (slot, f) in lift.phi.iter().enumerate() {
            if !f.is_zero() {
                agg.add_phi(slot, &(&alpha * f));
            }
        }
        for (i, j, aij) in &lift.a {
            agg.add_a(*i, *j, &(&alpha * aij));
        }
        agg.b = &agg.b + &(&alpha * &b_ver);
    }

    /* aggregate ring constraints */
    match input {
        Input::Principal(pst) => aggregate_principal(&mut h, pst, &plan, &mut agg),
        Input::Core(ist) => aggregate_core(&mut h, ist, &plan, key, &mut agg),
    }

    /* amortize (verifier side): re-derive challenges, fold phi */
    let phi_rows: Vec<Vec<Rq>> = (0..r)
        .map(|i| agg.phi[i * n..(i + 1) * n].to_vec())
        .collect();
    let c: Vec<Rq> = if !pf.tail {
        transcript::absorb_polyvec(&mut h, tags::U2, &pf.u2);
        transcript::challenges(&h, 0, r)
    } else {
        let mut c = Vec::with_capacity(r);
        transcript::absorb_poly(&mut h, tags::TAIL_H, &pf.u2.polys[0]);
        c.push(transcript::challenges(&h, 0, 1).pop().unwrap());
        for i in 1..r {
            transcript::absorb(
                &mut h,
                tags::TAIL_H,
                &[
                    &pf.u2.polys[2 * i - 1].to_bytes(),
                    &pf.u2.polys[2 * i].to_bytes(),
                ],
            );
            c.push(transcript::challenges(&h, 0, 1).pop().unwrap());
        }
        c
    };

    /* security of the commitments under the claimed norm */
    let norm = (pf.normsq as f64).sqrt();
    if !params::sis_secure(
        cpp.kappa,
        6.0 * OPNORM_T * SLACK * 2f64.powi(((cpp.f - 1) as u32 * cpp.b) as i32) * norm,
    ) {
        return Err("inner commitments not secure".into());
    }
    if !pf.tail && !params::sis_secure(cpp.kappa1, 2.0 * SLACK * norm) {
        return Err("outer commitments not secure".into());
    }

    let mut phi_folded = vec![Rq::zero(); n];
    for i in 0..r {
        for (pfo, f) in phi_folded.iter_mut().zip(phi_rows[i].iter()) {
            *pfo = &*pfo + &(f * &c[i]);
        }
    }

    Ok(Statement {
        r,
        n,
        m,
        tail: pf.tail,
        quadratic: pf.quadratic,
        cpp: *cpp,
        u1: pf.u1.clone(),
        u2: pf.u2.clone(),
        c,
        a: agg.a_vec(),
        phi: PolyVec::from_polys(phi_folded),
        b: agg.b.clone(),
        betasq: pf.normsq,
        h,
    })
}

fn rebuild_plan(n_in: &[usize], nu: &[usize]) -> Plan {
    let r: usize = nu.iter().sum();
    let mut n = 0usize;
    let mut acc = 0usize;
    for i in 0..n_in.len() {
        acc += n_in[i];
        if nu[i] > 0 {
            n = n.max(acc.div_ceil(nu[i]));
            acc = 0;
        }
    }
    let mut part_base = vec![0usize; n_in.len()];
    let mut chunk = 0usize;
    let mut off = 0usize;
    for i in 0..n_in.len() {
        part_base[i] = chunk * n + off;
        off += n_in[i];
        if nu[i] > 0 {
            chunk += nu[i];
            off = 0;
        }
    }
    Plan {
        n_in: n_in.to_vec(),
        nu: nu.to_vec(),
        r,
        n,
        part_base,
    }
}
