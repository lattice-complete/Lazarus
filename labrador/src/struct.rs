#[Derive(Debug, Clone)]
pub struct CommitmentParams {
    // amortized opening decomposition parts
    func: u64,
    // uniform decomposition parts
    func_u: u64,
    // quadratic garbage decomposition parts
    func_g: u64,
    // amortized opening decomposition basis
    basis: u64,
    // uniform decomposition basis
    basis_u: u64,
    // quadrstic garbage decomposition basis
    basis_g: u64,
    // inner commitment rank
    kappa_inner: u64,
    // outer commitment rank
    kappa_outer: u64,
}

#[Derive(Debug, Clone)]
pub struct Proof {
    // input witness multiplicity
    r: u64,
    // input witness rank
    n: u64,
    // input witness decomposition parts
    n_u: u64,
    // true for no outer commitment
    is_tail: bool,
    // commitment parameters
    commitment_params: CommitmentParams,
    // outer commitment 1 (kappa_outer)
    commitment_outer_1: u64, // TODO: change to the polynomial_mod_q type
    // JL matrix nonce
    jl_nonce: u64,
    // JL projection
    vec_p: [u64; 256],
    // int to pol lifting pols (LIFTS)
    bb: u64, // TODO: understand what this is and change to the polynomial_mod_q type
    // outer commitment 2 (kappa_outer)
    commitment_outer_2: u64, // TODO: change to the polynomial_mod_q type
    // output witness norm squared
    norm_sq: u64,
}

#[Derive(Debug, Clone)]
pub struct Statement {
    // total amortized multiplicity
    r: u64,
    // rank of amortized opening
    n: u64,
    // rank of outer commitment opnenings (inner commitments, garbage)
    m: u64,
    // true for no outer commitment
    is_tail: bool,
    // commitment parameters
    commitment_params: Vec<CommitmentParams>,
    // outer commitment 1 (kappa1) or inner commitments (r*kappa) if tail
    outer_commitment_1: u64, // TODO: change to the polynomial_ext type
    // outer commitment 2 (kappa1) or garbage terms ((r*r+r)/2) if tail
    outer_commitment_2: u64, // TODO: change to the polynomial_ext type
    // challenges (r)
    challenge: u64,
    // aggregated dot-product constraint
    constraint: Constraint,
    // norm bound
    beta_sq: u64,
    // hash
    hash: [u8; 16],
}

#[Derive(Debug, Clone)]
pub struct Constraint {
    // extension degree
    deg_ext: u64,
    // sparce matrix a
    a: [u64; 1], // TODO: change to the sparse matrix type
    phi: u64,    // TODO: change to the polynomial_ext type
    b: u64,      // TODO: change to the polynomial_ext type
}

#[Derive(Debug, Clone)]
pub struct Witness {
    // witness multiplicity
    r: u64,
    // witness rank
    n: u64,
    // output witness norm squared
    norm_sq: u64,
    s: u64, // TODO: change to polynomial type
}
