#[derive(Debug, Clone)]
pub struct CommitmentParams {
    // amortized opening decomposition parts
    pub func: u64,
    // uniform decomposition parts
    pub func_u: u64,
    // quadratic garbage decomposition parts
    pub func_g: u64,
    // amortized opening decomposition basis
    pub basis: u64,
    // uniform decomposition basis
    pub basis_u: u64,
    // quadrstic garbage decomposition basis
    pub basis_g: u64,
    // inner commitment rank
    pub kappa_inner: u64,
    // outer commitment rank
    pub kappa_outer: u64,
}

#[derive(Debug, Clone)]
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

impl Proof {
    pub fn new() -> Self {
        Proof {
            r: 0,
            n: 0,
            n_u: 0,
            is_tail: false,
            commitment_params: CommitmentParams {
                func: 0,
                func_u: 0,
                func_g: 0,
                basis: 0,
                basis_u: 0,
                basis_g: 0,
                kappa_inner: 0,
                kappa_outer: 0,
            },
            commitment_outer_1: 0,
            jl_nonce: 0,
            vec_p: [0; 256],
            bb: 0,
            commitment_outer_2: 0,
            norm_sq: 0,
        }
    }
}

#[derive(Debug, Clone)]
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
    pub commitment_params: Vec<CommitmentParams>,
    // outer commitment 1 (kappa1) or inner commitments (r*kappa) if tail
    pub outer_commitment_1: u64, // TODO: change to the polynomial_ext type
    // outer commitment 2 (kappa1) or garbage terms ((r*r+r)/2) if tail
    pub outer_commitment_2: u64, // TODO: change to the polynomial_ext type
    // challenges (r)
    pub challenge: u64,
    // aggregated dot-product constraint
    pub constraint: Constraint,
    // norm bound
    pub beta_sq: u64,
    // hash
    pub hash: [u8; 16],
}

#[derive(Debug, Clone)]
pub struct Constraint {
    // extension degree
    deg_ext: u64,
    // sparce matrix a
    a: [u64; 1], // TODO: change to the sparse matrix type
    phi: u64,    // TODO: change to the polynomial_ext type
    b: u64,      // TODO: change to the polynomial_ext type
}

#[derive(Debug, Clone)]
pub struct Witness {
    // witness multiplicity
    r: u64,
    // witness rank
    n: u64,
    // output witness norm squared
    norm_sq: u64,
    s: u64, // TODO: change to polynomial type
}
