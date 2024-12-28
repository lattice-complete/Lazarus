use crate::algebra::RqMatrix;
use crate::algebra::Zq;

type RqMatrix2D = Vec<Vec<RqMatrix>>;
type RqMatrix3D = Vec<Vec<Vec<RqMatrix>>>;

fn setup_matrices(
    size_n: Zq,
    size_r: Zq,
    t1: Zq,
    t2: Zq,
    kappa: Zq,
    kappa1: Zq,
    kappa2: Zq,
) -> (RqMatrix, RqMatrix2D, RqMatrix3D, RqMatrix3D) {
    // Initialize matrix A
    let a_matrix = RqMatrix::new(kappa, size_n);

    // Initialize matrix B
    let b_matrix: Vec<Vec<RqMatrix>> = (0..size_r.value())
        .map(|_i| {
            (0..t1.value())
                .map(|_k| RqMatrix::new(kappa1, kappa))
                .collect()
        })
        .collect();

    // Initialize matrix C
    let c_matrix: Vec<Vec<Vec<RqMatrix>>> = (0..size_r.value())
        .map(|_i| {
            (0..size_r.value())
                .map(|_j| {
                    (0..t2.value())
                        .map(|_k| RqMatrix::new(kappa2, Zq::from(1)))
                        .collect()
                })
                .collect()
        })
        .collect();

    // Initialize matrix D with the same shape as C
    let d_matrix: Vec<Vec<Vec<RqMatrix>>> = (0..size_r.value())
        .map(|_i| {
            (0..size_r.value())
                .map(|_j| {
                    (0..t2.value())
                        .map(|_k| RqMatrix::new(kappa2, Zq::from(1)))
                        .collect()
                })
                .collect()
        })
        .collect();

    (a_matrix, b_matrix, c_matrix, d_matrix)
}

pub fn setup(
    size_n: Zq,
    size_r: Zq,
    t1: Zq,
    t2: Zq,
    kappa: Zq,
    kappa1: Zq,
    kappa2: Zq,
) -> (RqMatrix, RqMatrix2D, RqMatrix3D, RqMatrix3D) {
    setup_matrices(size_n, size_r, t1, t2, kappa, kappa1, kappa2)
    // 0. setup
    // public parameters after setup: [a_ij^(k), a_ij^(l), phi^(k), phi^(l), b^(k), b0(l)']

    // 1. setup constraints
    // 1.1 get s_i and do norm check
    // 1.1.1 get s_i = s_1 - s_r; r is the number of witness s
    // each s_i is a vector of ring elements

    // 1.1.2 get beta norm bound, refer paper page 26, theorem 6.3

    // 1.1.3 do check: sum of s_i norm <= beta_square

    // 1.2 calculate b^(k)
    // 1.2.1 calculate dot product ss = a_ij * <s_i, s_j> for all i, j
    // a_ij is the quadratic coefficient, phi^(k) is the linear coefficient
    // 1.2.2 calculate phi_s = <phi^(k), s_i> for all i
    // 1.2.3 calculate b^(k) = sum(ss) + sum(phi_s)

    // 1.3 calculate b'^(l)
    // 1.3.1 calculate dot product ss = a_ij' * <s_i, s_j> for all i, j
    // a_ij' is the quadratic coefficient, phi^(l)' is the linear coefficient
    // 1.3.2 calculate phi_s = <phi^(l)', s_i> for all i
    // 1.3.3 calculate b'^(l) = sum(ss) + sum(phi_s)

    // L = |F'| = ceiling(128 / logQ)
}
