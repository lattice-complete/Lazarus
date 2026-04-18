// src/verifier.rs
use crate::core::{aggregation::*, decompose::*, norm::*, outer_commitment::*};
use crate::statement::Statement;
use crate::transcript::Transcript;
use algebra::{polynomial_ring::PolynomialRing, rq_matrix::RqMatrix, utils::*, zq::Zq};
use profiler_macro::time_profiler;

// What does the verifier already know? (Statement)
// 1. public paraments: [a_ij^(k), a_ij^(l), phi^(k), phi^(l), b^(k), b0(l)']
// 2. generate from V: [PI_i, psi^(k), omega^(k), vec<alpha>, vec<beta>, c_i]
// What are sent to the verifier?
// [u1, p, b^{''(k)},u2, z, t_i, g_ij, h_ij]

#[time_profiler]
#[allow(clippy::needless_range_loop)] // formulas index parallel arrays g[i][j], h[i][j]; iterators obscure the math
pub fn verify(
    st: Statement,
    tr: Transcript,
    a_matrix: &RqMatrix,
    b_matrix: &[Vec<RqMatrix>],
    c_matrix: &[Vec<Vec<RqMatrix>>],
    d_matrix: &[Vec<Vec<RqMatrix>>],
) {
    // same parameters as in the prover
    let _size_r = Zq::new(3); // r: Number of witness elements
    let size_n = Zq::new(5); // n
    let basis = Zq::new(10);
    let digits = Zq::new(3); // t1
    let t1 = digits;
    let t2 = digits;
    let _kappa = Zq::new(3); // Example size
    let kappa1 = Zq::from(5);
    let kappa2 = Zq::from(5);
    let lambda = Zq::new(128);
    let double_lambda = lambda * Zq::new(2);
    let log_q = Zq::new(32);
    let deg_bound_d = Zq::new(8); // random polynomial degree bound
    let new_beta = Zq::new(250); // Example value for beta

    let Statement {
        a_constraint,
        phi_constraint,
        b_constraint,
        a_constraint_ct,
        phi_constraint_ct,
        b_constraint_ct: _,
    } = st;

    let Transcript {
        u1,
        pai,
        p: _,
        psi,
        omega,
        b_ct_aggr,
        alpha,
        beta,
        u2,
        c,
        z,
        t,
        g,
        h,
    } = tr;
    let size_r = Zq::from(g.len());
    // 1. check g_ij ?= g_ji
    for i in 0..size_r.value() {
        for j in (i + 1)..size_r.value() {
            assert_eq!(
                g[i][j], g[j][i],
                "g_ij is not equal to g_ji at indices ({}, {})",
                i, j
            );
        }
    }
    // 2. check h_ij ?= h_ji
    for i in 0..size_r.value() {
        for j in (i + 1)..size_r.value() {
            assert_eq!(
                h[i][j], h[j][i],
                "h_ij is not equal to h_ji at indices ({}, {})",
                i, j
            );
        }
    }
    // 3. check if norm sum of z, t, g, h <= beta'^2
    // 3.1 decompose z, t, g, h to basis form
    let z_basis_form = poly_vec_decompose_and_aggregate(&z, basis, digits);
    let all_t_i_basis_form_aggregated = poly_matrix_decompose_and_aggregate(&t, basis, digits);
    let g_matrix_aggregated = poly_matrix_decompose_and_aggregate(&g, basis, t2);
    let h_gar_poly_basis_form_aggregated = poly_matrix_decompose_and_aggregate(&h, basis, digits);
    let norm_z = poly_matrix_norm_squared(&z_basis_form);
    let norm_t = poly_3d_norm_squared(&all_t_i_basis_form_aggregated);
    let norm_g = poly_3d_norm_squared(&g_matrix_aggregated);
    let norm_h = poly_3d_norm_squared(&h_gar_poly_basis_form_aggregated);
    let norm_sum = norm_z + norm_t + norm_g + norm_h;
    println!("Verifier: Check norms of decomposed inner commitments");
    assert!(norm_sum <= new_beta.pow(2));

    println!("Verifier: Check amortized opening of inner commitments");
    // 4. check if Az is valid
    let a_times_z: Vec<PolynomialRing> = matrix_poly_times_poly_vector(&a_matrix.values, &z);
    // calculate sum(ci * ti)
    let sum_c_times_t: Vec<PolynomialRing> = inner_product_poly_matrix_and_poly_vector(&t, &c);
    assert_eq!(a_times_z, sum_c_times_t);

    println!("Verifier: Check aggregated innerproduct constraints");
    // 5. check if <z, z> ?= sum(g_ij * c_i * c_j)
    let z_z_inner_product = inner_product_polynomial_ring_vector(&z, &z);

    let mut sum_g_ij_c_i_c_j = zero_poly();
    for i in 0..size_r.value() {
        for j in 0..size_r.value() {
            let g_ij = &g[i][j]; // Borrow g[i][j] instead of moving it
            let c_i = &c[i];
            let c_j = &c[j];
            sum_g_ij_c_i_c_j = sum_g_ij_c_i_c_j + (g_ij * c_i * c_j);
        }
    }

    assert_eq!(z_z_inner_product, sum_g_ij_c_i_c_j);

    println!("Verifier: Check aggregated linear constraints");
    // 6. check if sum(<phi_i, z> * c_i) ?= sum(h_ij * c_i * c_j)
    // aggregation parameters
    let size_k = Zq::new(lambda.value() / log_q.value());
    let constraint_num_l = Zq::new(5); // Define L
    let constraint_num_k = Zq::new(5);

    // 6.1 caculate b^{''(k)}
    // 6.1.1 calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
    let a_ct_aggr = compute_aggr_ct_constraint_a(
        &a_constraint_ct,
        &psi,
        size_k,
        size_r,
        constraint_num_l,
        deg_bound_d,
    );
    assert_eq!(a_ct_aggr.len(), size_k.value());
    assert_eq!(a_ct_aggr[0].len(), size_r.value());
    assert_eq!(a_ct_aggr[0][0].len(), size_r.value());
    // 6.1.2 calculate phi_i^{''(k)} =
    //       sum(psi_l^(k) * phi_i^{'(l)}) for all l = 1..L
    //       + sum(omega_j^(k) * sigma_{-1} * pi_i^{j)) for all j = 1..256
    let phi_ct_aggr = compute_aggr_ct_constraint_phi(
        &phi_constraint_ct,
        &pai,
        size_k,
        size_r,
        constraint_num_l,
        deg_bound_d,
        size_n,
        double_lambda,
        &psi,
        &omega,
    );
    assert_eq!(phi_ct_aggr.len(), size_k.value());
    assert_eq!(phi_ct_aggr[0].len(), size_r.value());

    // b_ct_aggr does not need to be calculated here, it's from prover
    assert_eq!(b_ct_aggr.len(), size_k.value());

    let a_aggr = compute_aggr_constraint_a(
        &a_constraint,
        &a_ct_aggr,
        constraint_num_k,
        &alpha,
        &beta,
        size_r,
        size_k,
    );
    assert_eq!(a_aggr.len(), size_r.value());
    assert_eq!(a_aggr[0].len(), size_r.value());

    let phi_aggr = compute_aggr_constraint_phi(
        &phi_constraint,
        &phi_ct_aggr,
        constraint_num_k,
        &alpha,
        &beta,
        size_r,
        size_n,
        deg_bound_d,
        size_k,
    );

    let b_aggr = compute_aggr_constraint_b(
        &b_constraint,
        &b_ct_aggr,
        constraint_num_k,
        &alpha,
        &beta,
        size_k,
    );

    // check if sum(<phi_i, z> * c_i) ?= sum(h_ij * c_i * c_j)
    // calculate sum(<phi_i, z> * c_i)
    let sum_phi_i_z_c_i = phi_aggr
        .iter()
        .zip(c.iter())
        .map(|(phi_i, c_i)| inner_product_polynomial_ring_vector(phi_i, &z) * c_i)
        .fold(zero_poly(), |acc, val| acc + val);
    // calculate sum(h_ij * c_i * c_j)
    let mut sum_h_ij_c_i_c_j = zero_poly();
    for i in 0..size_r.value() {
        for j in 0..size_r.value() {
            sum_h_ij_c_i_c_j = sum_h_ij_c_i_c_j + (&h[i][j] * &c[i] * &c[j]);
        }
    }
    assert_eq!(sum_phi_i_z_c_i * Zq::from(2), sum_h_ij_c_i_c_j);

    println!("Verifier: Compute aggregated relation");
    // 7. check if sum(a_ij * g_ij) + sum(h_ii) -b ?= 0
    check_aggr_relation(&a_aggr, &b_aggr, &g, &h);

    println!("Verifier: Check opening of outer commitments(todo)");
    // 8. check if u1 is valid
    let u1_check = calculate_outer_comm_u1(
        b_matrix,
        c_matrix,
        &g_matrix_aggregated,
        &all_t_i_basis_form_aggregated,
        kappa1,
        t1,
        t2,
        size_r,
        size_n,
    );
    assert_eq!(u1, u1_check);
    // 9. check if u2 is valid
    let u2_check = calculate_outer_comm_u2(
        d_matrix,
        &h_gar_poly_basis_form_aggregated,
        t2,
        kappa2,
        size_r,
        size_n,
        deg_bound_d,
    );
    assert_eq!(u2, u2_check);
}
