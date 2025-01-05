use crate::core::{
    aggregation::*, conjugation_automorphism::conjugation_automorphism, constraint::*,
    decompose::*, gaussian_generator::generate_gaussian_distribution, norm::*, outer_commitment::*,
};
use crate::statement::Statement;
use crate::transcript::Transcript;
use algebra::{
    polynomial_ring::PolynomialRing,
    rq_matrix::RqMatrix,
    utils::{
        generate_random_polynomial_ring, inner_product_poly_matrix_and_poly_vector,
        inner_product_polynomial_ring_vector, inner_product_zq_vector,
        matrix_poly_times_poly_vector, zero_poly,
    },
    zq::Zq,
};
use profiler_macro::time_profiler;
use rand::Rng;

#[time_profiler]
pub fn prove(
    a_matrix: &RqMatrix,
    b_matrix: &[Vec<RqMatrix>],
    c_matrix: &[Vec<Vec<RqMatrix>>],
    d_matrix: &[Vec<Vec<RqMatrix>>],
) -> (Statement, Transcript) {
    // s is a vector of size r. each s_i is a PolynomialRing<Zq> with n coefficients
    let size_r = Zq::new(3); // r: Number of witness elements
    let size_n = Zq::new(5); // n
    let basis = Zq::new(10);
    let digits = Zq::new(3); // t1
    let t1 = digits;
    let t2 = digits;
    let kappa = Zq::new(3); // Example size
    let kappa1 = Zq::from(5);
    let kappa2 = Zq::from(5);
    let lambda = Zq::new(128);
    let double_lambda = lambda * Zq::new(2);
    let log_q = Zq::new(32);
    let deg_bound_d = Zq::new(8); // random polynomial degree bound
    let beta = Zq::new(70); // Example value for beta
    let mut rng = rand::thread_rng();
    let constraint_num_k = Zq::new(5);
    let constraint_num_l = Zq::new(5); // Define L
    println!("Prover: Generate random witness");
    let witness_s: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|_| {
            (0..size_n.value())
                .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                .collect()
        })
        .collect();
    let sum_squared_norms = poly_matrix_norm_squared(&witness_s);
    println!("Prover: Check witness norm is in range");
    // Check the condition
    assert!(
        sum_squared_norms <= beta.pow(2),
        "The condition is not satisfied: sum of squared norms exceeds beta^2"
    );

    println!("Prover: Generate random DPCS (dot product constraint system) constraints");
    // In DPCS (dot product constraint system), there are k constraints, each constraint has a, phi, and b
    // Generate random a^(k)_{i,j}: k length vector of matrices, each matrix is r x r, and each element is a Zq
    // TODO: a_ij == a_ji, and aij = 0 for |i − j| > 1. Furthermore, a^(k)_{ij} = 0 unless i,j ≤ 2ν.
    // refer to paper page 15
    let a_constraint: Vec<Vec<Vec<PolynomialRing>>> = (0..constraint_num_k.value())
        .map(|_| {
            (0..size_r.value())
                .map(|_| {
                    (0..size_n.value())
                        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                        .collect()
                })
                .collect()
        })
        .collect();
    let phi_constraint: Vec<Vec<Vec<PolynomialRing>>> = (0..constraint_num_k.value())
        .map(|_| {
            (0..size_r.value())
                .map(|_| {
                    (0..size_n.value())
                        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                        .collect()
                })
                .collect()
        })
        .collect();

    let b_constraint: Vec<PolynomialRing> = (0..constraint_num_k.value())
        .map(|k| calculate_b_constraint(&witness_s, &a_constraint[k], &phi_constraint[k]))
        .collect();

    // In DPCS(dot product constraint system) for constant terms(ct), there are k constraints, each constraint has a, phi and b.
    // Generate random a^(l)_{i,j}: l length vector of matrix, matrix length is r x r, each element is a Zq
    // todo: aij == aji, refer to paper page 10
    let a_constraint_ct: Vec<Vec<Vec<PolynomialRing>>> = (0..constraint_num_l.value())
        .map(|_| {
            (0..size_r.value())
                .map(|_| {
                    (0..size_n.value())
                        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                        .collect()
                })
                .collect()
        })
        .collect();

    // Generate random phi^(k)_{i}: k length vector of matrix, matrix length is r x n, each element in matrix is a Zq
    let phi_constraint_ct: Vec<Vec<Vec<PolynomialRing>>> = (0..constraint_num_l.value())
        .map(|_| {
            (0..size_r.value())
                .map(|_| {
                    (0..size_n.value())
                        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                        .collect()
                })
                .collect()
        })
        .collect();
    assert_eq!(phi_constraint_ct.len(), constraint_num_l.value());
    assert_eq!(phi_constraint_ct[0].len(), size_r.value());
    assert_eq!(phi_constraint_ct[0][0].len(), size_n.value());

    // calculate b^(l)
    let b_constraint_poly: Vec<PolynomialRing> = (0..constraint_num_l.value())
        .map(|l| calculate_b_constraint(&witness_s, &a_constraint_ct[l], &phi_constraint_ct[l]))
        .collect();
    // only keep constant term
    let b_constraint_ct: Vec<Zq> = (0..constraint_num_l.value())
        .map(|l| b_constraint_poly[l].coefficients[0])
        .collect();

    println!("Prover: Do Ajtai commitment");
    // A: matrix size: kappa * n, each element is PolynomialRing(R_q)
    // calculate t_i = A * s_i for all i = 1..r
    // size of t_i = (kappa * n)R_q * 1R_q = kappa * n
    let t: Vec<Vec<PolynomialRing>> = witness_s
        .iter()
        .map(|s_i| matrix_poly_times_poly_vector(&a_matrix.values, s_i))
        .collect();
    assert!(t.len() == kappa.value());

    // ================================================
    // prover
    // 2. GOAL: calculate ajtai commitment (1st outer commitment)
    // 2.1 split t to t_i for all i
    // 2.1.1 get basis b1, refer to paper page 16, labrador c code line 142
    // s = β / sqrt(rnd)
    // r: s_len, n
    // for example:
    // t_0: [
    //  [8, 46, 61, 71, 33, 33, 18], -> t[0][0]
    //  [20, 54, 94, 93, 70, 33, 14], -> t[0][1]
    //  [24, 40, 100, 85, 121, 57, 56],  -> t[0][2]
    //  [14, 37, 91, 118, 159, 109, 72],  -> t[0][3]
    // ]
    // such as basis = 10
    // t1: length of t[i][j][k], such as length of t[0][0][0] = length of [8, 0, 0] = 3
    // Then:
    // t_0_basis_form: [
    // [[8, 0, 0], [6, 4, 0], [1, 6, 0], [1, 7, 0], [3, 3, 0], [3, 3, 0], [8, 1, 0]]
    // [[0, 2, 0], [4, 5, 0], [4, 9, 0], [3, 9, 0], [0, 7, 0], [3, 3, 0], [4, 1, 0]]
    // [[4, 2, 0], [0, 4, 0], [0, 0, 1], [5, 8, 0], [1, 2, 1], [7, 5, 0], [6, 5, 0]]
    // [[4, 1, 0], [7, 3, 0], [1, 9, 0], [8, 1, 1], [9, 5, 1], [9, 0, 1], [2, 7, 0]]
    // ]
    println!("Prover: Do decomposition");

    let all_t_i_basis_form_aggregated = poly_matrix_decompose_and_aggregate(&t, basis, digits);

    // 2.2.1 get basis b2 same as 2.1.1
    // Calculate garbage polynomial g_ij = <s_i, s_j>
    let g: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|i| {
            (0..size_r.value())
                .map(|j| {
                    let s_i = &witness_s[i];
                    let s_j = &witness_s[j];
                    inner_product_polynomial_ring_vector(s_i, s_j)
                })
                .collect::<Vec<PolynomialRing>>()
        })
        .collect();
    assert_eq!(g.len(), size_r.value());
    assert_eq!(g[0].len(), size_r.value());
    for i in 0..size_r.value() {
        for j in 0..size_r.value() {
            assert_eq!(
                g[i][j], g[j][i],
                "g_ij is not equal to g_ji at indices ({}, {})",
                i, j
            );
        }
    }

    let g_matrix_aggregated = poly_matrix_decompose_and_aggregate(&g, basis, t2);

    // 2.3 calculate u1
    let u1 = calculate_outer_comm_u1(
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
    println!("Prover: Send proof u1");

    // ================================================

    println!("Prover: Do JL projection");
    // 3. GOAL: JL projection
    let nd = size_n * deg_bound_d;
    // generate gaussian distribution matrices
    // there are size_r matrices, each matrix size is 256 * nd
    // TODO: should from verifier
    let pai = (0..size_r.value())
        .map(|_| generate_gaussian_distribution(nd))
        .collect::<Vec<Vec<Vec<Zq>>>>();

    assert_eq!(pai.len(), size_r.value());
    assert_eq!(pai[0].len(), double_lambda.value());
    assert_eq!(pai[0][0].len(), nd.value());
    // 3.1 PI_i is randomly chosen from \Chi { -1, 0, 1 }^{256 * nd}
    //      (Using Guassian Distribution)

    // 3.2 caculate p_j = sum(<pi_i^(j), s_i>) for all i-r
    /*
        - pi^(j) is the j-th row of gaussian_distribution_matrix, j = 1..256
        - concat s_i's coefficients, output a vector with length nd
        - <pi_i^(j), s_i> = pi_i^(j)[0] * s_i[0] + pi_i^(j)[1] * s_i[1] + ... + pi_i^(j)[nd-1] * s_i[nd-1], is the inner product of pi_i^(j) and s_i
        - p_j = sum(<pi_i^(j), s_i>) for all i = 1..r, is a Zq
        - vector p = [p_1, p_2, ..., p_256], is a vector with length 256(2λ), type Vec<Zq>
    */

    // concat all s_i's coefficients into a vector
    // such as: s = [s_0, s_1]
    // s_0: [PolynomialRing{coefficients: [1, 2, 3]}, PolynomialRing{coefficients: [4, 5, 6]}, PolynomialRing{coefficients: [7, 8, 9]}], output: ss_0 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    // s_1: [PolynomialRing{coefficients: [1, 2, 3]}, PolynomialRing{coefficients: [4, 5, 6]}, PolynomialRing{coefficients: [7, 8, 9]}], output: ss_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    // so ss = [ss_0, ss_1]
    let s_coeffs: Vec<Vec<Zq>> = witness_s
        .iter()
        .map(|s_i| {
            s_i.iter()
                .flat_map(|s_i_poly| s_i_poly.coefficients.clone())
                .collect()
        })
        .collect();
    assert_eq!(s_coeffs.len(), size_r.value());
    assert_eq!(s_coeffs[0].len(), nd.value());
    // implement p calculation, inner product of gaussian_distribution_matrices and s_coeffs
    let mut p: Vec<Zq> = Vec::with_capacity(double_lambda.value());
    for j in 0..double_lambda.value() {
        let mut sum = Zq::new(0);
        for i in 0..size_r.value() {
            let pai_element = &pai[i][j];
            let s_i = &s_coeffs[i];
            let inner_product = inner_product_zq_vector(pai_element, s_i);
            sum += inner_product;
        }
        p.push(sum);
    }

    assert_eq!(p.len(), double_lambda.value());

    // sanity check: verify p_j = ct(sum(<σ−1(pi_i^(j)), s_i>)) for all i = 1..r
    for (j, &p_j) in p.iter().enumerate() {
        let mut sum = PolynomialRing {
            coefficients: vec![Zq::from(0); deg_bound_d.value()],
        };
        for i in 0..size_r.value() {
            let pai_element = &pai[i][j];
            let s_i = &s_coeffs[i];
            let pai_poly = PolynomialRing {
                coefficients: pai_element.clone(),
            };
            let pai_poly_ca = conjugation_automorphism(&pai_poly);
            let s_i_poly = PolynomialRing {
                coefficients: s_i.clone(),
            };
            sum = sum + &pai_poly_ca * s_i_poly;
        }
        assert_eq!(sum.coefficients[0], p_j);
    }
    println!("Prover: Send proof p");

    // todo: send p to verifier(put in transcript)

    // 3.3 Verifier have to check: || p || <= \sqrt{128} * beta

    // ================================================

    println!("Prover: Do aggregation");
    // 4. GOAL: Aggregation
    // 4.1 psi^(k) is randomly chosen from Z_q^{L}
    // k = 1..λ/log2^q
    let size_k = lambda.ceil_div(log_q);
    let psi: Vec<Vec<Zq>> = (0..size_k.value())
        .map(|_| {
            (0..constraint_num_l.value())
                .map(|_| Zq::new(rng.gen_range(1..10)))
                .collect()
        })
        .collect();
    assert_eq!(psi.len(), size_k.value());
    assert_eq!(psi[0].len(), constraint_num_l.value());

    // 4.2 omega^(k) is randomly chosen from Z_q^{256}
    //      (Both using Guassian Distribution)
    let omega: Vec<Vec<Zq>> = (0..size_k.value())
        .map(|_| {
            (0..double_lambda.value())
                .map(|_| Zq::new(rng.gen_range(1..10)))
                .collect()
        })
        .collect();
    assert_eq!(omega.len(), size_k.value());
    assert_eq!(omega[0].len(), double_lambda.value());

    // 4.3 caculate b^{''(k)}
    // 4.3.1 calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
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
    // 4.3.2 calculate phi_i^{''(k)} =
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

    // 4.3.3 calculate b^{''(k)} = sum(a_ij^{''(k)} * <s_i, s_j>) + sum(<phi_i^{''(k)}, s_i>)
    let b_ct_aggr = compute_aggr_ct_constraint_b(
        &a_ct_aggr,
        &phi_ct_aggr,
        size_k,
        size_r,
        deg_bound_d,
        &witness_s,
    );
    assert_eq!(b_ct_aggr.len(), size_k.value());

    // todo: send b^{''(k)} to verifier

    // Verifier check: b_0^{''(k)} ?= <⟨omega^(k),p⟩> + sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L
    for k in 0..size_k.value() {
        let b_k_0_from_poly: Zq = b_ct_aggr[k].coefficients[0];
        // sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L
        let mut b_k_0_computed: Zq = (0..constraint_num_l.value())
            .map(|l| {
                let psi_k_l = psi[k][l];
                let b_l_0 = b_constraint_ct[l];
                psi_k_l * b_l_0
            })
            .sum();
        // <⟨omega^(k),p⟩>
        let omega_k = &omega[k];
        let inner_product_omega_k_p = inner_product_zq_vector(omega_k, &p);
        // add them together
        b_k_0_computed += inner_product_omega_k_p;
        assert_eq!(b_k_0_from_poly, b_k_0_computed);
    }

    // ================================================

    println!("Prover: Aggregate linear constraints");
    // 5. GOAL: Calculate u2 (2nd outer commitment)
    // 5.1 vec<alpha> and vec<beta> are randomly chosen from R_q^{K} and R_q^{128/logQ}
    let alpha: Vec<PolynomialRing> = (0..constraint_num_k.value())
        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
        .collect();
    let beta: Vec<PolynomialRing> = (0..size_k.value())
        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
        .collect();
    // 5.2 phi_i = sum(alpha_k * phi_i) + beta_k * phi_i^{''(k)}
    let phi_aggr: Vec<Vec<PolynomialRing>> = compute_aggr_constraint_phi(
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
    assert_eq!(phi_aggr.len(), size_r.value());
    assert_eq!(phi_aggr[0].len(), size_n.value());

    // 5.3 Calculate garbage polynomial h_ij = 1/2 * (<phi_i, s_j> + <phi_j, s_i>)
    let h: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|i| {
            (0..size_r.value())
                .map(|j| {
                    let phi_i = &phi_aggr[i];
                    let phi_j = &phi_aggr[j];
                    let s_i = &witness_s[i];
                    let s_j = &witness_s[j];
                    inner_product_polynomial_ring_vector(phi_i, s_j)
                        + inner_product_polynomial_ring_vector(phi_j, s_i)
                    // Notice we do not divide by 2 here as paper described, because there is no division in the ring, we multiply by 2 instead with other terms to make verifier check work
                })
                .collect::<Vec<PolynomialRing>>()
        })
        .collect();
    assert_eq!(h.len(), size_r.value());
    assert_eq!(h[0].len(), size_r.value());
    for i in 0..size_r.value() {
        for j in (i + 1)..size_r.value() {
            assert_eq!(
                h[i][j], h[j][i],
                "h_ij is not equal to h_ji at indices ({}, {})",
                i, j
            );
        }
    }

    let h_gar_poly_basis_form_aggregated = poly_matrix_decompose_and_aggregate(&h, basis, digits);

    // 5.4 u2 = sum D_ij * h_ij^(k) for all k = 1..(t1-1)
    let u2 = calculate_outer_comm_u2(
        d_matrix,
        &h_gar_poly_basis_form_aggregated,
        t2,
        kappa2,
        size_r,
        size_n,
        deg_bound_d,
    );
    println!("Prover: Send proof u2");

    // Send u2 to verifier
    // transcript.add(u2)

    // ================================================

    println!("Prover: Amortize proof");
    // 6. GOAL: calculate z (Amortized Opening)
    // 6.1 c_i is randomly chosen from C, i = 1..r
    // todo: get c from challenge space, refer to paper page 6
    let c: Vec<PolynomialRing> = (0..size_r.value())
        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
        .collect();
    // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
    let z: Vec<PolynomialRing> = inner_product_poly_matrix_and_poly_vector(&witness_s, &c);
    // Send z, t_i, g_ij, h_ij to verifier
    // transcript.add(z);
    // return transcript;

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
            let h_ij = &h[i][j]; // Borrow h[i][j] instead of moving it
            let c_i = &c[i];
            let c_j = &c[j];
            sum_h_ij_c_i_c_j = sum_h_ij_c_i_c_j + (h_ij * c_i * c_j);
        }
    }

    assert_eq!(sum_phi_i_z_c_i * Zq::from(2), sum_h_ij_c_i_c_j);
    println!("Prover: Send amortized proof");
    println!("Prover is finished");
    let st = Statement {
        a_constraint,
        phi_constraint,
        b_constraint,
        a_constraint_ct,
        phi_constraint_ct,
        b_constraint_ct,
    };
    let tr = Transcript {
        u1,
        pai,
        p,
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
    };
    (st, tr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup::setup;
    use crate::verifier::verify;
    use std::vec;

    #[test]
    fn test_setup_prover() {
        let size_r = Zq::new(3); // r: Number of witness elements
        let size_n = Zq::new(5); // n
        let digits = Zq::new(3); // t1
        let t1 = digits;
        let t2 = digits;
        let kappa = Zq::new(3); // Example size
        let kappa1 = Zq::from(5);
        let kappa2 = Zq::from(5);

        // 0. setup
        // matrices A, B, C, D are common reference string
        let (a_matrix, b_matrix, c_matrix, d_matrix) =
            setup(size_n, size_r, t1, t2, kappa, kappa1, kappa2);
        // todo: st should be publicly shared
        let (st, tr) = prove(&a_matrix, &b_matrix, &c_matrix, &d_matrix);
        verify(st, tr, &a_matrix, &b_matrix, &c_matrix, &d_matrix);
    }

    #[test]
    fn test_h_verify() {
        let size_r = Zq::new(3); // r: Number of witness elements
        let size_n = Zq::new(5); // n
        let deg_bound_d = Zq::new(8); // random polynomial degree bound
                                      // generate size_r * size_n witness_s
        let witness_s: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|_| {
                (0..size_n.value())
                    .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                    .collect()
            })
            .collect();

        // generate size_r * size_n phi_aggr
        let phi_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|_| {
                (0..size_n.value())
                    .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        let h: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|i| {
                (0..size_r.value())
                    .map(|j| {
                        let phi_i = &phi_aggr[i];
                        let phi_j = &phi_aggr[j];
                        let s_i = &witness_s[i];
                        let s_j = &witness_s[j];
                        inner_product_polynomial_ring_vector(phi_i, s_j)
                            + inner_product_polynomial_ring_vector(phi_j, s_i)
                    })
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        let c: Vec<PolynomialRing> = (0..size_r.value())
            .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
            .collect();
        // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
        let z: Vec<PolynomialRing> = inner_product_poly_matrix_and_poly_vector(&witness_s, &c);

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
    }

    #[test]
    fn test_inner_product_and_z_computation() {
        let deg_bound_d = Zq::new(8); // random polynomial degree bound
                                      // this test it to check: <z, z> ?= sum(g_ij * c_i * c_j)
        let witness_s = vec![
            vec![
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
            ],
            vec![
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
            ],
            vec![
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
                generate_random_polynomial_ring(deg_bound_d.value()),
            ],
        ];

        let size_r = Zq::from(witness_s.len());
        let size_n = Zq::from(witness_s[0].len());
        let g: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|i| {
                (0..size_r.value())
                    .map(|j| {
                        let s_i = &witness_s[i];
                        let s_j = &witness_s[j];
                        assert_eq!(
                            s_i.len(),
                            size_n.value(),
                            "s_i must have the same length as size_n"
                        );
                        assert_eq!(
                            s_j.len(),
                            size_n.value(),
                            "s_j must have the same length as size_n"
                        );
                        inner_product_polynomial_ring_vector(s_i, s_j)
                    })
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        let c: Vec<PolynomialRing> = (0..size_r.value())
            .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
            .collect();
        // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
        let z: Vec<PolynomialRing> = inner_product_poly_matrix_and_poly_vector(&witness_s, &c);

        // check if <z, z> ?= sum(g_ij * c_i * c_j)
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
    }

    // TODO(junochiu): check b_k?
    /*#[test]
    fn test_calculate_b_k() {
        let r = 3;
        let n = 4;

        let s: Vec<Vec<PolynomialRing>> = (0..r)
            .map(|_| {
                (0..n)
                    .map(|n_i| PolynomialRing {
                        coefficients: vec![Zq::from(n_i)],
                    })
                    .collect()
            })
            .collect();

        let a_constraint: Vec<Vec<PolynomialRing>> = (0..r)
            .map(|_| {
                (0..r)
                    .map(|r_i| PolynomialRing {
                        coefficients: vec![Zq::from(r_i)],
                    })
                    .collect()
            })
            .collect();
        let phi_constraint: Vec<Vec<PolynomialRing>> = (0..r)
            .map(|_| {
                (0..n)
                    .map(|n_i| PolynomialRing {
                        coefficients: vec![Zq::from(n_i)],
                    })
                    .collect()
            })
            .collect();
        let b_k = calculate_b_constraint(&s, &a_constraint, &phi_constraint);
        //assert_eq!(b_k, 1983);
    }*/

    #[test]
    fn test_a_new() {
        let size = Zq::from(3);
        let a = RqMatrix::new(size, size);
        assert_eq!(a.values.len(), size.value());
        assert_eq!(a.values[0].len(), size.value());
    }

    #[test]
    fn test_calculate_a_times_s_i() {
        let size = Zq::from(2);
        let a = RqMatrix::new(size, size);
        let s_i = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(3), Zq::from(4)],
            },
        ];
        // matrix a: [[(2+2x), (2+2x)], [(2+2x), (2+2x)]]
        // s_i: [(1+2x), (3+4x)]
        // a * s_i = [
        //   [(2+2x), (2+2x)] * [(1+2x), (3+4x)],
        //   [(2+2x), (2+2x)] * [(1+2x), (3+4x)]
        // ]
        // = [
        //   (2+2x) * (1+2x) + (2+2x) * (3+4x),
        //   (2+2x) * (1+2x) + (2+2x) * (3+4x)
        // ]
        // = [
        //   (2+4x+2x+4x^2) + (6+8x+6x+8x^2),
        //   (2+4x+2x+4x^2) + (6+8x+6x+8x^2)
        // ]
        // = [
        //   8+20x+12x^2,
        //   8+20x+12x^2
        // ]
        let result = matrix_poly_times_poly_vector(&a.values, &s_i);
        let expected = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(8), Zq::from(20), Zq::from(12)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(8), Zq::from(20), Zq::from(12)],
            },
        ];
        assert_eq!(result, expected);
    }

    // add test for polynomial addition and multiplication with overload
    #[test]
    fn test_polynomial_addition_and_multiplication() {
        let a = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let b = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let c = &a + &b;
        assert_eq!(c.coefficients, vec![Zq::from(5), Zq::from(7), Zq::from(9)]);
        let d = &a * &b;
        assert_eq!(
            d.coefficients,
            vec![
                Zq::from(4),
                Zq::from(13),
                Zq::from(28),
                Zq::from(27),
                Zq::from(18)
            ]
        );
    }
}
