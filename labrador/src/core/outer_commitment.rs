use algebra::{PolynomialRing, RqMatrix, Zq};

// 2.3 calculate u1
// 2.3.1 B & C is randomly chosen similar to A
// 2.3.2 calculate u1 = sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))
// B_ik: Rq^{kappa1 x kappa}, t_i: Rq^{kappa}, t_i^(k): Rq^{kappa}
// B_ik * t_i^(k): Rq^{kappa1}
// First summation: ∑ B_ik * t_i^(k), 1 ≤ i ≤ r, 0 ≤ k ≤ t1−1
// Initialize u1 with zeros with size kappa1, each element is a polynomial ring
// TODO(junochiu): to check if it make sense to solve the too many argument issue
#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)] // triple-indexed math c_matrix[i][j][k] · g_ij^(k); iterators obscure the formula
pub fn calculate_outer_comm_u1(
    b_matrix: &[Vec<RqMatrix>],
    c_matrix: &[Vec<Vec<RqMatrix>>],
    g_matrix_aggregated: &[Vec<Vec<PolynomialRing>>],
    all_t_i_basis_form_aggregated: &[Vec<Vec<PolynomialRing>>],
    kappa1: Zq,
    t1: Zq,
    t2: Zq,
    size_r: Zq,
    size_n: Zq,
) -> Vec<PolynomialRing> {
    let mut u1 = vec![
        PolynomialRing {
            coefficients: vec![Zq::from(0); size_n.value()]
        };
        kappa1.value()
    ];
    // Calculate u1 using the pre-generated b_matrix
    for i in 0..size_r.value() {
        for k in 0..t1.value() {
            let b_i_k = &b_matrix[i][k];
            let t_i_k = &all_t_i_basis_form_aggregated[i][k];
            // matrix<Rq> * vector<Rq> -> vector<Rq>
            let b_ik_times_t_ik = b_i_k
                .values
                .iter()
                .map(|row| {
                    row.iter().zip(t_i_k.iter()).map(|(b, t)| b * t).fold(
                        PolynomialRing {
                            coefficients: vec![Zq::from(0); size_n.value()],
                        },
                        |acc, val| acc + val,
                    )
                })
                .collect::<Vec<PolynomialRing>>();
            u1 = u1
                .iter()
                .zip(b_ik_times_t_ik.iter())
                .map(|(a, b)| a + b)
                .collect();
        }
    }

    // Second summation: ∑ C_ijk * g_ij^(k)
    // Calculate u1 using the pre-generated c_matrix
    for i in 0..size_r.value() {
        for j in i..size_r.value() {
            for k in 0..t2.value() {
                let c_i_j_k = &c_matrix[i][j][k];
                let g_i_j = &g_matrix_aggregated[i][j];
                let c_i_j_k_times_g_i_j = c_i_j_k
                    .values
                    .iter()
                    .map(|row| {
                        row.iter().zip(g_i_j.iter()).map(|(c, g)| c * g).fold(
                            PolynomialRing {
                                coefficients: vec![Zq::from(0); size_n.value()],
                            },
                            |acc, val| acc + val,
                        )
                    })
                    .collect::<Vec<PolynomialRing>>();
                u1 = u1
                    .iter()
                    .zip(c_i_j_k_times_g_i_j.iter())
                    .map(|(a, b)| a + b)
                    .collect();
            }
        }
    }

    u1
}

pub fn calculate_outer_comm_u2(
    d_matrix: &[Vec<Vec<RqMatrix>>],
    h_gar_poly_basis_form_aggregated: &[Vec<Vec<PolynomialRing>>],
    t2: Zq,
    kappa2: Zq,
    size_r: Zq,
    size_n: Zq,
    deg_bound_d: Zq,
) -> Vec<PolynomialRing> {
    (0..size_r.value())
        .flat_map(|i| {
            (i..size_r.value()).flat_map(move |j| (0..t2.value()).map(move |k| (i, j, k)))
        })
        .fold(
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(0); deg_bound_d.value()]
                };
                kappa2.value()
            ],
            |acc, (i, j, k)| {
                let d_i_j_k = &d_matrix[i][j][k];
                let h_i_j = &h_gar_poly_basis_form_aggregated[i][j];
                let d_i_j_k_times_h_i_j = d_i_j_k
                    .values
                    .iter()
                    .map(|row| {
                        row.iter().zip(h_i_j.iter()).map(|(c, h)| c * h).fold(
                            PolynomialRing {
                                coefficients: vec![Zq::from(0); size_n.value()],
                            },
                            |acc, val| acc + val,
                        )
                    })
                    .collect::<Vec<PolynomialRing>>();
                acc.iter()
                    .zip(d_i_j_k_times_h_i_j.iter())
                    .map(|(a, b)| a + b)
                    .collect()
            },
        )
}
