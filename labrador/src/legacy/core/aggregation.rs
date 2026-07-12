use crate::legacy::core::conjugation_automorphism::conjugation_automorphism;
use algebra::{inner_product_polynomial_ring_vector, zero_poly, PolynomialRing, Zq};

// 4.3.1 aggregation: calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
pub fn compute_aggr_ct_constraint_a(
    a_constraint_ct: &[Vec<Vec<PolynomialRing>>],
    psi: &[Vec<Zq>],
    size_k: Zq,
    size_r: Zq,
    constraint_num_l: Zq,
    deg_bound_d: Zq,
) -> Vec<Vec<Vec<PolynomialRing>>> {
    let a_ct_aggr: Vec<Vec<Vec<PolynomialRing>>> = (0..size_k.value())
        .map(|k| {
            let psi_k = &psi[k];
            (0..size_r.value())
                .map(|i| {
                    (0..size_r.value())
                        .map(|j| {
                            (0..constraint_num_l.value())
                                .map(|l| &a_constraint_ct[l][i][j] * psi_k[l])
                                .fold(
                                    PolynomialRing {
                                        coefficients: vec![Zq::from(0); deg_bound_d.value()],
                                    },
                                    |acc, x| acc + x,
                                )
                        })
                        .collect::<Vec<PolynomialRing>>()
                })
                .collect::<Vec<Vec<PolynomialRing>>>()
        })
        .collect();
    a_ct_aggr
}

// 4.3.2 aggregation: calculate phi_i^{''(k)} =
//       sum(psi_l^(k) * phi_i^{'(l)}) for all l = 1..L
//       + sum(omega_j^(k) * sigma_{-1} * pi_i^{j)) for all j = 1..256
// TODO(junochiu): check if it make sense to solve the too many argument issue
#[allow(clippy::too_many_arguments)]
pub fn compute_aggr_ct_constraint_phi(
    phi_constraint_ct: &[Vec<Vec<PolynomialRing>>],
    pai: &[Vec<Vec<Zq>>],
    size_k: Zq,
    size_r: Zq,
    constraint_num_l: Zq,
    deg_bound_d: Zq,
    size_n: Zq,
    double_lambda: Zq,
    psi: &[Vec<Zq>],
    omega: &[Vec<Zq>],
) -> Vec<Vec<Vec<PolynomialRing>>> {
    let phi_ct_aggr: Vec<Vec<Vec<PolynomialRing>>> = (0..size_k.value())
        .map(|k| {
            (0..size_r.value())
                .map(|i| {
                    // Part 1: sum(psi_l^(k) * phi_constraint_ct[l][i] for all l)
                    let part1: Vec<PolynomialRing> = (0..constraint_num_l.value())
                        .map(|l| {
                            let psi = psi[k][l];
                            phi_constraint_ct[l][i]
                                .iter()
                                .map(|p| p * psi)
                                .collect::<Vec<PolynomialRing>>()
                        })
                        .fold(
                            vec![
                                PolynomialRing {
                                    coefficients: vec![Zq::from(0); deg_bound_d.value()]
                                };
                                size_n.value()
                            ],
                            |acc, product| {
                                acc.iter().zip(product.iter()).map(|(a, b)| a + b).collect()
                            },
                        );

                    // Part 2: sum(omega_j^(k) * sigma_{-1} * pi_i^{j} for all j)
                    let part2: Vec<PolynomialRing> = (0..double_lambda.value())
                        .map(|j| {
                            let omega = omega[k][j];
                            pai[i][j]
                                .chunks(deg_bound_d.value())
                                .take(size_n.value())
                                .map(|chunk| {
                                    let pai_poly = PolynomialRing {
                                        coefficients: chunk.to_vec(),
                                    };
                                    let pai_poly_ca = conjugation_automorphism(&pai_poly);
                                    pai_poly_ca * omega
                                })
                                .collect::<Vec<PolynomialRing>>()
                        })
                        .fold(
                            vec![
                                PolynomialRing {
                                    coefficients: vec![Zq::from(0); 1]
                                };
                                size_n.value()
                            ],
                            |acc, chunks_ca| {
                                acc.iter()
                                    .zip(chunks_ca.iter())
                                    .map(|(a, b)| a + b)
                                    .collect()
                            },
                        );

                    // Sum part1 and part2 element-wise
                    part1
                        .iter()
                        .zip(part2.iter())
                        .map(|(a, b)| a + b)
                        .collect::<Vec<PolynomialRing>>()
                })
                .collect::<Vec<Vec<PolynomialRing>>>()
        })
        .collect();
    phi_ct_aggr
}

// 4.3.3 aggregation: calculate b^{''(k)} = sum(a_ij^{''(k)} * <s_i, s_j>) + sum(<phi_i^{''(k)}, s_i>)
pub fn compute_aggr_ct_constraint_b(
    a_ct_aggr: &[Vec<Vec<PolynomialRing>>],
    phi_ct_aggr: &[Vec<Vec<PolynomialRing>>],
    size_k: Zq,
    size_r: Zq,
    deg_bound_d: Zq,
    witness_s: &[Vec<PolynomialRing>],
) -> Vec<PolynomialRing> {
    (0..size_k.value())
        .map(|k| {
            (0..size_r.value())
                .map(|i| {
                    (0..size_r.value())
                        .map(|j| {
                            &a_ct_aggr[k][i][j]
                                * inner_product_polynomial_ring_vector(&witness_s[i], &witness_s[j])
                        })
                        .fold(
                            PolynomialRing {
                                coefficients: vec![Zq::from(0); deg_bound_d.value()],
                            },
                            |acc, x| acc + x,
                        )
                        + inner_product_polynomial_ring_vector(&phi_ct_aggr[k][i], &witness_s[i])
                })
                .fold(
                    PolynomialRing {
                        coefficients: vec![Zq::from(0); deg_bound_d.value()],
                    },
                    |acc, x| acc + x,
                )
        })
        .collect::<Vec<PolynomialRing>>()
}

// aggregation: a_i = sum(alpha_k * a_ij) + sum(beta_k * a_ij^{''(k)})
pub fn compute_aggr_constraint_a(
    a_constraint: &[Vec<Vec<PolynomialRing>>],
    a_ct_aggr: &[Vec<Vec<PolynomialRing>>],
    constraint_num_k: Zq,
    alpha: &[PolynomialRing],
    beta: &[PolynomialRing],
    size_r: Zq,
    size_k: Zq,
) -> Vec<Vec<PolynomialRing>> {
    (0..size_r.value())
        .map(|i| {
            (0..size_r.value())
                .map(|j| {
                    // Part 1: sum(alpha_k * a_ij)
                    let part1: PolynomialRing = (0..constraint_num_k.value())
                        .map(|k| &a_constraint[k][i][j] * &alpha[k])
                        .fold(zero_poly(), |acc, product| acc + product);

                    // Part 2: sum(beta_k * a_ij^{''(k)})
                    let part2: PolynomialRing = (0..size_k.value())
                        .map(|k| &a_ct_aggr[k][i][j] * &beta[k])
                        .fold(zero_poly(), |acc, product| acc + product);

                    // Sum part1 and part2 element-wise
                    part1 + part2
                })
                .collect::<Vec<PolynomialRing>>()
        })
        .collect::<Vec<Vec<PolynomialRing>>>()
}

// aggregation: phi_i = sum(alpha_k * phi_i) + sum(beta_k * phi_i^{''(k)})
// TODO(junochiu): check if it make sense to solve the too many argument issue
#[allow(clippy::too_many_arguments)]
pub fn compute_aggr_constraint_phi(
    phi_constraint: &[Vec<Vec<PolynomialRing>>],
    phi_ct_aggr: &[Vec<Vec<PolynomialRing>>],
    constraint_num_k: Zq,
    alpha: &[PolynomialRing],
    beta: &[PolynomialRing],
    size_r: Zq,
    size_n: Zq,
    deg_bound_d: Zq,
    size_k: Zq,
) -> Vec<Vec<PolynomialRing>> {
    let phi_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|i| {
            // Part 1: sum(alpha_k * phi_i)
            let part1: Vec<PolynomialRing> = (0..constraint_num_k.value())
                .map(|k| {
                    let alpha = &alpha[k];
                    phi_constraint[k][i]
                        .iter()
                        .map(|p| p * alpha)
                        .collect::<Vec<PolynomialRing>>()
                })
                .fold(
                    vec![
                        PolynomialRing {
                            coefficients: vec![Zq::from(0); deg_bound_d.value()]
                        };
                        size_n.value()
                    ],
                    |acc, product| acc.iter().zip(product.iter()).map(|(a, b)| a + b).collect(),
                );

            // Part 2: sum(beta_k * phi_i^{''(k)})
            let part2: Vec<PolynomialRing> = (0..size_k.value())
                .map(|k| {
                    let beta = &beta[k];
                    phi_ct_aggr[k][i]
                        .iter()
                        .map(|p| p * beta)
                        .collect::<Vec<PolynomialRing>>()
                })
                .fold(
                    vec![
                        PolynomialRing {
                            coefficients: vec![Zq::from(0); deg_bound_d.value()]
                        };
                        size_n.value()
                    ],
                    |acc, product| acc.iter().zip(product.iter()).map(|(a, b)| a + b).collect(),
                );
            // Sum part1 and part2 element-wise
            part1
                .iter()
                .zip(part2.iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<PolynomialRing>>()
        })
        .collect();
    phi_aggr
}

// aggregation: b_i = sum(alpha_k * b^(k)) + sum(beta_k * b^{''(k)})
pub fn compute_aggr_constraint_b(
    b_constraint: &[PolynomialRing],
    b_ct_aggr: &[PolynomialRing],
    constraint_num_k: Zq,
    alpha: &[PolynomialRing],
    beta: &[PolynomialRing],
    size_k: Zq,
) -> PolynomialRing {
    // Part 1: sum(alpha_k * b^(k))
    let part1: PolynomialRing = (0..constraint_num_k.value())
        .map(|k| &b_constraint[k] * &alpha[k])
        .fold(zero_poly(), |acc, product| acc + product);

    // Part 2: sum(beta_k * b^{''(k)})
    let part2: PolynomialRing = (0..size_k.value())
        .map(|k| &b_ct_aggr[k] * &beta[k])
        .fold(zero_poly(), |acc, product| acc + product);

    // Sum part1 and part2
    part1 + part2
}

pub fn check_aggr_relation(
    a_aggr: &[Vec<PolynomialRing>],
    b_aggr: &PolynomialRing,
    g: &[Vec<PolynomialRing>],
    h: &[Vec<PolynomialRing>],
) {
    let size_r = Zq::from(a_aggr.len());
    // 7. check if sum(a_ij * g_ij) + sum(h_ii) -b ?= 0
    // 7.1 calculate sum(a_ij * g_ij)
    let sum_a_ij_g_ij = a_aggr
        .iter()
        .zip(g.iter())
        .map(|(a_i, g_i)| {
            a_i.iter()
                .zip(g_i.iter())
                .map(|(a_ij, g_ij)| a_ij * g_ij)
                .fold(zero_poly(), |acc, val| acc + val)
        })
        .fold(zero_poly(), |acc, val| acc + val);

    // 7.2 calculate sum(h_ii)
    let sum_h_ii = (0..size_r.value()).fold(zero_poly(), |acc, i| acc + &h[i][i]);

    // 2 times sum
    let b_aggr2 = b_aggr * Zq::from(2);
    let sum_a_ij_g_ij2 = sum_a_ij_g_ij * Zq::from(2);

    // 7.3 check if sum(a_ij * g_ij) + sum(h_ii) -b ?= 0
    assert_eq!(sum_a_ij_g_ij2 + sum_h_ii, b_aggr2);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::legacy::core::{
        constraint::calculate_b_constraint, gaussian_generator::generate_gaussian_distribution,
    };
    use algebra::utils::generate_random_polynomial_ring;
    use rand::Rng;
    #[test]
    fn test_aggr_relation_full_example() {
        let size_r = Zq::new(3); // r: Number of witness elements
        let size_n = Zq::new(5); // n
        let deg_bound_d = Zq::new(8); // random polynomial degree bound
        let lambda = Zq::new(128);
        let double_lambda = lambda * Zq::new(2);
        let constraint_num_l = Zq::new(5);
        let constraint_num_k = Zq::new(5);
        let log_q = Zq::new(2);
        let mut rng = rand::thread_rng();
        // generate size_r * size_n witness_s
        let witness_s: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|_| {
                (0..size_n.value())
                    .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                    .collect()
            })
            .collect();

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

        println!("Prover: Do JL projection");
        // 3. GOAL: JL projection
        let nd = size_n * deg_bound_d;
        // generate gaussian distribution matrices
        // there are size_r matrices, each matrix size is 256 * nd
        let pai = (0..size_r.value())
            .map(|_| generate_gaussian_distribution(nd))
            .collect::<Vec<Vec<Vec<Zq>>>>();

        let size_k = Zq::new(lambda.value() / log_q.value());
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

        let alpha: Vec<PolynomialRing> = (0..constraint_num_k.value())
            .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
            .collect();
        let beta: Vec<PolynomialRing> = (0..size_k.value())
            .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
            .collect();

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

        check_aggr_relation(&a_aggr, &b_aggr, &g, &h);
    }

    #[test]
    fn test_aggr_relation() {
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

        // generate size_r * size_r a_aggr
        let a_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|_| {
                (0..size_r.value())
                    .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        // generate size_r * size_n phi_aggr
        let phi_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|_| {
                (0..size_n.value())
                    .map(|_| {
                        // generate_random_polynomial_ring(deg_bound_d.value())
                        generate_random_polynomial_ring(64)
                    })
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        let b_aggr: PolynomialRing = calculate_b_constraint(&witness_s, &a_aggr, &phi_aggr);

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

        // Calculate b^(k)
        let mut quad_sum = zero_poly();
        let mut linear_sum = zero_poly();
        for i in 0..size_r.value() {
            for j in 0..size_r.value() {
                // calculate inner product of s[i] and s[j], will return a single PolynomialRing
                let elem_s_i = &witness_s[i];
                let elem_s_j = &witness_s[j];
                // Calculate inner product and update b
                let inner_product_si_sj = inner_product_polynomial_ring_vector(elem_s_i, elem_s_j);
                let a_constr = &a_aggr[i][j];
                quad_sum = quad_sum + (inner_product_si_sj * a_constr);
            }
            // calculate inner product of s[i] and phi
            let inner_product_si_phi =
                inner_product_polynomial_ring_vector(&witness_s[i], &phi_aggr[i]);
            println!("inner_product_si_phi: {:?}", inner_product_si_phi);
            println!("h[i][i]: {:?}", h[i][i]);
            assert_eq!(&inner_product_si_phi * Zq::from(2), h[i][i]);
            linear_sum = linear_sum + inner_product_si_phi;
        }

        // use function to check
        check_aggr_relation(&a_aggr, &b_aggr, &g, &h);

        // ================================================
        // manually check

        // 7. check if sum(a_ij * g_ij) + sum(h_ii) -b ?= 0
        // 7.1 calculate sum(a_ij * g_ij)
        let sum_a_ij_g_ij = a_aggr
            .iter()
            .zip(g.iter())
            .map(|(a_i, g_i)| {
                a_i.iter()
                    .zip(g_i.iter())
                    .map(|(a_ij, g_ij)| a_ij * g_ij)
                    .fold(zero_poly(), |acc, val| acc + val)
            })
            .fold(zero_poly(), |acc, val| acc + val);

        // 7.2 calculate sum(h_ii)
        let sum_h_ii = (0..size_r.value()).fold(zero_poly(), |acc, i| acc + &h[i][i]);

        // 2 times sum
        let quad_sum2 = quad_sum * Zq::from(2);
        let linear_sum2 = linear_sum * Zq::from(2);
        let b_aggr2 = &b_aggr * Zq::from(2);
        let sum_a_ij_g_ij2 = sum_a_ij_g_ij * Zq::from(2);
        assert_eq!(linear_sum2, sum_h_ii);
        assert_eq!(quad_sum2, sum_a_ij_g_ij2);
        assert_eq!(quad_sum2 + linear_sum2, b_aggr2);
        // 7.3 check if sum(a_ij * g_ij) + sum(h_ii) -b ?= 0
        assert_eq!(sum_a_ij_g_ij2 + sum_h_ii, b_aggr2);
    }
}
