use crate::algebra::{PolynomialRing, RqMatrix, Zq};
use crate::setup::setup;
use profiler_macro::time_profiler;
use rand::Rng;

// a: Vec<PolynomialRing>, b: PolynomialRing
// calculate c = a * b, c_i = a_i * b
// c: Vec<PolynomialRing>
fn poly_vec_times_poly(a: &Vec<PolynomialRing>, b: &PolynomialRing) -> Vec<PolynomialRing> {
    a.iter().map(|a_i| a_i * b).collect()
}
// a: Vec<PolynomialRing>, b: Vec<PolynomialRing>
// calculate c = a + b, c_i = a_i + b_i
// c: Vec<PolynomialRing>
fn poly_vec_add_poly_vec(a: &Vec<PolynomialRing>, b: &Vec<PolynomialRing>) -> Vec<PolynomialRing> {
    a.iter().zip(b.iter()).map(|(a_i, b_i)| a_i + b_i).collect()
}

// inner product of 2 vectors of PolynomialRing
fn inner_product_polynomial_ring_vector(
    a: &[PolynomialRing],
    b: &[PolynomialRing],
) -> PolynomialRing {
    assert_eq!(
        a.len(),
        b.len(),
        "inner_product_polynomial_ring_vector: a and b must have the same length"
    );
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| a * b)
        .collect::<Vec<PolynomialRing>>()
        .into_iter()
        .reduce(|acc, x| acc + x)
        .unwrap()
}

fn inner_product_zq_vector(a: &Vec<Zq>, b: &Vec<Zq>) -> Zq {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
}

// a: Vec<Vec<PolynomialRing>>, b: Vec<PolynomialRing>
// calculate c = sum(c_i), c_i = poly_vec_times_poly(a_i, b_i)
// c: Vec<PolynomialRing>
fn inner_product_poly_matrix_and_poly_vector(
    poly_matrix: &Vec<Vec<PolynomialRing>>,
    poly_vector: &Vec<PolynomialRing>,
) -> Vec<PolynomialRing> {
    assert_eq!(poly_matrix.len(), poly_vector.len(), "inner_product_poly_matrix_and_poly_vector: poly_matrix and poly_vector must have the same length");
    poly_matrix
        .iter()
        .zip(poly_vector.iter())
        .map(|(poly_matrix_row, poly_vector_element)| {
            poly_vec_times_poly(&poly_matrix_row, &poly_vector_element)
        })
        .fold(
            vec![zero_poly(); poly_matrix[0].len()],
            |acc, x: Vec<PolynomialRing>| poly_vec_add_poly_vec(&acc, &x),
        )
}
// Function to calculate b^(k)
fn calculate_b_constraint(
    s: &[Vec<PolynomialRing>],
    a_constraint: &[Vec<PolynomialRing>],
    phi_constraint: &[Vec<PolynomialRing>],
) -> PolynomialRing {
    let mut b: PolynomialRing = PolynomialRing {
        coefficients: vec![Zq::from(0)],
    };
    let s_len_usize = s.len();

    // Calculate b^(k)
    for i in 0..s_len_usize {
        for j in 0..s_len_usize {
            // calculate inner product of s[i] and s[j], will return a single PolynomialRing
            let elem_s_i = &s[i];
            let elem_s_j = &s[j];
            // Calculate inner product and update b
            let inner_product_si_sj = inner_product_polynomial_ring_vector(elem_s_i, elem_s_j);
            let a_constr = &a_constraint[i][j];
            b = b + (inner_product_si_sj * a_constr);
        }
        // calculate inner product of s[i] and phi
        let inner_product_si_phi = inner_product_polynomial_ring_vector(&s[i], &phi_constraint[i]);
        b = b + inner_product_si_phi;
    }

    b
}

// calculate matrix times vector of PolynomialRing
fn matrix_poly_times_poly_vector(
    poly_matrix: &Vec<Vec<PolynomialRing>>,
    poly_vec: &Vec<PolynomialRing>,
) -> Vec<PolynomialRing> {
    poly_matrix
        .iter()
        .map(|row| inner_product_polynomial_ring_vector(row, poly_vec))
        .collect::<Vec<PolynomialRing>>()
}

// convert number to basis
// 42 = 0 * 2^7 + 1 * 2^6 + 0 * 2^5 + 1 * 2^4 + 0 * 2^3 + 1 * 2^2 + 0 * 2^1 + 0 * 2^0
// first digit: 42 / 2 = 21, result_i = 0
// second digit: 21 / 2 = 10, result_i = 1
// third digit: 10 / 2 = 5, result_i = 0
// forth digit: 5 / 2 = 2, result_i = 1
// fifth digit: 2 / 2 = 1, result_i = 0
// sixth digit: 1 / 2 = 0, result_i = 1

fn num_to_basis(num: Zq, basis: Zq, digits: Zq) -> Vec<Zq> {
    let mut result = Vec::new();
    let mut remainder = num;

    let zero = Zq::from(0);
    let base = basis;

    for _ in 0..digits.value() {
        let digit = remainder % base;
        result.push(digit);
        remainder = Zq::from(remainder.value() / base.value());
    }

    while result.len() < digits.value() {
        // push 0 to the highest position
        result.push(zero);
    }

    result
}

// convert ring polynomial to basis
fn ring_polynomial_to_basis(poly: &PolynomialRing, basis: Zq, digits: Zq) -> Vec<Vec<Zq>> {
    poly.coefficients
        .iter()
        .map(|&coeff| num_to_basis(coeff, basis, digits))
        .collect()
}

fn generate_gaussian_distribution(nd: Zq) -> Vec<Vec<Zq>> {
    let nd_usize: usize = nd.value();
    let modulus: usize = Zq::modulus();
    let mut rng = rand::thread_rng();
    let mut matrix = vec![vec![Zq::from(0); nd_usize]; 256]; // Initialize a 256 x nd matrix

    for row in matrix.iter_mut() {
        for cell in row.iter_mut() {
            let random_value: f32 = rng.gen(); // Generate a random float between 0 and 1
            *cell = if random_value < 0.25 {
                // todo: should we use symmetric distribution from -q/2 to q/2?
                Zq::from(modulus - 1) // 1/4 probability
            } else if random_value < 0.75 {
                Zq::from(0) // 1/2 probability
            } else {
                Zq::from(1) // 1/4 probability
            };
        }
    }

    matrix
}

// Conjugation Automorphism σ_{-1}
// for polynomial ring a = 1+2x+3x^2, since x^64 = -1, apply this method to a, will get 1+2*(Zq.modulus()-1) * x^(64-1) +3*(Zq.modulus()-1) * x^(64-2)
//todo: Aut(Rq) ∼= Z×2d what is this???
fn conjugation_automorphism(poly: &PolynomialRing) -> PolynomialRing {
    let modulus_minus_one = Zq::from(Zq::modulus() - 1);
    let transformed_coeffs: Vec<Zq> = (0..PolynomialRing::DEGREE_BOUND)
        .map(|i| {
            if i < poly.coefficients.len() {
                if i == 0 {
                    poly.coefficients[i]
                } else {
                    poly.coefficients[i] * modulus_minus_one
                }
            } else {
                Zq::from(0)
            }
        })
        .collect();
    // reverse the coefficients except constant term
    let reversed_coefficients = transformed_coeffs
        .iter()
        .take(1)
        .cloned()
        .chain(transformed_coeffs.iter().skip(1).rev().cloned())
        .collect::<Vec<Zq>>();
    PolynomialRing {
        coefficients: reversed_coefficients,
    }
}

fn generate_random_polynomial_ring(deg_bound_d: usize) -> PolynomialRing {
    let mut rng = rand::thread_rng();
    PolynomialRing {
        coefficients: (0..deg_bound_d)
            .map(|_| Zq::from(rng.gen_range(1..5)))
            .collect(),
    }
}

// aggregate basis form of a vector of PolynomialRing
fn aggregate_poly_vec_basis_form(poly_basis_form: &Vec<Vec<Vec<Zq>>>) -> Vec<Vec<PolynomialRing>> {
    poly_basis_form
        .iter()
        .map(|poly_i_j_basis_form| {
            let num_loop_needed = poly_i_j_basis_form.first().map_or(0, |v| v.len());
            (0..num_loop_needed)
                .map(|k| {
                    let coefficients = poly_i_j_basis_form
                        .iter()
                        .map(|basis_part| basis_part.get(k).cloned().unwrap_or(Zq::from(0)))
                        .collect();
                    PolynomialRing { coefficients }
                })
                .collect()
        })
        .collect()
}

fn poly_vec_decompose_to_basis(
    poly: &Vec<PolynomialRing>,
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<Zq>>> {
    poly.iter()
        .map(|poly_i| ring_polynomial_to_basis(poly_i, basis, digits))
        .collect::<Vec<Vec<Vec<Zq>>>>()
}

fn poly_matrix_decompose_to_basis(
    poly: &Vec<Vec<PolynomialRing>>,
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<Vec<Zq>>>> {
    poly.iter()
        .map(|poly_i| poly_vec_decompose_to_basis(poly_i, basis, digits))
        .collect::<Vec<Vec<Vec<Vec<Zq>>>>>()
}

fn poly_vec_decompose_and_aggregate(
    poly: &Vec<PolynomialRing>,
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<PolynomialRing>> {
    // Decompose each PolynomialRing into basis parts
    let poly_basis_form = poly_vec_decompose_to_basis(poly, basis, digits);
    aggregate_poly_vec_basis_form(&poly_basis_form)
}

fn poly_matrix_decompose_and_aggregate(
    poly: &Vec<Vec<PolynomialRing>>,
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<PolynomialRing>>> {
    // Decompose h_ij into basis t_1 parts
    let poly_basis_form = poly_matrix_decompose_to_basis(poly, basis, digits);

    // Pick elements at each position across all inner vectors and aggregate them
    poly_basis_form
        .iter()
        .map(aggregate_poly_vec_basis_form)
        .collect()
}

// Calculate the sum of squared norms for a single PolynomialRing instance.
fn poly_norm_squared(poly: &PolynomialRing) -> Zq {
    poly.coefficients
        .iter()
        .fold(Zq::new(0), |acc, coeff| acc + coeff.pow(2))
}

// Calculate the sum of squared norms for a vector of PolynomialRing instances.
fn poly_vec_norm_squared(polys: &Vec<PolynomialRing>) -> Zq {
    polys
        .iter()
        .fold(Zq::new(0), |acc, poly| acc + poly_norm_squared(poly))
}

// Calculate the sum of squared norms for a matrix of PolynomialRing instances.
fn poly_matrix_norm_squared(poly_matrix: &Vec<Vec<PolynomialRing>>) -> Zq {
    poly_matrix.iter().fold(Zq::new(0), |acc, vector| {
        acc + poly_vec_norm_squared(vector)
    })
}

// Calculate the sum of squared norms for a 3D vector of PolynomialRing instances.
fn poly_3d_norm_squared(polymat3d: &Vec<Vec<Vec<PolynomialRing>>>) -> Zq {
    polymat3d.iter().fold(Zq::new(0), |acc, poly_matrix| {
        acc + poly_matrix_norm_squared(poly_matrix)
    })
}

// statement
struct St {
    a_constraint: Vec<Vec<Vec<PolynomialRing>>>,
    phi_constraint: Vec<Vec<Vec<PolynomialRing>>>,
    b_constraint: Vec<PolynomialRing>,
    a_constraint_ct: Vec<Vec<Vec<PolynomialRing>>>,
    phi_constraint_ct: Vec<Vec<Vec<PolynomialRing>>>,
    b_constraint_ct: Vec<Zq>,
}

struct Tr {
    u1: Vec<PolynomialRing>,        // Replace with the actual type
    pai: Vec<Vec<Vec<Zq>>>,         // Replace with the actual type
    p: Vec<Zq>,                     // Replace with the actual type
    psi: Vec<Vec<Zq>>,              // Replace with the actual type
    omega: Vec<Vec<Zq>>,            // Replace with the actual type
    b_ct_aggr: Vec<PolynomialRing>, // Replace with the actual type
    alpha: Vec<PolynomialRing>,     // Replace with the actual type
    beta: Vec<PolynomialRing>,      // Replace with the actual type
    u2: Vec<PolynomialRing>,
    c: Vec<PolynomialRing>,
    z: Vec<PolynomialRing>,
    t: Vec<Vec<PolynomialRing>>, // Replace with the actual type
    g: Vec<Vec<PolynomialRing>>, // Replace with the actual type
    h: Vec<Vec<PolynomialRing>>,
}

// 4.3.1 aggregation: calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
fn compute_aggr_ct_constraint_a(
    a_constraint_ct: &Vec<Vec<Vec<PolynomialRing>>>,
    psi: &Vec<Vec<Zq>>,
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
fn compute_aggr_ct_constraint_phi(
    phi_constraint_ct: &Vec<Vec<Vec<PolynomialRing>>>,
    pai: &Vec<Vec<Vec<Zq>>>,
    size_k: Zq,
    size_r: Zq,
    constraint_num_l: Zq,
    deg_bound_d: Zq,
    size_n: Zq,
    double_lambda: Zq,
    psi: &Vec<Vec<Zq>>,
    omega: &Vec<Vec<Zq>>,
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
fn compute_aggr_ct_constraint_b(
    a_ct_aggr: &Vec<Vec<Vec<PolynomialRing>>>,
    phi_ct_aggr: &Vec<Vec<Vec<PolynomialRing>>>,
    size_k: Zq,
    size_r: Zq,
    deg_bound_d: Zq,
    witness_s: &Vec<Vec<PolynomialRing>>,
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
fn compute_aggr_constraint_a(
    a_constraint: &Vec<Vec<Vec<PolynomialRing>>>,
    a_ct_aggr: &Vec<Vec<Vec<PolynomialRing>>>,
    constraint_num_k: Zq,
    alpha: &Vec<PolynomialRing>,
    beta: &Vec<PolynomialRing>,
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
fn compute_aggr_constraint_phi(
    phi_constraint: &Vec<Vec<Vec<PolynomialRing>>>,
    phi_ct_aggr: &Vec<Vec<Vec<PolynomialRing>>>,
    constraint_num_k: Zq,
    alpha: &Vec<PolynomialRing>,
    beta: &Vec<PolynomialRing>,
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
fn compute_aggr_constraint_b(
    b_constraint: &Vec<PolynomialRing>,
    b_ct_aggr: &Vec<PolynomialRing>,
    constraint_num_k: Zq,
    alpha: &Vec<PolynomialRing>,
    beta: &Vec<PolynomialRing>,
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

fn zero_poly() -> PolynomialRing {
    PolynomialRing {
        coefficients: vec![Zq::from(0); 1],
    }
}

fn check_aggr_relation(
    a_aggr: &Vec<Vec<PolynomialRing>>,
    b_aggr: &PolynomialRing,
    g: &Vec<Vec<PolynomialRing>>,
    h: &Vec<Vec<PolynomialRing>>,
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

// 2.3 calculate u1
// 2.3.1 B & C is randomly chosen similar to A
// 2.3.2 calculate u1 = sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))
// B_ik: Rq^{kappa1 x kappa}, t_i: Rq^{kappa}, t_i^(k): Rq^{kappa}
// B_ik * t_i^(k): Rq^{kappa1}
// First summation: ∑ B_ik * t_i^(k), 1 ≤ i ≤ r, 0 ≤ k ≤ t1−1
// Initialize u1 with zeros with size kappa1, each element is a polynomial ring
fn calculate_outer_comm_u1(
    b_matrix: &Vec<Vec<RqMatrix>>,
    c_matrix: &Vec<Vec<Vec<RqMatrix>>>,
    g_matrix_aggregated: &Vec<Vec<Vec<PolynomialRing>>>,
    all_t_i_basis_form_aggregated: &Vec<Vec<Vec<PolynomialRing>>>,
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

// calculate u2 = sum D_ij * h_ij^(k) for all k = 1..(t1-1)
fn calculate_outer_comm_u2(
    d_matrix: &Vec<Vec<Vec<RqMatrix>>>,
    h_gar_poly_basis_form_aggregated: &Vec<Vec<Vec<PolynomialRing>>>,
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

#[time_profiler()]
pub fn prove(
    a_matrix: &RqMatrix,
    b_matrix: &Vec<Vec<RqMatrix>>,
    c_matrix: &Vec<Vec<Vec<RqMatrix>>>,
    d_matrix: &Vec<Vec<Vec<RqMatrix>>>,
) -> (St, Tr) {
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
        &b_matrix,
        &c_matrix,
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
                .map(|s_i_poly| s_i_poly.coefficients.clone())
                .flatten()
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
            let inner_product = inner_product_zq_vector(&pai_element, &s_i);
            sum = sum + inner_product;
        }
        p.push(sum);
    }

    assert_eq!(p.len(), double_lambda.value());

    // sanity check: verify p_j = ct(sum(<σ−1(pi_i^(j)), s_i>)) for all i = 1..r
    for j in 0..double_lambda.value() {
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
        assert_eq!(sum.coefficients[0], p[j]);
    }
    println!("Prover: Send proof p");

    // todo: send p to verifier(put in transcript)

    // 3.3 Verifier have to check: || p || <= \sqrt{128} * beta

    // ================================================

    println!("Prover: Do aggregation");
    // 4. GOAL: Aggregation
    // 4.1 psi^(k) is randomly chosen from Z_q^{L}
    // k = 1..λ/log2^q
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
        let inner_product_omega_k_p = inner_product_zq_vector(&omega_k, &p);
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
                    let inner_product_ij = inner_product_polynomial_ring_vector(&phi_i, &s_j)
                        + inner_product_polynomial_ring_vector(&phi_j, &s_i);
                    // Notice we do not divide by 2 here as paper described, because there is no division in the ring, we multiply by 2 instead with other terms to make verifier check work
                    inner_product_ij
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
        &d_matrix,
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
        .map(|(phi_i, c_i)| inner_product_polynomial_ring_vector(&phi_i, &z) * c_i)
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
    let st = St {
        a_constraint,
        phi_constraint,
        b_constraint,
        a_constraint_ct,
        phi_constraint_ct,
        b_constraint_ct,
    };
    let tr = Tr {
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
    return (st, tr);
}

fn verify(
    st: St,
    tr: Tr,
    a_matrix: &RqMatrix,
    b_matrix: &Vec<Vec<RqMatrix>>,
    c_matrix: &Vec<Vec<Vec<RqMatrix>>>,
    d_matrix: &Vec<Vec<Vec<RqMatrix>>>,
) {
    // same parameters as in the prover
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
    let new_beta = Zq::new(250); // Example value for beta

    let St {
        a_constraint,
        phi_constraint,
        b_constraint,
        a_constraint_ct,
        phi_constraint_ct,
        b_constraint_ct,
    } = st;

    let Tr {
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
        .map(|(phi_i, c_i)| inner_product_polynomial_ring_vector(&phi_i, &z) * c_i)
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
        &b_matrix,
        &c_matrix,
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
        &d_matrix,
        &h_gar_poly_basis_form_aggregated,
        t2,
        kappa2,
        size_r,
        size_n,
        deg_bound_d,
    );
    assert_eq!(u2, u2_check);
}

#[cfg(test)]
mod tests {
    use super::*;
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
            .map(|i| {
                (0..size_n.value())
                    .map(|j| generate_random_polynomial_ring(deg_bound_d.value()))
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
                        let inner_product_ij = inner_product_polynomial_ring_vector(&phi_i, &s_j)
                            + inner_product_polynomial_ring_vector(&phi_j, &s_i);
                        inner_product_ij
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
            .map(|(phi_i, c_i)| inner_product_polynomial_ring_vector(&phi_i, &z) * c_i)
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
                        inner_product_polynomial_ring_vector(&s_i, &s_j)
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
                        let inner_product_ij = inner_product_polynomial_ring_vector(&phi_i, &s_j)
                            + inner_product_polynomial_ring_vector(&phi_j, &s_i);
                        inner_product_ij
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
            .map(|i| {
                (0..size_r.value())
                    .map(|j| generate_random_polynomial_ring(deg_bound_d.value()))
                    .collect::<Vec<PolynomialRing>>()
            })
            .collect();

        // generate size_r * size_n phi_aggr
        let phi_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
            .map(|i| {
                (0..size_n.value())
                    .map(|j| {
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
                        inner_product_polynomial_ring_vector(&s_i, &s_j)
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
                        let inner_product_ij = inner_product_polynomial_ring_vector(&phi_i, &s_j)
                            + inner_product_polynomial_ring_vector(&phi_j, &s_i);
                        inner_product_ij
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
                let inner_product_si_sj =
                    inner_product_polynomial_ring_vector(&elem_s_i, &elem_s_j);
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
                        inner_product_polynomial_ring_vector(&s_i, &s_j)
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

    #[test]
    fn test_poly_vec_times_poly() {
        // Arrange
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = PolynomialRing {
            coefficients: vec![Zq::from(2), Zq::from(3), Zq::from(4)],
        };

        // Act
        let result = poly_vec_times_poly(&a, &b);
        // expected[0] = a[0] * b = (1 + 2x + 3x^2) * (2 + 3x + 4x^2) = 2 + 7x + 16x^2 + 17x^3 + 12x^4
        // expected[1] = a[1] * b = (4 + 5x + 6x^2) * (2 + 3x + 4x^2) = 8 + 22x + 43x^2 + 38x^3 + 24x^4
        // Assert
        let expected = vec![
            PolynomialRing {
                coefficients: vec![
                    Zq::from(2),
                    Zq::from(7),
                    Zq::from(16),
                    Zq::from(17),
                    Zq::from(12),
                ],
            },
            PolynomialRing {
                coefficients: vec![
                    Zq::from(8),
                    Zq::from(22),
                    Zq::from(43),
                    Zq::from(38),
                    Zq::from(24),
                ],
            },
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_poly_vec_add_poly_vec() {
        // Arrange
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];

        // Act
        let result = poly_vec_add_poly_vec(&a, &b);

        // Assert
        let expected = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(8), Zq::from(10), Zq::from(12)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(14), Zq::from(16), Zq::from(18)],
            },
        ];
        assert_eq!(result, expected);
    }

    #[test]
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
        // assert_eq!(b_k, 1983);
    }

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

    #[test]
    fn test_num_to_basis() {
        let num = Zq::from(42);
        let basis = Zq::from(2);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(3);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(1),
                Zq::from(0),
                Zq::from(2),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(6);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(4),
                Zq::from(4),
                Zq::from(2),
                Zq::from(0),
                Zq::from(0),
                Zq::from(0)
            ]
        );

        let num = Zq::from(100);
        let basis = Zq::from(10);
        let digits = Zq::from(6);
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(
            binary,
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(0)
            ]
        );
    }

    #[test]
    fn test_basis_to_num_vector() {
        let basis = Zq::from(10);
        let digits = Zq::from(3);
        let vec1 = [8, 46, 61, 71, 33, 33, 18];
        let vec2 = [20, 54, 94, 93, 70, 33, 14];
        let vec3 = [24, 40, 100, 85, 121, 57, 56];
        let vec4 = [14, 37, 91, 118, 159, 109, 72];
        for vec in [vec1, vec2, vec3, vec4] {
            let mut temp_vec: Vec<Vec<Zq>> = Vec::new();
            for i in vec {
                let num = num_to_basis(Zq::from(i), basis, digits);
                temp_vec.push(num);
            }
        }
    }

    #[test]
    fn test_ring_polynomial_to_basis() {
        let poly = PolynomialRing {
            coefficients: vec![Zq::from(42), Zq::from(100), Zq::from(100)],
        };
        let basis = Zq::from(2);
        let digits = Zq::from(8);
        let expected_result = vec![
            vec![
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
            ],
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(1),
                Zq::from(0),
            ],
            vec![
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(0),
                Zq::from(0),
                Zq::from(1),
                Zq::from(1),
                Zq::from(0),
            ],
        ];
        let result = ring_polynomial_to_basis(&poly, basis, digits);
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_inner_product_polynomial_ring() {
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];

        let result = inner_product_polynomial_ring_vector(&a, &b);

        // Expected result calculation:
        // (1 + 2x + 3x^2) * (7 + 8x + 9x^2) = 7 + 22x + 46x^2 + 42x^3 + 27x^4
        // (4 + 5x + 6x^2) * (10 + 11x + 12x^2) = 40 +  96x + 163x^2 + 126x^3 + 72x^4
        // Sum: 47 + 116x + 209x^2 + 168x^3 + 99x^4

        let expected = PolynomialRing {
            coefficients: vec![
                Zq::from(47),
                Zq::from(116),
                Zq::from(209),
                Zq::from(168),
                Zq::from(99),
            ],
        };

        assert_eq!(result.coefficients, expected.coefficients);
    }

    #[test]
    fn test_inner_product_zq_vector() {
        // Arrange
        let a = vec![Zq::from(1), Zq::from(2), Zq::from(3)];
        let b = vec![Zq::from(4), Zq::from(5), Zq::from(6)];

        // Act
        let result = inner_product_zq_vector(&a, &b);

        // Assert
        let expected =
            (Zq::from(1) * Zq::from(4)) + (Zq::from(2) * Zq::from(5)) + (Zq::from(3) * Zq::from(6));
        assert_eq!(
            result, expected,
            "inner_product_zq_vector did not return the correct result"
        );
    }

    #[test]
    fn test_inner_product_poly_matrix_and_poly_vector() {
        // Arrange
        let poly_matrix = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(2)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(4)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(5), Zq::from(6)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(7), Zq::from(8)],
                },
            ],
        ];
        let poly_vector = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(9), Zq::from(10)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(11), Zq::from(12)],
            },
        ];

        // Expected Calculation:
        // u = sum D_ij * h_ij^(k) for all k = 1..(t1-1)
        // For this test case:
        // Result[0] = (1+2x) * (9+10x) + (5+6x) * (11+12x) = 64 + 154x + 92x^2
        // Result[1] = (3+4x) * (9+10x) + (7+8x) * (11+12x) = 104 + 238x + 136x^2
        let expected = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(64), Zq::from(154), Zq::from(92)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(104), Zq::from(238), Zq::from(136)],
            },
        ];

        // Act
        let result = inner_product_poly_matrix_and_poly_vector(&poly_matrix, &poly_vector);

        // Assert
        assert_eq!(result, expected, "The inner product of the polynomial matrix and vector did not produce the expected result.");
    }

    #[test]
    fn test_generate_gaussian_distribution() {
        let nd = Zq::from(10);
        let matrix = generate_gaussian_distribution(nd);
        assert_eq!(matrix.len(), 256);
        assert_eq!(matrix[0].len(), nd.value());
        assert_eq!(matrix[1].len(), nd.value());
        assert!(matrix.iter().all(|row| row
            .iter()
            .all(|&val| val.value == Zq::Q - 1 || val.value == 0 || val.value == 1)));
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

    #[test]
    fn test_conjugation_automorphism() {
        // Create example PolynomialRings a and b
        let a = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let b = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };

        // Compute <a, b>
        let inner_ab = inner_product_zq_vector(&a.coefficients, &b.coefficients);
        assert_eq!(inner_ab.value(), 32);
        // Compute σ_{-1}(a)
        let sigma_inv_a = conjugation_automorphism(&a);
        // Compute <σ_{-1}(a), b>
        let inner_sigma_inv_a_b = &sigma_inv_a * &b;
        // = -12x^{65}-28x^{64}-23x^{63}-12x^{62}+6x^{2}+5x+4
        // = 23x^{63}-12x^{62}+6x^{2}+17x+32 (since x^64 = -1)
        assert_eq!(inner_sigma_inv_a_b.coefficients.len(), 64);
        assert_eq!(inner_sigma_inv_a_b.coefficients[0], Zq::from(32));
        assert_eq!(inner_sigma_inv_a_b.coefficients[1], Zq::from(17));
        assert_eq!(inner_sigma_inv_a_b.coefficients[2], Zq::from(6));
        assert_eq!(inner_sigma_inv_a_b.coefficients[62], Zq::from(Zq::Q - 12));
        assert_eq!(inner_sigma_inv_a_b.coefficients[63], Zq::from(Zq::Q - 23));
        for i in 3..62 {
            assert_eq!(inner_sigma_inv_a_b.coefficients[i], Zq::from(0));
        }

        // Get the constant term of <σ_{-1}(a), b>
        let ct_inner_sigma_inv_a_b = inner_sigma_inv_a_b.coefficients[0];

        // Assert that <a, b> == ct <σ_{-1}(a), b>
        assert_eq!(
            inner_ab, ct_inner_sigma_inv_a_b,
            "<a, b> should equal the constant term of <σ-1(a), b>"
        );
    }

    #[test]
    fn test_decompose_poly_to_basis_form() {
        // Arrange: Create sample input polynomial rings
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(123), Zq::from(456), Zq::from(789)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(12), Zq::from(45), Zq::from(78)],
        };
        let poly_input = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly1.clone(), poly2.clone()],
        ];
        let basis = Zq::from(10);
        let digits = Zq::from(3);

        // Act: Call the function to decompose the polynomial
        let result = poly_matrix_decompose_and_aggregate(&poly_input, basis, digits);

        let expected_row1 = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(6), Zq::from(9)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(0), Zq::from(0), Zq::from(0)],
                },
            ],
        ];

        let expected = vec![expected_row1.clone(), expected_row1.clone()];
        assert_eq!(
            result, expected,
            "The decomposition did not match the expected output."
        );
    }

    #[test]
    fn test_decompose_poly_ring_vector_to_basis_form() {
        // Arrange: Create sample input vector of PolynomialRing
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(123), Zq::from(456), Zq::from(789)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(12), Zq::from(45), Zq::from(78)],
        };
        let poly_vec = vec![poly1.clone(), poly2.clone()];
        let basis = Zq::from(10);
        let digits = Zq::from(3);

        // Act: Call the decomposition function
        let result = poly_vec_decompose_and_aggregate(&poly_vec, basis, digits);

        // Define expected output
        let expected = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(6), Zq::from(9)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(2), Zq::from(5), Zq::from(8)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(4), Zq::from(7)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(0), Zq::from(0), Zq::from(0)],
                },
            ],
        ];

        // Assert
        assert_eq!(
            result, expected,
            "The decomposition did not match the expected output."
        );
    }

    #[test]
    fn test_poly_norm_squared() {
        let poly = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let expected = Zq::from(14); // 1^2 + 2^2 + 3^2 = 14
        let result = poly_norm_squared(&poly);
        assert_eq!(
            result, expected,
            "poly_norm_squared should return the sum of squared coefficients"
        );
    }

    #[test]
    fn test_poly_matrix_norm_squared() {
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let poly_matrix = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly2.clone(), poly1.clone()],
        ];
        // poly_norm_squared(poly1) = 1 + 4 + 9 = 14
        // poly_norm_squared(poly2) = 16 + 25 + 36 = 77
        // Total sum: 14 + 77 + 77 + 14 = 182
        let expected = Zq::from(182);
        let result = poly_matrix_norm_squared(&poly_matrix);
        assert_eq!(result, expected, "poly_matrix_norm_squared should return the sum of squared norms of all polynomials in the matrix");
    }

    #[test]
    fn test_poly_vec_norm_squared() {
        // Arrange
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let vec_polys = vec![poly1.clone(), poly2.clone()];

        // Act
        let result = poly_vec_norm_squared(&vec_polys);

        // Assert
        // poly1 norm: 1^2 + 2^2 + 3^2 = 14
        // poly2 norm: 4^2 + 5^2 + 6^2 = 77
        let expected = Zq::from(14) + Zq::from(77);
        assert_eq!(
            result, expected,
            "poly_vec_norm_squared did not return the correct sum of squared norms"
        );
    }

    #[test]
    fn test_poly_3d_norm_squared() {
        // Arrange
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
        };
        let poly_matrix = vec![
            vec![poly1.clone(), poly2.clone()],
            vec![poly1.clone(), poly2.clone()],
        ];
        let polymat3d = vec![poly_matrix.clone(), poly_matrix.clone()];

        // Act
        let result = poly_3d_norm_squared(&polymat3d);

        // Assert
        // Each poly_matrix contains two vectors of polynomials, each vector has 2 polynomials with norms 14 and 77
        // Each matrix: 2 vectors * (14 + 77) = 2 * 91 = 182
        // Total: 2 matrices * 182 = 364
        let expected = Zq::from(364);
        assert_eq!(
            result, expected,
            "poly_3d_norm_squared did not return the correct sum of squared norms"
        );
    }

    #[test]
    fn test_inner_product_polynomial_ring_vector() {
        // Define sample PolynomialRing vectors
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];
        // (1 + 2x + 3x^2) * (7 + 8x + 9x^2) + (4 + 5x + 6x^2) * (10 + 11x + 12x^2)
        // = (7 + 22x + 46x^2 + 42x^3 + 27x^4) + (40 + 94x + 163x^2 + 126x^3 + 72x^4)
        // = 47 + 116x + 209x^2 + 168x^3 + 99x^4
        let result = inner_product_polynomial_ring_vector(&a, &b);
        let expected = PolynomialRing {
            coefficients: vec![
                Zq::from(47),  // 47
                Zq::from(116), // 116
                Zq::from(209), // 209
                Zq::from(168), // 168
                Zq::from(99),  // 99
            ],
        };

        assert_eq!(result, expected);
    }
}
