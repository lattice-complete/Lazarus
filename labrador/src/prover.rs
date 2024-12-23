use crate::algebra::{PolynomialRing, RqMatrix, Zq};
use crate::setup::setup;
use profiler_macro::time_profiler;
use rand::Rng;

// inner product of 2 vectors of PolynomialRing
fn inner_product_polynomial_ring_vector(
    a: &Vec<PolynomialRing>,
    b: &Vec<PolynomialRing>,
) -> PolynomialRing {
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

// Function to calculate b^(k)
fn calculate_b_constraint(
    s: &Vec<Vec<PolynomialRing>>,
    a_constraint: &Vec<Vec<PolynomialRing>>,
    phi_constraint: &Vec<Vec<PolynomialRing>>,
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
            let inner_product_si_sj = inner_product_polynomial_ring_vector(&elem_s_i, &elem_s_j);
            let a_constr = &a_constraint[i][j];
            b = b + (inner_product_si_sj * a_constr);
        }
        // calculate inner product of s[i] and phi
        for (x, y) in s[i].iter().zip(phi_constraint[i].iter()) {
            b = b + (x * y);
        }
    }

    b
}

// calculate matrix times vector of PolynomialRing
fn matrix_times_vector_poly(a: &RqMatrix, s_i: &Vec<PolynomialRing>) -> Vec<PolynomialRing> {
    a.values
        .iter()
        .map(|row| {
            row.iter()
                .zip(s_i.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<PolynomialRing>>()
        })
        .collect::<Vec<Vec<PolynomialRing>>>()
        .into_iter()
        .flatten()
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
    let one = Zq::from(1);
    let mut base = basis;

    for _ in 0..digits.value() {
        let digit = remainder.clone() % base.clone();
        result.push(digit);
        remainder = remainder.clone() / base.clone();
    }

    while result.len() < digits.value() as usize {
        // push 0 to the highest position
        result.push(zero.clone());
    }

    result
}

// convert ring polynomial to basis
fn ring_polynomial_to_basis(poly: &PolynomialRing, basis: Zq, digits: Zq) -> Vec<Vec<Zq>> {
    poly.coefficients
        .iter()
        .map(|coeff| num_to_basis(coeff.clone(), basis.clone(), digits.clone()))
        .collect()
}

fn generate_gaussian_distribution(nd: Zq) -> Vec<Vec<Zq>> {
    let nd_usize: usize = nd.value() as usize;
    let modulus: usize = Zq::modulus();
    let mut rng = rand::thread_rng();
    let mut matrix = vec![vec![Zq::from(0); nd_usize]; 256]; // Initialize a 256 x nd matrix

    for i in 0..256 {
        for j in 0..nd_usize {
            let random_value: f32 = rng.gen(); // Generate a random float between 0 and 1
            matrix[i][j] = if random_value < 0.25 {
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
fn conjugation_automorphism(poly: &PolynomialRing) -> PolynomialRing {
    let modulus_minus_one = Zq::from(Zq::modulus() - 1);
    let transformed_coeffs: Vec<Zq> = (0..PolynomialRing::DEGREE_BOUND)
        .map(|i| {
            if i < poly.coefficients.len() {
                if i == 0 {
                    poly.coefficients[i].clone()
                } else {
                    poly.coefficients[i].clone() * modulus_minus_one
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
            .map(|_| Zq::from(rng.gen_range(0..10)))
            .collect(),
    }
}

fn decompose_poly_to_basis_form(
    poly: &Vec<Vec<PolynomialRing>>,
    basis: Zq,
    digits: Zq,
) -> Vec<Vec<Vec<PolynomialRing>>> {
    // Decompose h_ij into basis t_1 parts
    let poly_basis_form: Vec<Vec<Vec<Vec<Zq>>>> = poly
        .iter()
        .map(|poly_i| {
            poly_i
                .iter()
                .map(|poly_i_j| ring_polynomial_to_basis(poly_i_j, basis, digits))
                .collect::<Vec<Vec<Vec<Zq>>>>()
        })
        .collect::<Vec<Vec<Vec<Vec<Zq>>>>>();

    // Pick elements at each position across all inner vectors and aggregate them
    let mut poly_basis_form_aggregated: Vec<Vec<Vec<PolynomialRing>>> = Vec::new();
    for (_i, poly_i_basis_form) in poly_basis_form.iter().enumerate() {
        let mut row_results: Vec<Vec<PolynomialRing>> = Vec::new();
        for (_j, poly_i_j_basis_form) in poly_i_basis_form.iter().enumerate() {
            let mut row_results_j: Vec<PolynomialRing> = Vec::new();
            // Get the number of basis parts and the number of loops needed
            let num_basis_needed = poly_i_j_basis_form.len();
            let num_loop_needed = poly_i_j_basis_form.first().map_or(0, |v| v.len());
            for k in 0..num_loop_needed {
                let mut row_k: Vec<Zq> = Vec::new();
                for basis_needed in 0..num_basis_needed {
                    if let Some(num_to_be_pushed) =
                        poly_i_j_basis_form.get(basis_needed).and_then(|v| v.get(k))
                    {
                        row_k.push(num_to_be_pushed.clone());
                    } else {
                        row_k.push(Zq::from(0));
                    }
                }
                row_results_j.push(PolynomialRing {
                    coefficients: row_k,
                });
            } // finish poly_i_j_basis_form calculation
            row_results.push(row_results_j);
        }
        poly_basis_form_aggregated.push(row_results);
    }
    poly_basis_form_aggregated
}

#[time_profiler()]
pub fn prove() {
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
    let beta = Zq::new(50); // Example value for beta
                            // 0. setup
                            // matrices A, B, C, D are common reference string
    let (a_matrix, b_matrix, c_matrix, d_matrix) =
        setup(size_n, size_r, t1, t2, kappa, kappa1, kappa2);

    let witness_s: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|_| {
            (0..size_n.value())
                .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
                .collect()
        })
        .collect();
    println!("s: {:?}", witness_s);
    // Calculate the sum of squared norms
    let mut sum_squared_norms = Zq::new(0);
    for vector in &witness_s {
        let norm_squared: Zq = vector
            .iter()
            .map(|elem| elem.coefficients[0].pow(2))
            .fold(Zq::new(0), |acc, val| acc + val);
        sum_squared_norms = sum_squared_norms + norm_squared;
    }
    println!("sum_squared_norms: {}", sum_squared_norms.value());
    println!("beta^2: {}", beta.pow(2));
    // Check the condition
    assert!(
        sum_squared_norms <= beta.pow(2),
        "The condition is not satisfied: sum of squared norms exceeds beta^2"
    );

    let mut rng = rand::thread_rng();
    let constraint_num_k = Zq::new(6);
    // In DPCS (dot product constraint system), there are k constraints, each constraint has a, phi, and b
    // Generate random a^(k)_{i,j}: k length vector of matrices, each matrix is r x r, and each element is a Zq
    // TODO: Ensure a_ij == a_ji
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
    println!("b_constraint: {:?}", b_constraint);

    // In DPCS(dot product constraint system) for constant terms(ct), there are k constraints, each constraint has a, phi and b.
    // Generate random a^(l)_{i,j}: l length vector of matrix, matrix length is r x r, each element is a Zq
    // todo: aij == aji, refer to paper page 10
    let constraint_num_l = Zq::new(5); // Define L
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
    println!("a_constraint_ct: {:?}", a_constraint_ct);
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
    println!("phi_constraint: {:?}", phi_constraint_ct);
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
    println!("b_constraint_ct: {:?}", b_constraint_ct);

    // let size_n = 5;
    // A: matrix size: kappa * n, each element is PolynomialRing(R_q)
    // calculate t_i = A * s_i for all i = 1..r
    // size of t_i = (kappa * n)R_q * 1R_q = kappa * n
    let t: Vec<Vec<PolynomialRing>> = witness_s
        .iter()
        .map(|s_i| {
            let t_i = matrix_times_vector_poly(&a_matrix, s_i);
            println!("size of t_i: {:?}", t_i.len());
            t_i
        })
        .collect();
    println!("Calculated all t_i: {:?}", t);
    // print size of all_t_i
    println!("size of all_t_i: {:?}", t.len());
    // check size of all_t_i is kappa
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

    let all_t_i_basis_form_aggregated = decompose_poly_to_basis_form(&t, basis, digits);
    println!(
        "all_t_i_basis_form_aggregated: {:?}",
        all_t_i_basis_form_aggregated
    );

    // 2.2.1 get basis b2 same as 2.1.1
    // Calculate g_ij = <s_i, s_j>
    let num_s = Zq::new(witness_s.len());

    // Calculate garbage polynomial g_ij = <s_i, s_j>
    let g_gar_poly: Vec<Vec<PolynomialRing>> = (0..size_r.value())
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
    println!("g_gar_poly: {:?}", g_gar_poly);
    assert_eq!(g_gar_poly.len(), size_r.value());
    assert_eq!(g_gar_poly[0].len(), size_r.value());
    for i in 0..size_r.value() {
        for j in (i + 1)..size_r.value() {
            assert_eq!(
                g_gar_poly[i][j], g_gar_poly[j][i],
                "g_ij is not equal to g_ji at indices ({}, {})",
                i, j
            );
        }
    }

    let g_matrix_aggregated = decompose_poly_to_basis_form(&g_gar_poly, basis, t2);
    println!("g_matrix_aggregated: {:?}", g_matrix_aggregated);

    // 2.3 calculate u1
    // 2.3.1 B & C is randomly chosen similar to A
    // 2.3.2 calculate u1 = sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))
    // B_ik: Rq^{kappa1 x kappa}, t_i: Rq^{kappa}, t_i^(k): Rq^{kappa}
    // B_ik * t_i^(k): Rq^{kappa1}
    // First summation: ∑ B_ik * t_i^(k), 1 ≤ i ≤ r, 0 ≤ k ≤ t1−1
    // Initialize u1 with zeros with size kappa1, each element is a polynomial ring
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
    println!("u1: {:?}", u1);

    // Second summation: ∑ C_ijk * g_ij^(k)
    // Calculate u1 using the pre-generated c_matrix
    for i in 0..size_r.value() {
        for j in i..size_r.value() {
            for k in 0..t2.value() {
                println!("i: {}, j: {}, k: {}", i, j, k);
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
    println!("u1: {:?}", u1);

    // ================================================

    // 3. GOAL: JL projection
    let nd = size_n * deg_bound_d;
    // generate gaussian distribution matrices
    // there are size_r matrices, each matrix size is 256 * nd
    // TODO: should from verifier
    let gaussian_distribution_matrices = (0..size_r.value())
        .map(|_| generate_gaussian_distribution(nd))
        .collect::<Vec<Vec<Vec<Zq>>>>();
    println!(
        "gaussian_distribution_matrices: {:?}",
        gaussian_distribution_matrices
    );
    assert_eq!(gaussian_distribution_matrices.len(), size_r.value());
    assert_eq!(
        gaussian_distribution_matrices[0].len(),
        double_lambda.value()
    );
    assert_eq!(gaussian_distribution_matrices[0][0].len(), nd.value());
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
    println!("s_coeffs: {:?}", s_coeffs);
    // implement p calculation, inner product of gaussian_distribution_matrices and s_coeffs
    let mut p: Vec<Zq> = Vec::with_capacity(double_lambda.value());
    for j in 0..double_lambda.value() {
        let mut sum = Zq::new(0);
        for i in 0..size_r.value() {
            let pai = &gaussian_distribution_matrices[i][j];
            let s_i = &s_coeffs[i];
            let inner_product = inner_product_zq_vector(&pai, &s_i);
            sum = sum + inner_product;
        }
        p.push(sum);
    }
    println!("p: {:?}", p);
    assert_eq!(p.len(), double_lambda.value());

    // sanity check: verify p_j = ct(sum(<σ−1(pi_i^(j)), s_i>)) for all i = 1..r
    for j in 0..double_lambda.value() {
        let mut sum = PolynomialRing {
            coefficients: vec![Zq::from(0); deg_bound_d.value()],
        };
        for i in 0..size_r.value() {
            let pai = &gaussian_distribution_matrices[i][j];
            let s_i = &s_coeffs[i];
            let pai_poly = PolynomialRing {
                coefficients: pai.clone(),
            };
            let pai_poly_ca = conjugation_automorphism(&pai_poly);
            let s_i_poly = PolynomialRing {
                coefficients: s_i.clone(),
            };
            sum = sum + &pai_poly_ca * s_i_poly;
        }
        println!("sum: {:?}", sum);
        assert_eq!(sum.coefficients[0], p[j]);
    }

    // todo: send p to verifier(put in transcript)

    // 3.3 Verifier have to check: || p || <= \sqrt{128} * beta

    // ================================================

    // 4. GOAL: Aggregation
    // 4.1 psi^(k) is randomly chosen from Z_q^{L}
    // k = 1..λ/log2^q
    let size_k = lambda / log_q;
    let psi_challenge: Vec<Vec<Zq>> = (0..size_k.value())
        .map(|_| {
            (0..constraint_num_l.value())
                .map(|_| Zq::new(rng.gen_range(0..10)))
                .collect()
        })
        .collect();
    assert_eq!(psi_challenge.len(), size_k.value());
    assert_eq!(psi_challenge[0].len(), constraint_num_l.value());

    // 4.2 omega^(k) is randomly chosen from Z_q^{256}
    //      (Both using Guassian Distribution)
    let omega_challenge: Vec<Vec<Zq>> = (0..size_k.value())
        .map(|_| {
            (0..double_lambda.value())
                .map(|_| Zq::new(rng.gen_range(0..10)))
                .collect()
        })
        .collect();
    assert_eq!(omega_challenge.len(), size_k.value());
    assert_eq!(omega_challenge[0].len(), double_lambda.value());

    // 4.3 caculate b^{''(k)}
    // 4.3.1 calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
    let a_constraint_ct_aggr: Vec<Vec<Vec<PolynomialRing>>> = (0..size_k.value())
        .map(|k| {
            let psi_k = &psi_challenge[k];
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
    println!("a_constraint_ct_aggr: {:?}", a_constraint_ct_aggr);
    assert_eq!(a_constraint_ct_aggr.len(), size_k.value());
    assert_eq!(a_constraint_ct_aggr[0].len(), size_r.value());
    assert_eq!(a_constraint_ct_aggr[0][0].len(), size_r.value());
    // 4.3.2 calculate phi_i^{''(k)} =
    //       sum(psi_l^(k) * phi_i^{'(l)}) for all l = 1..L
    //       + sum(omega_j^(k) * sigma_{-1} * pi_i^{j)) for all j = 1..256
    let phi_ct_aggr: Vec<Vec<Vec<PolynomialRing>>> = (0..size_k.value())
        .map(|k| {
            (0..size_r.value())
                .map(|i| {
                    // Part 1: sum(psi_l^(k) * phi_constraint_ct[l][i] for all l)
                    let part1: Vec<PolynomialRing> = (0..constraint_num_l.value())
                        .map(|l| {
                            let psi = psi_challenge[k][l];
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
                            let omega = omega_challenge[k][j];
                            gaussian_distribution_matrices[i][j]
                                .chunks(deg_bound_d.value())
                                .take(size_n.value())
                                .map(|chunk| {
                                    let pai_poly = PolynomialRing {
                                        coefficients: chunk.to_vec(),
                                    };
                                    conjugation_automorphism(&pai_poly) * omega
                                })
                                .collect::<Vec<PolynomialRing>>()
                        })
                        .fold(
                            vec![
                                PolynomialRing {
                                    coefficients: vec![Zq::from(0); deg_bound_d.value()]
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
    println!("phi_ct_aggr: {:?}", phi_ct_aggr);
    assert_eq!(phi_ct_aggr.len(), size_k.value());
    assert_eq!(phi_ct_aggr[0].len(), size_r.value());

    // 4.3.3 calculate b^{''(k)} = sum(a_ij^{''(k)} * <s_i, s_j>) + sum(<phi_i^{''(k)}, s_i>)
    let b_aggr: Vec<PolynomialRing> = (0..size_k.value())
        .map(|k| {
            (0..size_r.value())
                .map(|i| {
                    (0..size_r.value())
                        .map(|j| {
                            &a_constraint_ct_aggr[k][i][j]
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
        .collect();
    println!("b_aggr: {:?}", b_aggr);
    assert_eq!(b_aggr.len(), size_k.value());

    // todo: send b^{''(k)} to verifier

    // Verifier check: b_0^{''(k)} ?= <⟨omega^(k),p⟩> + sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L
    for k in 0..size_k.value() {
        let b_k_0_from_poly: Zq = b_aggr[k].coefficients[0];
        // sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L
        let mut b_k_0_computed: Zq = (0..constraint_num_l.value())
            .map(|l| {
                let psi_k_l = psi_challenge[k][l];
                let b_l_0 = b_constraint_ct[l];
                psi_k_l * b_l_0
            })
            .sum();
        // <⟨omega^(k),p⟩>
        let omega_k = &omega_challenge[k];
        let inner_product_omega_k_p = inner_product_zq_vector(&omega_k, &p);
        // add them together
        b_k_0_computed += inner_product_omega_k_p;
        // print k
        println!("k: {}", k);
        assert_eq!(b_k_0_from_poly, b_k_0_computed);
    }

    // ================================================

    // 5. GOAL: Calculate u2 (2nd outer commitment)
    // 5.1 vec<alpha> and vec<beta> are randomly chosen from R_q^{K} and R_q^{128/logQ}
    let alpha_challenge: Vec<PolynomialRing> = (0..constraint_num_k.value())
        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
        .collect();
    let beta_challenge: Vec<PolynomialRing> = (0..size_k.value())
        .map(|_| generate_random_polynomial_ring(deg_bound_d.value()))
        .collect();
    // 5.2 phi_i = sum(alpha_k * phi_i) + beta_k * phi_i^{''(k)}
    let phi_aggr: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|i| {
            // Part 1: sum(alpha_k * phi_i)
            let part1: Vec<PolynomialRing> = (0..constraint_num_k.value())
                .map(|k| {
                    let alpha = &alpha_challenge[k];
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
                    let beta = &beta_challenge[k];
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
    println!("phi_aggr: {:?}", phi_aggr);
    assert_eq!(phi_aggr.len(), size_r.value());
    assert_eq!(phi_aggr[0].len(), size_n.value());
    // 5.3 Calculate garbage polynomial h_ij = 1/2 * (<phi_i, s_j> + <phi_j, s_i>)
    let h_gar_poly: Vec<Vec<PolynomialRing>> = (0..size_r.value())
        .map(|i| {
            (0..size_r.value())
                .map(|j| {
                    let phi_i = &phi_aggr[i];
                    let phi_j = &phi_aggr[j];
                    let s_i = &witness_s[i];
                    let s_j = &witness_s[j];
                    let inner_product_ij = inner_product_polynomial_ring_vector(&phi_i, &s_j)
                        + inner_product_polynomial_ring_vector(&phi_j, &s_i);
                    inner_product_ij / Zq::from(2)
                })
                .collect::<Vec<PolynomialRing>>()
        })
        .collect();
    println!("h_gar_poly: {:?}", h_gar_poly);
    assert_eq!(h_gar_poly.len(), size_r.value());
    assert_eq!(h_gar_poly[0].len(), size_r.value());
    for i in 0..size_r.value() {
        for j in (i + 1)..size_r.value() {
            assert_eq!(
                h_gar_poly[i][j], h_gar_poly[j][i],
                "h_ij is not equal to h_ji at indices ({}, {})",
                i, j
            );
        }
    }

    let h_gar_poly_basis_form_aggregated = decompose_poly_to_basis_form(&h_gar_poly, basis, digits);
    println!(
        "h_gar_poly_basis_form_aggregated: {:?}",
        h_gar_poly_basis_form_aggregated
    );

    // 5.4 u2 = sum D_ij * h_ij^(k) for all k = 1..(t1-1)
    let u2 = (0..size_r.value())
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
                println!("i: {}, j: {}, k: {}", i, j, k);
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
        );

    println!("u2: {:?}", u2);

    // Send u2 to verifier
    // transcript.add(u2)

    // ================================================

    // 6. GOAL: calculate z (Amortized Opening)
    // 6.1 c_i is randomly chosen from C, i = 1..r
    // todo: get c from challenge space, refer to paper page 6
    let c_challenge: Vec<Vec<Zq>> = (0..size_r.value())
        .map(|_| {
            (0..size_n.value())
                .map(|_| Zq::new(rng.gen_range(0..10)))
                .collect()
        })
        .collect();
    // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
    let z: Vec<PolynomialRing> = (0..size_r.value())
        .map(|i| {
            let c_i = &c_challenge[i];
            let s_i = &witness_s[i];
            c_i.iter().zip(s_i.iter()).map(|(c, s)| s * *c).fold(
                PolynomialRing {
                    coefficients: vec![Zq::from(0); size_n.value()],
                },
                |acc, x| acc + x,
            )
        })
        .collect();
    println!("z: {:?}", z);

    // Send z, t_i, g_ij, h_ij to verifier
    // transcript.add(z);
    // return transcript;
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_setup_prover() {
        prove();
    }

    #[test]
    fn test_multiply_by_polynomial_ring() {
        let poly1 = PolynomialRing {
            coefficients: vec![Zq::from(1), Zq::from(2)],
        };
        let poly2 = PolynomialRing {
            coefficients: vec![Zq::from(3), Zq::from(4)],
        };
        let result = poly1 * poly2;
        assert_eq!(
            result.coefficients,
            vec![Zq::from(3), Zq::from(10), Zq::from(8)]
        ); // 1*3, 1*4 + 2*3, 2*4
    }
    #[test]
    fn test_polynomial_ring_mul_overflow() {
        // Create two polynomials that will cause overflow when multiplied
        // For example, (X^63 + 1) * (X^63 + 1) = X^126 + 2X^63 + 1
        // Modulo X^64 + 1, X^64 = -1, so X^126 = X^(2*64 -2) = X^-2 = X^62
        // Thus, X^126 + 2X^63 +1 mod X^64+1 = (-1)*X^62 + 2X^63 +1

        // Initialize poly1 as X^63 + 1
        let mut poly1_coeffs = vec![Zq::from(0); 64];
        poly1_coeffs[0] = Zq::from(1); // Constant term
        poly1_coeffs[63] = Zq::from(1); // X^63 term
        let poly1 = PolynomialRing {
            coefficients: poly1_coeffs,
        };

        // Multiply poly1 by itself
        let product = poly1.clone() * poly1.clone();

        // Expected coefficients after reduction modulo X^64 + 1:
        // coefficients[0] = 1
        // coefficients[62] = Zq::modulus() - 1  (since -1 mod q)
        // coefficients[63] = 2
        // All other coefficients should be 0
        let mut expected_coeffs = vec![Zq::from(1)];
        for _ in 1..62 {
            expected_coeffs.push(Zq::from(0));
        }
        expected_coeffs.push(Zq::from(Zq::modulus() - 1)); // X^62 term
        expected_coeffs.push(Zq::from(2)); // X^63 term

        // Assert that the product has the correct degree bound
        assert_eq!(
            product.coefficients.len(),
            64,
            "Product should be truncated to DEGREE_BOUND"
        );

        // Assert that the coefficients match the expected values
        assert_eq!(
            product.coefficients, expected_coeffs,
            "Overflow handling in multiplication is incorrect"
        );
    }

    #[test]
    fn test_calculate_b_k() {
        let r = 3;
        let n = 4;

        let s: Vec<Vec<PolynomialRing>> = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(13), Zq::from(14), Zq::from(15)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(16), Zq::from(17), Zq::from(18)],
                },
            ],
        ];

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
        println!("b_k: {:?}", b_k);
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
        let result = matrix_times_vector_poly(&a, &s_i);
        assert_eq!(result.len(), a.values.len() * s_i.len()); // Check that the result length is correct
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
                println!("num: {:?}", num);
                temp_vec.push(num);
            }
            println!("temp_vec for vec {:?}: {:?}", vec, temp_vec);
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
    fn test_generate_gaussian_distribution() {
        let nd = Zq::from(10);
        let matrix = generate_gaussian_distribution(nd);
        println!("matrix: {:?}", matrix);
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
    fn test_zq_addition() {
        let a = Zq::new(10);
        let b = Zq::new(20);
        let result = a + b;
        assert_eq!(result.value, 30);
    }

    #[test]
    fn test_zq_subtraction() {
        let a = Zq::new(10);
        let b = Zq::new(5);
        let result = a - b;
        assert_eq!(result.value, 5);
    }

    #[test]
    fn test_zq_multiplication() {
        let a = Zq::new(6);
        let b = Zq::new(7);
        let result = a * b;
        assert_eq!(result.value, 42);
    }

    #[test]
    fn test_zq_overflow() {
        let a = Zq::new(Zq::Q - 1);
        let b = Zq::new(2);
        let result = a + b;
        assert_eq!(result.value, 1); // (2^32 - 1) + 2 mod 2^32 = 1
    }

    #[test]
    fn test_zq_new() {
        let value = 4294967297; // Q + 1
        let zq = Zq::new(value);
        assert_eq!(zq.value, 1);
    }

    #[test]
    fn test_zq_division() {
        let a = Zq::new(20);
        let b = Zq::new(5);
        let result = a / b;
        assert_eq!(result.value, 4);
    }

    #[test]
    fn test_zq_remainder() {
        let a = Zq::new(10);
        let b = Zq::new(3);
        let result = a % b;
        assert_eq!(result.value, 1);
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
        println!("inner_ab: {:?}", inner_ab);
        assert_eq!(inner_ab.value(), 32);
        // Compute σ_{-1}(a)
        let sigma_inv_a = conjugation_automorphism(&a);
        println!("sigma_inv_a: {:?}", sigma_inv_a);
        // Compute <σ_{-1}(a), b>
        let inner_sigma_inv_a_b = &sigma_inv_a * &b;
        println!("inner_sigma_inv_a_b: {:?}", inner_sigma_inv_a_b);

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
        let poly_input = vec![vec![poly1.clone(), poly2.clone()]];
        let basis = Zq::from(10);
        let digits = Zq::from(3);

        // Act: Call the function to decompose the polynomial
        let result = decompose_poly_to_basis_form(&poly_input, basis, digits);

        let expected = vec![vec![
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
        ]];
        assert_eq!(
            result, expected,
            "The decomposition did not match the expected output."
        );
    }
}
