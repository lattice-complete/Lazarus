use algebra::{polynomial_ring::PolynomialRing, utils::*, zq::Zq};

// Function to calculate b^(k)
pub fn calculate_b_constraint(
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
