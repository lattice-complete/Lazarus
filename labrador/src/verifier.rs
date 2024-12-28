// src/verifier.rs

use profiler_macro::time_profiler;

// What does the verifier already know? (Statement)
// 1. public paraments: [a_ij^(k), a_ij^(l), phi^(k), phi^(l), b^(k), b0(l)']
// 2. generate from V: [PI_i, psi^(k), omega^(k), vec<alpha>, vec<beta>, c_i]
// What are sent to the verifier?
// [u1, p, b^{''(k)},u2, z, t_i, g_ij, h_ij]
#[time_profiler()]
pub fn verify() {
    println!("Verifying something...");
    // 1. g_ij ?= g_ji
    // 2. h_ij ?= h_ji
    // 3. Check norm_square < beta_square:
    // 3.1 ||z^(i)||^2 + sum(t_i^(k)) * sum(||g_ij||^2) + sum(||h_ij||^2) < beta_square
    // ...
}
