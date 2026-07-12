use algebra::{polynomial_ring::PolynomialRing, zq::Zq};
pub struct Transcript {
    pub u1: Vec<PolynomialRing>,        // Replace with the actual type
    pub pai: Vec<Vec<Vec<Zq>>>,         // Replace with the actual type
    pub p: Vec<Zq>,                     // Replace with the actual type
    pub psi: Vec<Vec<Zq>>,              // Replace with the actual type
    pub omega: Vec<Vec<Zq>>,            // Replace with the actual type
    pub b_ct_aggr: Vec<PolynomialRing>, // Replace with the actual type
    pub alpha: Vec<PolynomialRing>,     // Replace with the actual type
    pub beta: Vec<PolynomialRing>,      // Replace with the actual type
    pub u2: Vec<PolynomialRing>,
    pub c: Vec<PolynomialRing>,
    pub z: Vec<PolynomialRing>,
    pub t: Vec<Vec<PolynomialRing>>, // Replace with the actual type
    pub g: Vec<Vec<PolynomialRing>>, // Replace with the actual type
    pub h: Vec<Vec<PolynomialRing>>,
}
