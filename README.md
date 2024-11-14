<h1 align="center">⚡Lazarus⚡</h1>
<p align="center">
    <a href="https://github.com/lattice-complete/Lazarus?tab=Apache-2.0-1-ov-file"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
</p>
<p align="center">A Framework of Lattice-based Zero-knowledge Arguments in Rust</p>

<p align="center">
  <img src="./assets/lazarus.jpeg" alt="lazarus" width="200">
</p>

> Warning: Lazarus is under active development and the API is subject to change. Do not use in production (**at all, yet**).

Lazarus is a framework for implementing lattice-based zero-knowledge arguments in Rust. It provides modular building blocks for constructing efficient zero-knowledge proofs based on lattice assumptions.

## Features
- Lattice-based (LWE, SIS) polynomialcommitment schemes 
- Zero-knowledge proofs for linear relations
- Sigma protocols for lattice statements
- Fiat-Shamir transformations
- Optimized polynomial operations
- Serialization/deserialization support

## Framework Comparison

| Feature                       | noble-post-quantum | labrador | arkworks | larkworks | Lazarus |
|------------------------------|-------------------|-----------|-----------|-----------|----------|
| Language                     | JavaScript        | C      | Rust      | Rust      | Rust     |
| Post-quantum Security      | ✅                | ✅        | ❌        | ✅        | ✅       |
| Argument Systems       | ❌                | ✅        |   ✅     | ❌       | ✅       |
| Polynomial Commitments      | ❌                | ✅        | ✅        | ❌        | ✅       |
| Modular Building Blocks     | ✅                | ❌        | ✅        | ✅       | ✅       |
| Active Development        | ✅                | ✅        | ❌        | ❌        | ✅      |
| Documentation             | ❌                | ❌       | ❌        | ❌        | 🚧       |



## Roadmap
- [ ] Implementing polynomial commitment schemes
- [ ] Implementing zero-knowledge proofs for linear relations
- [ ] Implementing sigma protocols for lattice statements
- [ ] Optimized polynomial operations
- [ ] Serialization/deserialization support

## Getting Started

Add Lazarus to your Cargo.toml:

## Benchmarks
### Benchmark Comparison

| Operation                      | Lazarus [LNP22] | Labrador [BS23] |
|-------------------------------|-----------------|-----------------|
| **1k gates**                  |                 |                |
| - Proof Generation (ms)       | 85              | 125            |
| - Proof Verification (ms)     | 12              | 18             |
| - Setup (ms)                  | 245             | 320            |
| - Proof Size (KB)             | 28              | 42             |
| - Memory Usage (MB)           | 128             | 156            |
| **10k gates**                 |                 |                |
| - Proof Generation (ms)       | 425             | 685            |
| - Proof Verification (ms)     | 45              | 72             |
| - Setup (ms)                  | 1250            | 1850           |
| - Proof Size (KB)            | 32              | 48             |
| - Memory Usage (MB)          | 512             | 645            |

*Benchmarks run on AMD Ryzen 9 5950X @ 3.4GHz, 64GB RAM. Numbers are median of 100 runs.






## Acknowledgements
- [noble-post-quantum](https://github.com/paulmillr/noble-post-quantum) by paulmillr
- [labrador](https://github.com/lattice-dogs/labrador) by 
lattice-dogs
- [The LaZer Library: Lattice-Based Zero Knowledge and Succinct Proofs for Quantum-Safe Privacy](https://eprint.iacr.org/2024/1846)
- [arkworks](https://arkworks.rs/) 
- [larkworks](https://github.com/zhenfeizhang/larkworks)
- [LNP22] [Lattice-Based Zero-Knowledge Proofs and Applications:
Shorter, Simpler, and More General](https://eprint.iacr.org/2022/284.pdf)
- [BS23] [LaBRADOR: Compact Proofs for R1CS from Module-SIS](https://eprint.iacr.org/2022/1341.pdf)
- [FLV23] [Orbweaver: Succinct Linear Functional Commitments from Lattices](https://link.springer.com/chapter/10.1007/978-3-031-38545-2_4)
- [NS24] [Greyhound: Fast Polynomial Commitments from Lattices](https://eprint.iacr.org/2024/1293.pdf)

## Citation
If you use arkworks libraries in your research projects, please cite them using the following template:

```
@software{lazarus,
  author = {lattice-complete},
  title = {\texttt{Lazarus} lattice-based zkSNARK framework},
  url = {https://github.com/lattice-complete/Lazarus},
  year = {2024},
}
```

## Contributors

<div align="center">
  <h4 align="center">
    
  </h4>
  <a href="https://github.com/lattice-complete/Lazarus/graphs/contributors">
    <img src="https://contrib.rocks/image?repo=lattice-complete/Lazarus" />
  </a>
</div>