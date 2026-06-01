# lau-complex-analysis

> Complex analysis library: holomorphic functions, contour integration, residue theory, conformal mapping, and agent frequency analysis

## What This Does

Complex analysis library: holomorphic functions, contour integration, residue theory, conformal mapping, and agent frequency analysis. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-complex-analysis
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_complex_analysis::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub struct ArgumentPrinciple;
    pub fn winding_number<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
    pub fn count_zeros<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
    pub fn zeros_minus_poles<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
pub struct RoucheTheorem;
    pub fn check_condition<F, G>(
    pub fn count_zeros<F, G, H>(
pub struct RoucheResult 
pub struct HolomorphicResult 
pub struct HolomorphicCheck 
    pub fn new(tolerance: f64) -> Self 
    pub fn check_at<F>(&self, f: F, x: f64, y: f64) -> HolomorphicResult
    pub fn check_region<F>(&self, f: F, points: &[(f64, f64)]) -> Vec<HolomorphicResult>
    pub fn is_entire<F>(&self, f: F, range: (f64, f64), steps: usize) -> bool
    pub fn derivative<F>(f: &F, z: Complex64) -> Complex64
pub struct MobiusTransformation 
    pub fn new(a: Complex64, b: Complex64, c: Complex64, d: Complex64) -> Self 
    pub fn apply(&self, z: Complex64) -> Complex64 
    pub fn inverse(&self) -> Self 
    pub fn compose(&self, other: &MobiusTransformation) -> Self 
    pub fn identity() -> Self 
    pub fn translation(w: Complex64) -> Self 
    pub fn scaling(lambda: Complex64) -> Self 
    pub fn rotation(theta: f64) -> Self 
    pub fn inversion() -> Self 
    pub fn from_three_points(
    pub fn determinant(&self) -> Complex64 
    pub fn fixed_points(&self) -> (Option<Complex64>, Option<Complex64>) 
    pub fn is_elliptic(&self) -> bool 
    pub fn apply_batch(&self, points: &[Complex64]) -> Vec<Complex64> 
pub enum Contour 
    pub fn evaluate(&self, t: f64) -> (Complex64, Complex64) 
    pub fn unit_circle() -> Self 
    pub fn circle(center: (f64, f64), radius: f64) -> Self 
    pub fn reverse(&self) -> Self 
pub struct ContourIntegrator 
    pub fn new(n_points: usize) -> Self 
    pub fn integrate<F>(&self, f: F, contour: &Contour) -> Complex64
    pub fn integrate_with_n<F>(&self, f: F, contour: &Contour, n: usize) -> Complex64
    pub fn integrate_abs<F>(&self, f: F, contour: &Contour) -> f64
pub enum SingularityType 
pub struct Singularity 
pub struct Residue;
    pub fn compute<F>(f: F, a: Complex64, radius: f64) -> Complex64
    pub fn simple_pole<F>(f: F, a: Complex64) -> Complex64
    pub fn pole_of_order<F>(f: F, a: Complex64, m: u32) -> Complex64
    pub fn classify<F>(f: F, a: Complex64) -> Singularity
pub struct ResidueTheorem;
    pub fn apply<F>(
    pub fn verify<F>(
pub fn point_inside_circle(point: Complex64, center: (f64, f64), radius: f64) -> bool 
pub trait ComplexExt 
pub fn exp(z: Complex64) -> Complex64 
pub fn sin(z: Complex64) -> Complex64 
pub fn cos(z: Complex64) -> Complex64 
pub fn log(z: Complex64) -> Complex64 
pub fn pow(z: Complex64, w: Complex64) -> Complex64 
pub fn c(re: f64, im: f64) -> Complex64 
pub struct FrequencyResult 
pub struct ZTransform;
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**55 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
