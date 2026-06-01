# lau-complex-analysis

A Rust library for **complex analysis** — holomorphic functions, contour integration, Cauchy's theorem, residue calculus, Taylor/Laurent series, conformal mapping, and frequency-domain agent analysis.

## What This Does

This crate provides the core building blocks of single-variable complex analysis, implemented from scratch in pure Rust:

- **Complex arithmetic** — polar form, branches of the logarithm, roots, exponential, trigonometric functions
- **Holomorphicity checking** — numerical verification of the Cauchy-Riemann equations at arbitrary points
- **Contour integration** — line segments, circles, arcs, polygons, and fully generic parametric curves, with adaptive quadrature
- **Cauchy's integral formula** — evaluating analytic functions and all their derivatives via contour integrals
- **Power series** — Taylor and Laurent expansion with automatic radius-of-convergence estimation
- **Residue theory** — residue computation at poles of any order, residue theorem for closed contours
- **Argument principle & Rouché's theorem** — counting zeros and poles inside contours
- **Conformal mapping** — Möbius transformations with circle-preserving geometry
- **Agent frequency analysis** — applying complex Fourier techniques to agent signal decomposition

## Key Idea

Complex analysis is the natural language for two-dimensional potential theory, signal processing, and holomorphic dynamics. This library treats complex-valued functions as first-class citizens: every contour is a parameterized curve `z(t)` with its derivative `z'(t)`, and every integral is computed via numerical quadrature over the parameter domain `[0, 1]`.

The **agent frequency** module applies these tools to multi-agent systems — decomposing agent behavior signals into frequency components, finding dominant periodic modes, and analyzing spectral structure using the complex exponential basis.

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-complex-analysis = "0.1"
```

Requires **Rust 2021 edition**.

### Dependencies

| Crate | Purpose |
|-------|---------|
| `num-complex` | `Complex64` type |
| `nalgebra` | Linear algebra (agent frequency analysis) |
| `serde` | Serialization of contours, series, residues |

## Quick Start

### Holomorphicity check

```rust
use num_complex::Complex64;
use lau_complex_analysis::HolomorphicCheck;

let check = HolomorphicCheck::new(1e-6);

// f(z) = z² is entire (holomorphic everywhere)
let result = check.check_at(|z| z * z, 1.0, 2.0);
assert!(result.is_holomorphic);

// f(z) = conj(z) is NOT holomorphic
let result = check.check_at(|z| z.conj(), 1.0, 2.0);
assert!(!result.is_holomorphic);
```

### Contour integration

```rust
use num_complex::Complex64;
use lau_complex_analysis::{Contour, ContourIntegrator};

// Integrate f(z) = 1/z around the unit circle → 2πi
let circle = Contour::Circle { center: (0.0, 0.0), radius: 1.0 };
let integrator = ContourIntegrator::new(10000); // 10000 quadrature points
let result = integrator.integrate(&circle, |z| 1.0 / z);
// result ≈ 2πi
```

### Residue theorem

```rust
use lau_complex_analysis::{Residue, ResidueTheorem};

// f(z) = 1/(z² + 1) has simple poles at z = ±i
let f = |z: Complex64| 1.0 / (z * z + 1.0);

// Residue at z = i is 1/(2i) = -i/2
let res_at_i = Residue::simple_pole(f, Complex64::new(0.0, 1.0));

// Sum of residues × 2πi gives the contour integral
let theorem = ResidueTheorem::new();
let integral = theorem.evaluate(f, &[Complex64::new(0.0, 1.0)], 1000);
```

### Taylor series

```rust
use lau_complex_analysis::TaylorSeries;

// Expand eᶻ around z₀ = 0
let series = TaylorSeries::new(|z: Complex64| z.exp(), Complex64::new(0.0, 0.0), 10);
let coeffs = series.coefficients(); // [1, 1, 1/2, 1/6, ...]
```

### Möbius transformation

```rust
use lau_complex_analysis::MobiusTransformation;
use num_complex::Complex64;

// Map the upper half-plane to the unit disk
let mob = MobiusTransformation::new(
    Complex64::new(1.0, 0.0),  // a
    Complex64::new(-1.0, 0.0), // b  (z₀ = -1)
    Complex64::new(1.0, 0.0),  // c
    Complex64::new(1.0, 0.0),  // d
);

let z = Complex64::new(0.0, 1.0); // i (upper half-plane)
let w = mob.apply(z);
// |w| < 1 — inside the unit disk
```

## API Reference

### `complex` — Extended Complex Arithmetic

| Item | Description |
|------|-------------|
| `ComplexExt` trait | `to_polar`, `from_polar`, `ln`, `ln_branch`, `exp`, `csin`, `ccos`, `powc`, `nth_root`, `all_nth_roots` |

### `holomorphic` — Cauchy-Riemann Verification

| Item | Description |
|------|-------------|
| `HolomorphicCheck` | Numerically checks `∂u/∂x = ∂v/∂y` and `∂u/∂y = -∂v/∂x` at a point |
| `HolomorphicResult` | Struct with partial derivatives, CR error, and boolean verdict |

### `integration` — Contour Integration

| Item | Description |
|------|-------------|
| `Contour` enum | `LineSegment`, `Circle`, `Arc`, `Polygon`, `Parametric` |
| `ContourIntegrator` | Evaluates `∮ f(z) dz` using trapezoidal quadrature |

### `cauchy` — Cauchy's Integral Formula

| Item | Description |
|------|-------------|
| `CauchyIntegral` | `f(z₀) = (1/2πi) ∮ f(z)/(z-z₀) dz` and derivatives `f⁽ⁿ⁾(z₀)` |

### `series` — Taylor & Laurent Series

| Item | Description |
|------|-------------|
| `TaylorSeries` | Coefficients via Cauchy's formula, radius of convergence |
| `LaurentSeries` | Positive and negative power coefficients, annular convergence |

### `residue` — Residue Calculus

| Item | Description |
|------|-------------|
| `Residue` | Compute residues at simple poles and higher-order poles |
| `ResidueTheorem` | `∮ f(z) dz = 2πi Σ Res(f, zₖ)` over closed contours |

### `argument` — Argument Principle & Rouché

| Item | Description |
|------|-------------|
| `ArgumentPrinciple` | `N - P = (1/2πi) ∮ f'(z)/f(z) dz` — count zeros minus poles |
| `RoucheTheorem` | If `|g| < |f|` on contour, `f` and `f+g` have the same number of zeros |

### `conformal` — Conformal Mapping

| Item | Description |
|------|-------------|
| `MobiusTransformation` | `w = (az+b)/(cz+d)` with inverse, composition, circle preservation |

### `agent_frequency` — Agent Signal Analysis

| Item | Description |
|------|-------------|
| `AgentFrequencyAnalysis` | FFT-based decomposition of agent behavior signals into frequency components, spectral analysis, dominant mode extraction |

## How It Works

### Contour representation

Every contour implements `evaluate(t) → (z, dz/dt)` for `t ∈ [0, 1]`. This uniform parameterization makes the integration machinery generic:

```
∮_C f(z) dz = ∫₀¹ f(z(t)) · z'(t) dt
```

The trapezoidal rule with configurable sample count provides the numerical quadrature. For closed contours (circles, polygons), periodicity of the integrand makes this spectrally accurate.

### Residue computation

For a **simple pole** at `z₀`, the residue is computed by:

```
Res(f, z₀) = lim_{z→z₀} (z - z₀) · f(z)
```

evaluated numerically with a small perturbation. For **higher-order poles**, a Laurent coefficient extraction is used.

### Cauchy-Riemann check

Given `f(z) = u(x,y) + iv(x,y)`, the library computes partial derivatives numerically via central differences:

```
∂u/∂x ≈ [u(x+h, y) - u(x-h, y)] / 2h
```

and checks `∂u/∂x = ∂v/∂y` and `∂u/∂y = -∂v/∂x` within tolerance.

## The Math

**Cauchy's integral formula** states that for `f` holomorphic inside and on a simple closed contour `C`:

```
f(z₀) = (1/2πi) ∮_C f(z)/(z - z₀) dz
```

Differentiating `n` times:

```
f⁽ⁿ⁾(z₀) = n!/(2πi) ∮_C f(z)/(z - z₀)ⁿ⁺¹ dz
```

The **residue theorem** generalizes this: for a meromorphic function `f` with isolated singularities `z₁, ..., zₙ` inside `C`:

```
∮_C f(z) dz = 2πi · Σₖ Res(f, zₖ)
```

The **argument principle** relates the winding number of `f(C)` around the origin to the count of zeros and poles:

```
Z - P = (1/2πi) ∮_C f'(z)/f(z) dz = (1/2π) Δ arg f(z)
```

**Möbius transformations** `w = (az+b)/(cz+d)` with `ad - bc ≠ 0` map circles to circles and preserve angles. They form the automorphism group of the Riemann sphere.

**Laurent series** extend Taylor series to annular domains around isolated singularities:

```
f(z) = Σ_{n=-∞}^{∞} aₙ(z - z₀)ⁿ
```

The coefficient `a₋₁` is the residue at `z₀`.

## License

MIT
