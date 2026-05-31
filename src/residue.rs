//! Residue theorem: pole detection, residue computation, contour integrals.

use num_complex::Complex64;
use crate::integration::{Contour, ContourIntegrator};
use crate::series::LaurentSeries;
use serde::{Deserialize, Serialize};

/// Type of singularity.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SingularityType {
    Removable,
    SimplePole,
    PoleOfOrder(u32),
    EssentialSingularity,
    BranchPoint,
}

/// Information about a singularity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Singularity {
    pub point: Complex64,
    pub singularity_type: SingularityType,
    pub residue: Complex64,
    pub order: Option<u32>,
}

/// Residue computation.
pub struct Residue;

impl Residue {
    /// Compute the residue of f at point a using the Laurent series approach.
    /// Res(f, a) = (1/2πi) ∮_C f(z) dz around a small circle.
    pub fn compute<F>(f: F, a: Complex64, radius: f64) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let integrator = ContourIntegrator::new(10000);
        let contour = Contour::circle((a.re, a.im), radius);
        let result = integrator.integrate(f, &contour);
        result / (2.0 * std::f64::consts::PI * Complex64::i())
    }

    /// Compute residue of a simple pole using:
    /// Res(f, a) = lim_{z→a} (z-a) f(z)
    pub fn simple_pole<F>(f: F, a: Complex64) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let h = 1e-10;
        (f(a + h) * h + f(a - h) * (-h)) / (2.0 * h)
    }

    /// Compute residue of a pole of order m:
    /// Res(f, a) = 1/(m-1)! lim_{z→a} d^(m-1)/dz^(m-1) [(z-a)^m f(z)]
    pub fn pole_of_order<F>(f: F, a: Complex64, m: u32) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        if m == 1 {
            return Self::simple_pole(f, a);
        }
        // Numerical (m-1)th derivative of g(z) = (z-a)^m * f(z)
        let g = |z: Complex64| (z - a).powf(m as f64) * f(z);
        // Compute (m-1)th derivative via Cauchy integral
        let r = 0.01;
        let n = 10000;
        let dt = 2.0 * std::f64::consts::PI / n as f64;
        let mp1 = m as f64; // we need the (m-1)th derivative
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..n {
            let theta = k as f64 * dt;
            let z = a + r * Complex64::from_polar(1.0, theta);
            let dz = r * Complex64::from_polar(1.0, theta) * Complex64::i();
            sum += g(z) / (z - a).powf(mp1) * dz * dt;
        }
        let factorial = (1..m).fold(1.0f64, |acc, k| acc * k as f64);
        sum / (2.0 * std::f64::consts::PI * Complex64::i()) / factorial
    }

    /// Detect the type of singularity at a point.
    pub fn classify<F>(f: F, a: Complex64) -> Singularity
    where
        F: Fn(Complex64) -> Complex64 + Clone,
    {
        // Compute Laurent coefficients near a
        let ls = LaurentSeries::compute(f.clone(), a, -5, 5, 0.001, 0.1);
        let residue = ls.coefficient(0).unwrap_or(Complex64::new(0.0, 0.0));

        // Find the lowest non-zero negative coefficient
        let mut pole_order: Option<u32> = None;
        for (i, c) in ls.coefficients.iter().enumerate() {
            let n = ls.n_min + i as i32;
            if n < 0 && c.norm() > 1e-10 {
                pole_order = Some((-n) as u32);
                break;
            }
        }

        let singularity_type = match pole_order {
            None => SingularityType::Removable,
            Some(1) => SingularityType::SimplePole,
            Some(m) if m <= 5 => SingularityType::PoleOfOrder(m),
            _ => SingularityType::EssentialSingularity,
        };

        Singularity {
            point: a,
            singularity_type,
            residue,
            order: pole_order,
        }
    }
}

/// Apply the residue theorem.
pub struct ResidueTheorem;

impl ResidueTheorem {
    /// Compute ∮_C f(z) dz = 2πi Σ Res(f, a_k)
    /// where the sum is over singularities a_k inside C.
    pub fn apply<F>(
        f: F,
        contour: &Contour,
        singularities: &[Complex64],
        probe_radius: f64,
    ) -> Complex64
    where
        F: Fn(Complex64) -> Complex64 + Clone,
    {
        let two_pi_i = 2.0 * std::f64::consts::PI * Complex64::i();
        let sum: Complex64 = singularities
            .iter()
            .map(|&a| Residue::compute(f.clone(), a, probe_radius))
            .sum();
        two_pi_i * sum
    }

    /// Numerically verify: compute the contour integral directly and compare
    /// with residue sum.
    pub fn verify<F>(
        f: F,
        contour: &Contour,
        singularities: &[Complex64],
        probe_radius: f64,
    ) -> (Complex64, Complex64)
    where
        F: Fn(Complex64) -> Complex64 + Clone,
    {
        let integrator = ContourIntegrator::new(50000);
        let numerical = integrator.integrate(f.clone(), contour);
        let residue_result = Self::apply(f, contour, singularities, probe_radius);
        (numerical, residue_result)
    }
}

/// Helper: check if a point is inside a contour (rough check for circles).
pub fn point_inside_circle(point: Complex64, center: (f64, f64), radius: f64) -> bool {
    let dx = point.re - center.0;
    let dy = point.im - center.1;
    dx * dx + dy * dy < radius * radius
}
