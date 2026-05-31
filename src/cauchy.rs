//! Cauchy's integral theorem and formula.

use num_complex::Complex64;
use crate::integration::{Contour, ContourIntegrator};

/// Cauchy integral formula computations.
pub struct CauchyIntegral {
    integrator: ContourIntegrator,
}

impl Default for CauchyIntegral {
    fn default() -> Self {
        Self::new(50000)
    }
}

impl CauchyIntegral {
    pub fn new(n_points: usize) -> Self {
        Self {
            integrator: ContourIntegrator::new(n_points),
        }
    }

    /// Cauchy's integral formula: f(a) = (1/2πi) ∮_C f(z)/(z-a) dz
    /// Returns f(a) evaluated via the integral.
    pub fn cauchy_formula<F>(&self, f: F, contour: &Contour, a: Complex64) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let integrand = |z: Complex64| f(z) / (z - a);
        let result = self.integrator.integrate(integrand, contour);
        result / (2.0 * std::f64::consts::PI * Complex64::i())
    }

    /// Generalized Cauchy formula for derivatives:
    /// f^(n)(a) = n!/(2πi) ∮_C f(z)/(z-a)^(n+1) dz
    pub fn cauchy_derivative<F>(&self, f: F, contour: &Contour, a: Complex64, n: u32) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let np1 = (n + 1) as f64;
        let integrand = |z: Complex64| f(z) / (z - a).powf(np1);
        let result = self.integrator.integrate(integrand, contour);
        let factorial = (1..=n).fold(1.0f64, |acc, k| acc * k as f64);
        result * factorial / (2.0 * std::f64::consts::PI * Complex64::i())
    }

    /// Cauchy's integral theorem: if f is holomorphic on and inside C,
    /// then ∮_C f(z) dz = 0. Returns the integral (should be ~0).
    pub fn cauchy_theorem<F>(&self, f: F, contour: &Contour) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        self.integrator.integrate(f, contour)
    }
}
