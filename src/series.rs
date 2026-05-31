//! Taylor and Laurent series: coefficient computation and convergence.

use num_complex::Complex64;
use crate::cauchy::CauchyIntegral;
use crate::integration::Contour;
use serde::{Deserialize, Serialize};

/// Taylor series expansion around a point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaylorSeries {
    /// Center of expansion.
    pub center: Complex64,
    /// Coefficients a_n such that f(z) = Σ a_n (z-c)^n.
    pub coefficients: Vec<Complex64>,
    /// Estimated radius of convergence.
    pub radius_of_convergence: f64,
}

impl TaylorSeries {
    /// Compute Taylor coefficients via Cauchy's formula:
    /// a_n = f^(n)(c) / n! = (1/2πi) ∮ f(z)/(z-c)^(n+1) dz
    pub fn compute<F>(f: F, center: Complex64, n_terms: usize, contour_radius: f64) -> Self
    where
        F: Fn(Complex64) -> Complex64 + Clone,
    {
        let cauchy = CauchyIntegral::new(30000);
        let contour = Contour::circle((center.re, center.im), contour_radius);

        let coefficients: Vec<Complex64> = (0..n_terms)
            .map(|n| {
                cauchy.cauchy_derivative(f.clone(), &contour, center, n as u32)
            })
            .collect();

        let radius_of_convergence = Self::estimate_radius(&coefficients);

        TaylorSeries {
            center,
            coefficients,
            radius_of_convergence,
        }
    }

    /// Compute Taylor coefficients numerically using finite differences.
    /// More reliable for simple functions.
    pub fn compute_numerical<F>(f: F, center: Complex64, n_terms: usize) -> Self
    where
        F: Fn(Complex64) -> Complex64,
    {
        let h = 1e-5;
        let coefficients: Vec<Complex64> = (0..n_terms)
            .map(|n| {
                // Numerical nth derivative via Cauchy integral on small circle
                let r = 0.1; // small radius for numerical stability
                let m = 2000;
                let dt = 2.0 * std::f64::consts::PI / m as f64;
                let np1 = (n + 1) as f64;
                let mut sum = Complex64::new(0.0, 0.0);
                for k in 0..m {
                    let theta = k as f64 * dt;
                    let z = center + r * Complex64::from_polar(1.0, theta);
                    let dz = r * Complex64::from_polar(1.0, theta) * Complex64::i();
                    sum += f(z) / (z - center).powf(np1) * dz * dt;
                }
                sum / (2.0 * std::f64::consts::PI * Complex64::i())
            })
            .collect();

        let radius_of_convergence = Self::estimate_radius(&coefficients);

        TaylorSeries {
            center,
            coefficients,
            radius_of_convergence,
        }
    }

    /// Evaluate the Taylor series at a point.
    pub fn evaluate(&self, z: Complex64) -> Complex64 {
        let diff = z - self.center;
        let mut result = Complex64::new(0.0, 0.0);
        let mut power = Complex64::new(1.0, 0.0);
        for coeff in &self.coefficients {
            result += coeff * power;
            power *= diff;
        }
        result
    }

    /// Estimate radius of convergence using root test.
    fn estimate_radius(coefficients: &Vec<Complex64>) -> f64 {
        if coefficients.is_empty() {
            return f64::INFINITY;
        }
        let n = coefficients.len();
        let limsup = coefficients
            .iter()
            .skip(1)
            .enumerate()
            .map(|(i, c)| c.norm().powf(1.0 / (i + 1) as f64))
            .fold(0.0f64, f64::max);

        if limsup == 0.0 {
            f64::INFINITY
        } else {
            1.0 / limsup
        }
    }
}

/// Laurent series expansion around a point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LaurentSeries {
    /// Center of expansion.
    pub center: Complex64,
    /// Coefficients from order n_min to n_max.
    /// coefficients[i] corresponds to order (n_min + i).
    pub coefficients: Vec<Complex64>,
    /// Minimum order (negative for principal part).
    pub n_min: i32,
    /// Maximum order.
    pub n_max: i32,
    /// Inner radius of annulus of convergence.
    pub inner_radius: f64,
    /// Outer radius of annulus of convergence.
    pub outer_radius: f64,
}

impl LaurentSeries {
    /// Compute Laurent coefficients numerically.
    /// a_n = (1/2πi) ∮_C f(z)/(z-c)^(n+1) dz for all n.
    pub fn compute<F>(
        f: F,
        center: Complex64,
        n_min: i32,
        n_max: i32,
        inner_radius: f64,
        outer_radius: f64,
    ) -> Self
    where
        F: Fn(Complex64) -> Complex64,
    {
        let mid_radius = (inner_radius + outer_radius) / 2.0;
        let m = 5000;
        let dt = 2.0 * std::f64::consts::PI / m as f64;

        let coefficients: Vec<Complex64> = (n_min..=n_max)
            .map(|n| {
                let np1 = (n + 1) as f64;
                let mut sum = Complex64::new(0.0, 0.0);
                for k in 0..m {
                    let theta = k as f64 * dt;
                    let z = center + mid_radius * Complex64::from_polar(1.0, theta);
                    let dz = mid_radius * Complex64::from_polar(1.0, theta) * Complex64::i();
                    let denom = (z - center).powf(np1);
                    sum += f(z) / denom * dz * dt;
                }
                sum / (2.0 * std::f64::consts::PI * Complex64::i())
            })
            .collect();

        LaurentSeries {
            center,
            coefficients,
            n_min,
            n_max,
            inner_radius,
            outer_radius,
        }
    }

    /// Evaluate the Laurent series at a point.
    pub fn evaluate(&self, z: Complex64) -> Complex64 {
        let diff = z - self.center;
        let mut result = Complex64::new(0.0, 0.0);
        for (i, coeff) in self.coefficients.iter().enumerate() {
            let n = self.n_min + i as i32;
            result += coeff * diff.powf(n as f64);
        }
        result
    }

    /// Get the principal part (terms with negative powers).
    pub fn principal_part(&self) -> Vec<(i32, Complex64)> {
        self.coefficients
            .iter()
            .enumerate()
            .filter(|(i, _)| (self.n_min + *i as i32) < 0)
            .map(|(i, c)| (self.n_min + i as i32, *c))
            .collect()
    }

    /// Get the coefficient a_n.
    pub fn coefficient(&self, n: i32) -> Option<Complex64> {
        let idx = n - self.n_min;
        if idx >= 0 && (idx as usize) < self.coefficients.len() {
            Some(self.coefficients[idx as usize])
        } else {
            None
        }
    }
}
