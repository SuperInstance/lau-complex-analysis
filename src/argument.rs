//! Argument principle and Rouché's theorem.

use num_complex::Complex64;
use crate::integration::{Contour, ContourIntegrator};
use serde::{Deserialize, Serialize};

/// Argument principle: (1/2πi) ∮_C f'(z)/f(z) dz = Z - P
/// where Z = zeros, P = poles inside C.
pub struct ArgumentPrinciple;

impl ArgumentPrinciple {
    /// Compute the winding number of f(C) around the origin,
    /// i.e., Z - P (zeros minus poles inside C).
    pub fn winding_number<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
    where
        F: Fn(Complex64) -> Complex64,
        G: Fn(Complex64) -> Complex64,
    {
        let integrator = ContourIntegrator::new(50000);
        let integrand = |z: Complex64| f_prime(z) / f(z);
        let result = integrator.integrate(integrand, contour);
        let winding = result / (2.0 * std::f64::consts::PI * Complex64::i());
        winding.re.round() as i32
    }

    /// Count zeros of f inside C (assuming f has no poles).
    pub fn count_zeros<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
    where
        F: Fn(Complex64) -> Complex64,
        G: Fn(Complex64) -> Complex64,
    {
        Self::winding_number(f, f_prime, contour)
    }

    /// Count zeros minus poles inside C.
    pub fn zeros_minus_poles<F, G>(f: F, f_prime: G, contour: &Contour) -> i32
    where
        F: Fn(Complex64) -> Complex64,
        G: Fn(Complex64) -> Complex64,
    {
        Self::winding_number(f, f_prime, contour)
    }
}

/// Rouché's theorem: if |f(z) - g(z)| < |f(z)| on C, then f and g have
/// the same number of zeros inside C.
pub struct RoucheTheorem;

impl RoucheTheorem {
    /// Check if Rouché's condition holds on the contour.
    /// Returns (condition_holds, max|f-g|, min|f|).
    pub fn check_condition<F, G>(
        f: F,
        g: G,
        contour: &Contour,
        n_samples: usize,
    ) -> RoucheResult
    where
        F: Fn(Complex64) -> Complex64,
        G: Fn(Complex64) -> Complex64,
    {
        let dt = 1.0 / n_samples as f64;
        let mut max_diff = 0.0f64;
        let mut min_f = f64::INFINITY;

        for k in 0..n_samples {
            let t = k as f64 * dt;
            let (z, _) = contour.evaluate(t);
            let fz = f(z);
            let gz = g(z);
            let diff = (fz - gz).norm();
            let fnorm = fz.norm();
            if diff > max_diff {
                max_diff = diff;
            }
            if fnorm < min_f {
                min_f = fnorm;
            }
        }

        RoucheResult {
            condition_holds: max_diff < min_f,
            max_difference: max_diff,
            min_f_norm: min_f,
        }
    }

    /// If Rouché condition holds, return the number of zeros of f inside C.
    pub fn count_zeros<F, G, H>(
        f: F,
        f_prime: G,
        _g: H,
        contour: &Contour,
        n_samples: usize,
    ) -> Option<i32>
    where
        F: Fn(Complex64) -> Complex64 + Clone,
        G: Fn(Complex64) -> Complex64,
        H: Fn(Complex64) -> Complex64,
    {
        let result = Self::check_condition(f.clone(), _g, contour, n_samples);
        if result.condition_holds {
            Some(ArgumentPrinciple::count_zeros(f, f_prime, contour))
        } else {
            None
        }
    }
}

/// Result of Rouché condition check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RoucheResult {
    pub condition_holds: bool,
    pub max_difference: f64,
    pub min_f_norm: f64,
}
