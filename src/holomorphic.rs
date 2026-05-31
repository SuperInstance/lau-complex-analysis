//! Holomorphic function checking via Cauchy-Riemann equations.

use num_complex::Complex64;
use crate::complex::ComplexExt;

/// Result of checking holomorphicity at a point.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct HolomorphicResult {
    pub point: (f64, f64),
    pub u_x: f64,
    pub u_y: f64,
    pub v_x: f64,
    pub v_y: f64,
    pub cr1_error: f64,
    pub cr2_error: f64,
    pub is_holomorphic: bool,
    pub tolerance: f64,
}

/// Check holomorphicity using Cauchy-Riemann equations.
pub struct HolomorphicCheck {
    pub tolerance: f64,
}

impl Default for HolomorphicCheck {
    fn default() -> Self {
        Self { tolerance: 1e-6 }
    }
}

impl HolomorphicCheck {
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Check if f(z) = u(x,y) + i*v(x,y) is holomorphic at a point.
    /// f takes a complex number and returns a complex number.
    /// Uses numerical differentiation.
    pub fn check_at<F>(&self, f: F, x: f64, y: f64) -> HolomorphicResult
    where
        F: Fn(Complex64) -> Complex64,
    {
        let h = 1e-7;
        let z = Complex64::new(x, y);
        let z_px = Complex64::new(x + h, y);
        let z_py = Complex64::new(x, y + h);

        let fz = f(z);
        let fpx = f(z_px);
        let fpy = f(z_py);

        // u(x,y) = Re(f(z)), v(x,y) = Im(f(z))
        let u_x = (fpx.re - fz.re) / h;
        let u_y = (fpy.re - fz.re) / h;
        let v_x = (fpx.im - fz.im) / h;
        let v_y = (fpy.im - fz.im) / h;

        // Cauchy-Riemann: u_x = v_y and u_y = -v_x
        let cr1_error = (u_x - v_y).abs();
        let cr2_error = (u_y + v_x).abs();

        HolomorphicResult {
            point: (x, y),
            u_x,
            u_y,
            v_x,
            v_y,
            cr1_error,
            cr2_error,
            is_holomorphic: cr1_error < self.tolerance && cr2_error < self.tolerance,
            tolerance: self.tolerance,
        }
    }

    /// Check holomorphicity at multiple points.
    pub fn check_region<F>(&self, f: F, points: &[(f64, f64)]) -> Vec<HolomorphicResult>
    where
        F: Fn(Complex64) -> Complex64,
    {
        points.iter().map(|&(x, y)| self.check_at(&f, x, y)).collect()
    }

    /// Check if f is entire (holomorphic everywhere) by sampling a grid.
    pub fn is_entire<F>(&self, f: F, range: (f64, f64), steps: usize) -> bool
    where
        F: Fn(Complex64) -> Complex64,
    {
        let (lo, hi) = range;
        for i in 0..=steps {
            for j in 0..=steps {
                let x = lo + (hi - lo) * i as f64 / steps as f64;
                let y = lo + (hi - lo) * j as f64 / steps as f64;
                if !self.check_at(&f, x, y).is_holomorphic {
                    return false;
                }
            }
        }
        true
    }

    /// Compute the derivative f'(z) numerically if f is holomorphic.
    pub fn derivative<F>(f: &F, z: Complex64) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let h = 1e-8;
        (f(z + h) - f(z - h)) / (2.0 * h)
    }
}
