//! Extended complex number operations: polar form, powers, roots, exp, log.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Extension trait for complex number operations.
pub trait ComplexExt {
    /// Convert to polar form (r, theta).
    fn to_polar(self) -> (f64, f64);
    /// Create from polar form.
    fn from_polar(r: f64, theta: f64) -> Complex64;
    /// Principal value of the complex logarithm.
    fn ln(self) -> Complex64;
    /// Complex logarithm with branch cut.
    fn ln_branch(self, branch: f64) -> Complex64;
    /// Complex exponential.
    fn exp(self) -> Complex64;
    /// Complex sine.
    fn csin(self) -> Complex64;
    /// Complex cosine.
    fn ccos(self) -> Complex64;
    /// Complex power: self^(w).
    fn powc(self, w: Complex64) -> Complex64;
    /// Principal nth root.
    fn nth_root(self, n: u32) -> Complex64;
    /// All nth roots.
    fn all_nth_roots(self, n: u32) -> Vec<Complex64>;
    /// Argument in (-pi, pi].
    fn arg(self) -> f64;
    /// Modulus.
    fn modulus(self) -> f64;
    /// Complex conjugate.
    fn conj(self) -> Complex64;
}

impl ComplexExt for Complex64 {
    fn to_polar(self) -> (f64, f64) {
        (self.norm(), self.arg())
    }

    fn from_polar(r: f64, theta: f64) -> Complex64 {
        Complex64::new(r * theta.cos(), r * theta.sin())
    }

    fn ln(self) -> Complex64 {
        Complex64::new(self.norm().ln(), self.arg())
    }

    fn ln_branch(self, branch: f64) -> Complex64 {
        let mut theta = self.arg();
        while theta <= branch {
            theta += 2.0 * std::f64::consts::PI;
        }
        while theta > branch + 2.0 * std::f64::consts::PI {
            theta -= 2.0 * std::f64::consts::PI;
        }
        Complex64::new(self.norm().ln(), theta)
    }

    fn exp(self) -> Complex64 {
        let r = self.re.exp();
        Complex64::new(r * self.im.cos(), r * self.im.sin())
    }

    fn csin(self) -> Complex64 {
        Complex64::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh(),
        )
    }

    fn ccos(self) -> Complex64 {
        Complex64::new(
            self.re.cos() * self.im.cosh(),
            -self.re.sin() * self.im.sinh(),
        )
    }

    fn powc(self, w: Complex64) -> Complex64 {
        if self.norm() == 0.0 {
            return Complex64::new(0.0, 0.0);
        }
        (self.ln() * w).exp()
    }

    fn nth_root(self, n: u32) -> Complex64 {
        let (r, theta) = self.to_polar();
        Complex64::from_polar(r.powf(1.0 / n as f64), theta / n as f64)
    }

    fn all_nth_roots(self, n: u32) -> Vec<Complex64> {
        let (r, theta) = self.to_polar();
        let rn = r.powf(1.0 / n as f64);
        (0..n)
            .map(|k| {
                Complex64::from_polar(
                    rn,
                    (theta + 2.0 * std::f64::consts::PI * k as f64) / n as f64,
                )
            })
            .collect()
    }

    fn arg(self) -> f64 {
        self.arg()
    }

    fn modulus(self) -> f64 {
        self.norm()
    }

    fn conj(self) -> Complex64 {
        Complex64::new(self.re, -self.im)
    }
}

/// Compute e^z.
pub fn exp(z: Complex64) -> Complex64 {
    z.exp()
}

/// Compute sin(z).
pub fn sin(z: Complex64) -> Complex64 {
    z.csin()
}

/// Compute cos(z).
pub fn cos(z: Complex64) -> Complex64 {
    z.ccos()
}

/// Compute log(z) (principal value).
pub fn log(z: Complex64) -> Complex64 {
    z.ln()
}

/// Compute z^w.
pub fn pow(z: Complex64, w: Complex64) -> Complex64 {
    z.powc(w)
}

/// Shorthand for creating complex numbers.
pub fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}
