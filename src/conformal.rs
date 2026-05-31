//! Conformal mapping: Möbius transformations.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use crate::complex::ComplexExt;

/// Möbius transformation: T(z) = (az + b) / (cz + d).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MobiusTransformation {
    pub a: Complex64,
    pub b: Complex64,
    pub c: Complex64,
    pub d: Complex64,
}

impl MobiusTransformation {
    /// Create a new Möbius transformation T(z) = (az + b) / (cz + d).
    pub fn new(a: Complex64, b: Complex64, c: Complex64, d: Complex64) -> Self {
        // Normalize so that ad - bc != 0
        let det = a * d - b * c;
        assert!(det.norm() > 1e-15, "Möbius transformation must be non-degenerate (ad - bc != 0)");
        Self { a, b, c, d }
    }

    /// Evaluate T(z).
    pub fn apply(&self, z: Complex64) -> Complex64 {
        (self.a * z + self.b) / (self.c * z + self.d)
    }

    /// Compute the inverse transformation T^(-1).
    pub fn inverse(&self) -> Self {
        let det = self.a * self.d - self.b * self.c;
        MobiusTransformation::new(
            self.d / det,
            -self.b / det,
            -self.c / det,
            self.a / det,
        )
    }

    /// Compose two Möbius transformations: T1 ∘ T2.
    pub fn compose(&self, other: &MobiusTransformation) -> Self {
        MobiusTransformation::new(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d,
        )
    }

    /// Identity transformation: T(z) = z.
    pub fn identity() -> Self {
        Self::new(
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        )
    }

    /// Translation: T(z) = z + w.
    pub fn translation(w: Complex64) -> Self {
        Self::new(
            Complex64::new(1.0, 0.0),
            w,
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        )
    }

    /// Scaling: T(z) = λz.
    pub fn scaling(lambda: Complex64) -> Self {
        Self::new(
            lambda,
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        )
    }

    /// Rotation: T(z) = e^(iθ) z.
    pub fn rotation(theta: f64) -> Self {
        Self::scaling(Complex64::from_polar(1.0, theta))
    }

    /// Inversion: T(z) = 1/z.
    pub fn inversion() -> Self {
        Self::new(
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
        )
    }

    /// Map three points z1, z2, z3 to w1, w2, w3.
    pub fn from_three_points(
        z1: Complex64, z2: Complex64, z3: Complex64,
        w1: Complex64, w2: Complex64, w3: Complex64,
    ) -> Self {
        // Cross-ratio approach: (z - z1)(z2 - z3) / ((z - z3)(z2 - z1)) = (w - w1)(w2 - w3) / ((w - w3)(w2 - w1))
        // We build T(z) that maps z -> w via the cross-ratio equality.

        // T(z) = (Az + B) / (Cz + D) where we solve for A,B,C,D
        // Using the standard formula for cross-ratio mapping:
        let cross_z = |z: Complex64| (z - z1) * (z2 - z3) / ((z - z3) * (z2 - z1));
        let cross_w = |w: Complex64| (w - w1) * (w2 - w3) / ((w - w3) * (w2 - w1));

        // We need T such that cross_w(T(z)) = cross_z(z)
        // This gives us: T = cross_w^(-1) ∘ cross_z
        // cross_w(w) = q means w = w3*(w2-w1)*q - w1*(w2-w3) / ((w2-w1)*q - (w2-w3))
        // Let's just solve the linear system from 3 point conditions.

        // T(zi) = wi for i=1,2,3
        // (a*z1 + b) = w1*(c*z1 + d)
        // (a*z2 + b) = w2*(c*z2 + d)
        // (a*z3 + b) = w3*(c*z3 + d)
        //
        // Set d = 1 (normalize), solve for a, b, c.

        // From condition 1: a*z1 + b = w1*c*z1 + w1
        // From condition 2: a*z2 + b = w2*c*z2 + w2
        // From condition 3: a*z3 + b = w3*c*z3 + w3

        // Subtracting: a*(z2-z1) + 0 = c*(w2*z2 - w1*z1) + (w2 - w1)
        // And:        a*(z3-z1) + 0 = c*(w3*z3 - w1*z1) + (w3 - w1)

        // 2 equations, 2 unknowns (a, c), then b from eq1.

        use std::f64;
        let d = Complex64::new(1.0, 0.0);

        // a*(z2-z1) = c*(w2*z2 - w1*z1) + (w2-w1)
        // a*(z3-z1) = c*(w3*z3 - w1*z1) + (w3-w1)

        let pz = z2 - z1;
        let qz = z3 - z1;
        let rc = w2 * z2 - w1 * z1;
        let sc = w3 * z3 - w1 * z1;
        let rw = w2 - w1;
        let sw = w3 - w1;

        // a*pz - c*rc = rw
        // a*qz - c*sc = sw
        // Cramer's rule:
        let det = pz * sc - qz * rc;
        if det.norm() < 1e-15 {
            // Degenerate case
            return Self::identity();
        }
        let a = (rw * sc - sw * rc) / det;
        let c = (rw * qz - sw * pz) / (-det);
        let b = w1 * c * z1 + w1 - a * z1;

        Self::new(a, b, c, d)
    }

    /// Compute the determinant ad - bc.
    pub fn determinant(&self) -> Complex64 {
        self.a * self.d - self.b * self.c
    }

    /// Find the fixed points of the transformation.
    pub fn fixed_points(&self) -> (Option<Complex64>, Option<Complex64>) {
        // T(z) = z => az + b = cz^2 + dz => cz^2 + (d-a)z - b = 0
        let a_coeff = self.c;
        let b_coeff = self.d - self.a;
        let c_coeff = -self.b;

        if a_coeff.norm() < 1e-15 {
            // Linear: b_coeff * z = -c_coeff => z = -c_coeff / b_coeff
            if b_coeff.norm() < 1e-15 {
                return (None, None);
            }
            return (Some(-c_coeff / b_coeff), None);
        }

        let disc = b_coeff * b_coeff - 4.0 * a_coeff * c_coeff;
        let sqrt_disc = disc.powc(Complex64::new(0.5, 0.0));
        let z1 = (-b_coeff + sqrt_disc) / (2.0 * a_coeff);
        let z2 = (-b_coeff - sqrt_disc) / (2.0 * a_coeff);
        (Some(z1), Some(z2))
    }

    /// Check if the transformation is elliptic (conjugate to a rotation).
    pub fn is_elliptic(&self) -> bool {
        let tr = self.a + self.d;
        let det = self.determinant();
        let sigma = tr * tr / det;
        sigma.re > 0.0 && sigma.im.abs() < 1e-10 && sigma.re < 4.0
    }

    /// Apply to a batch of points.
    pub fn apply_batch(&self, points: &[Complex64]) -> Vec<Complex64> {
        points.iter().map(|&z| self.apply(z)).collect()
    }
}
