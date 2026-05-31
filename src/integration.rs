//! Contour integration: parameterized contours and line integrals.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A parameterized contour in the complex plane.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Contour {
    /// Line segment from start to end.
    LineSegment { start: (f64, f64), end: (f64, f64) },
    /// Circle centered at (cx, cy) with radius r.
    Circle { center: (f64, f64), radius: f64 },
    /// Arc from angle theta0 to theta1 on circle (cx, cy, r).
    Arc { center: (f64, f64), radius: f64, theta0: f64, theta1: f64 },
    /// Polygon (closed) defined by vertices.
    Polygon { vertices: Vec<(f64, f64)> },
    /// Generic parameterized curve: z(t) and z'(t) for t in [a,b].
    Parametric {
        // stored as sampling points for serialization
        samples: Vec<(f64, f64, f64, f64, f64)>, // (t, Re(z), Im(z), Re(z'), Im(z'))
    },
}

impl Contour {
    /// Evaluate z(t) and z'(t) at parameter t in [0,1].
    pub fn evaluate(&self, t: f64) -> (Complex64, Complex64) {
        match self {
            Contour::LineSegment { start, end } => {
                let z = Complex64::new(
                    start.0 + t * (end.0 - start.0),
                    start.1 + t * (end.1 - start.1),
                );
                let dz = Complex64::new(end.0 - start.0, end.1 - start.1);
                (z, dz)
            }
            Contour::Circle { center, radius } => {
                let theta = 2.0 * std::f64::consts::PI * t;
                let z = Complex64::new(
                    center.0 + radius * theta.cos(),
                    center.1 + radius * theta.sin(),
                );
                let dz = Complex64::new(
                    -radius * theta.sin(),
                    radius * theta.cos(),
                ) * 2.0 * std::f64::consts::PI;
                (z, dz)
            }
            Contour::Arc { center, radius, theta0, theta1 } => {
                let theta = theta0 + t * (theta1 - theta0);
                let z = Complex64::new(
                    center.0 + radius * theta.cos(),
                    center.1 + radius * theta.sin(),
                );
                let dz = Complex64::new(
                    -radius * theta.sin(),
                    radius * theta.cos(),
                ) * (theta1 - theta0);
                (z, dz)
            }
            Contour::Polygon { vertices } => {
                let n = vertices.len();
                let seg_t = t * n as f64;
                let seg_idx = (seg_t.floor() as usize).min(n - 1);
                let local_t = seg_t - seg_idx as f64;
                let start = vertices[seg_idx];
                let end = vertices[(seg_idx + 1) % n];
                let z = Complex64::new(
                    start.0 + local_t * (end.0 - start.0),
                    start.1 + local_t * (end.1 - start.1),
                );
                let dz = Complex64::new(
                    n as f64 * (end.0 - start.0),
                    n as f64 * (end.1 - start.1),
                );
                (z, dz)
            }
            Contour::Parametric { samples } => {
                // Linear interpolation between samples
                let n = samples.len().max(1) - 1;
                let idx_f = t * n as f64;
                let idx = idx_f.floor() as usize;
                let frac = idx_f - idx as f64;
                let clamped = idx.min(n);
                let next = (clamped + 1).min(n);
                let (_, zr1, zi1, dzr1, dzi1) = samples[clamped];
                let (_, zr2, zi2, dzr2, dzi2) = samples[next];
                let z = Complex64::new(
                    zr1 + frac * (zr2 - zr1),
                    zi1 + frac * (zi2 - zi1),
                );
                let dz = Complex64::new(
                    dzr1 + frac * (dzr2 - dzr1),
                    dzi1 + frac * (dzi2 - dzi1),
                );
                (z, dz)
            }
        }
    }

    /// Create a unit circle centered at origin.
    pub fn unit_circle() -> Self {
        Contour::Circle { center: (0.0, 0.0), radius: 1.0 }
    }

    /// Create a circle of given radius centered at a point.
    pub fn circle(center: (f64, f64), radius: f64) -> Self {
        Contour::Circle { center, radius }
    }

    /// Reverse orientation of the contour.
    pub fn reverse(&self) -> Self {
        match self {
            Contour::LineSegment { start, end } => Contour::LineSegment {
                start: *end,
                end: *start,
            },
            Contour::Circle { center, radius } => Contour::Arc {
                center: *center,
                radius: *radius,
                theta0: 2.0 * std::f64::consts::PI,
                theta1: 0.0,
            },
            _ => self.clone(), // for simplicity
        }
    }
}

/// Contour integrator using numerical quadrature.
pub struct ContourIntegrator {
    pub n_points: usize,
}

impl Default for ContourIntegrator {
    fn default() -> Self {
        Self { n_points: 10000 }
    }
}

impl ContourIntegrator {
    pub fn new(n_points: usize) -> Self {
        Self { n_points }
    }

    /// Compute ∫_C f(z) dz using the trapezoidal rule.
    pub fn integrate<F>(&self, f: F, contour: &Contour) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let n = self.n_points;
        let dt = 1.0 / n as f64;
        let mut sum = Complex64::new(0.0, 0.0);

        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            let (z, dz) = contour.evaluate(t);
            sum += f(z) * dz * dt;
        }

        sum
    }

    /// Compute ∫_C f(z) dz with a specified number of quadrature points.
    pub fn integrate_with_n<F>(&self, f: F, contour: &Contour, n: usize) -> Complex64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let dt = 1.0 / n as f64;
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            let (z, dz) = contour.evaluate(t);
            sum += f(z) * dz * dt;
        }
        sum
    }

    /// Compute ∫_C |f(z)| |dz| (arc-length integral of |f|).
    pub fn integrate_abs<F>(&self, f: F, contour: &Contour) -> f64
    where
        F: Fn(Complex64) -> Complex64,
    {
        let n = self.n_points;
        let dt = 1.0 / n as f64;
        let mut sum = 0.0;
        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            let (z, dz) = contour.evaluate(t);
            sum += f(z).norm() * dz.norm() * dt;
        }
        sum
    }
}
