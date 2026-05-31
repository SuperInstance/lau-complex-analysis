//! Agent frequency analysis using complex methods.
//!
//! Analyzes oscillatory behavior of agents using Fourier/Z-transform
//! techniques in the complex plane.

use num_complex::Complex64;
use crate::complex::ComplexExt;
use crate::series::LaurentSeries;
use serde::{Deserialize, Serialize};

/// Frequency analysis result for an agent signal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyResult {
    /// Dominant frequencies detected (in radians).
    pub frequencies: Vec<f64>,
    /// Magnitude of each frequency component.
    pub magnitudes: Vec<f64>,
    /// Phase of each frequency component.
    pub phases: Vec<f64>,
    /// Estimated spectral radius (stability indicator).
    pub spectral_radius: f64,
    /// Whether the agent behavior is stable (spectral_radius < 1).
    pub is_stable: bool,
}

/// Complex Z-transform for discrete agent signals.
pub struct ZTransform;

impl ZTransform {
    /// Compute the Z-transform of a discrete signal at a point z.
    /// X(z) = Σ x[n] * z^(-n) for n = 0..N-1.
    pub fn compute(signal: &[f64], z: Complex64) -> Complex64 {
        signal
            .iter()
            .enumerate()
            .map(|(n, &x)| Complex64::new(x, 0.0) * z.powc(Complex64::new(-(n as f64), 0.0)))
            .sum()
    }

    /// Evaluate on the unit circle (DTFT).
    pub fn dtft(signal: &[f64], omega: f64) -> Complex64 {
        let z = Complex64::from_polar(1.0, omega).powc(Complex64::new(-1.0, 0.0));
        Self::compute(signal, z)
    }

    /// Compute the magnitude spectrum on the unit circle.
    pub fn magnitude_spectrum(signal: &[f64], n_freqs: usize) -> Vec<(f64, f64)> {
        (0..n_freqs)
            .map(|k| {
                let omega = 2.0 * std::f64::consts::PI * k as f64 / n_freqs as f64;
                let xz = Self::dtft(signal, omega);
                (omega, xz.norm())
            })
            .collect()
    }
}

/// Agent frequency analysis: detect oscillatory patterns in agent behavior.
pub struct AgentFrequencyAnalysis {
    pub n_frequency_bins: usize,
}

impl Default for AgentFrequencyAnalysis {
    fn default() -> Self {
        Self { n_frequency_bins: 256 }
    }
}

impl AgentFrequencyAnalysis {
    pub fn new(n_frequency_bins: usize) -> Self {
        Self { n_frequency_bins }
    }

    /// Analyze a discrete signal for dominant frequencies and stability.
    pub fn analyze(&self, signal: &[f64]) -> FrequencyResult {
        let spectrum = ZTransform::magnitude_spectrum(signal, self.n_frequency_bins);

        // Find peaks (simple local maximum detection)
        let mut peaks: Vec<(usize, f64, f64)> = Vec::new();
        for i in 1..spectrum.len() - 1 {
            if spectrum[i].1 > spectrum[i - 1].1 && spectrum[i].1 > spectrum[i + 1].1 && spectrum[i].1 > 0.01 {
                peaks.push((i, spectrum[i].0, spectrum[i].1));
            }
        }

        // Sort by magnitude descending
        peaks.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        // Take top frequencies (up to 10)
        let top: Vec<_> = peaks.into_iter().take(10).collect();
        let frequencies: Vec<f64> = top.iter().map(|&(_, omega, _)| omega).collect();
        let magnitudes: Vec<f64> = top.iter().map(|&(_, _, mag)| mag).collect();

        // Compute phases for dominant frequencies
        let phases: Vec<f64> = frequencies
            .iter()
            .map(|&omega| ZTransform::dtft(signal, omega).arg())
            .collect();

        // Estimate spectral radius via pole analysis
        // Fit an autoregressive model and find roots
        let spectral_radius = self.estimate_spectral_radius(signal);

        FrequencyResult {
            frequencies,
            magnitudes,
            phases,
            spectral_radius,
            is_stable: spectral_radius < 1.0,
        }
    }

    /// Estimate the spectral radius using the autocorrelation method.
    fn estimate_spectral_radius(&self, signal: &[f64]) -> f64 {
        if signal.len() < 4 {
            return 0.0;
        }

        // Simple AR(2) model via Yule-Walker
        let n = signal.len() as f64;
        let mean = signal.iter().sum::<f64>() / n;

        let centered: Vec<f64> = signal.iter().map(|&x| x - mean).collect();

        // Autocorrelation at lags 0, 1, 2
        let r0: f64 = centered.iter().map(|x| x * x).sum::<f64>() / n;
        if r0.abs() < 1e-15 {
            return 0.0;
        }
        let r1: f64 = centered[..centered.len() - 1]
            .iter()
            .zip(&centered[1..])
            .map(|(a, b)| a * b)
            .sum::<f64>()
            / n;
        let r2: f64 = centered[..centered.len() - 2]
            .iter()
            .zip(&centered[2..])
            .map(|(a, b)| a * b)
            .sum::<f64>()
            / n;

        // Solve: r1 = a1*r0 + a2*r1, r2 = a1*r1 + a2*r0
        // [r0 r1] [a1]   [r1]
        // [r1 r0] [a2] = [r2]
        let det = r0 * r0 - r1 * r1;
        if det.abs() < 1e-15 {
            return r1.abs() / r0;
        }
        let a1 = (r0 * r1 - r1 * r2) / det;
        let a2 = (r0 * r2 - r1 * r1) / det;

        // Roots of z^2 - a1*z - a2 = 0
        let disc = a1 * a1 + 4.0 * a2;
        if disc >= 0.0 {
            let sqrt_disc = disc.sqrt();
            let r1_root = (a1 + sqrt_disc) / 2.0;
            let r2_root = (a1 - sqrt_disc) / 2.0;
            r1_root.abs().max(r2_root.abs())
        } else {
            // Complex roots: magnitude = sqrt(-a2)
            (-a2).sqrt().abs()
        }
    }

    /// Detect if an agent signal has periodic behavior.
    pub fn detect_periodicity(&self, signal: &[f64]) -> Option<usize> {
        let spectrum = ZTransform::magnitude_spectrum(signal, self.n_frequency_bins);

        // Find the strongest peak
        let mut best_idx = 0;
        let mut best_mag = 0.0f64;
        for (i, &(_, mag)) in spectrum.iter().enumerate() {
            if i > 0 && i < spectrum.len() - 1 && mag > best_mag {
                best_mag = mag;
                best_idx = i;
            }
        }

        if best_mag < 0.01 {
            return None;
        }

        // Convert frequency to period
        let omega = spectrum[best_idx].0;
        if omega < 1e-10 {
            return None;
        }
        let period = (2.0 * std::f64::consts::PI / omega).round() as usize;
        if period > 0 && period < signal.len() {
            Some(period)
        } else {
            None
        }
    }

    /// Compute the complex cepstrum for analyzing echo patterns.
    pub fn complex_cepstrum(&self, signal: &[f64]) -> Vec<Complex64> {
        let n = signal.len();
        (0..n)
            .map(|k| {
                let omega = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                let xz = ZTransform::dtft(signal, omega);
                if xz.norm() > 1e-15 {
                    xz.ln() / n as f64
                } else {
                    Complex64::new(0.0, 0.0)
                }
            })
            .collect()
    }

    /// Analyze stability of an agent's transfer function.
    /// Given poles of the transfer function, check if all are inside the unit circle.
    pub fn stability_analysis(&self, poles: &[Complex64]) -> StabilityReport {
        let max_radius = poles.iter().map(|p| p.norm()).fold(0.0f64, f64::max);
        let unstable_poles: Vec<Complex64> = poles.iter().filter(|p| p.norm() >= 1.0).cloned().collect();

        StabilityReport {
            poles: poles.to_vec(),
            spectral_radius: max_radius,
            is_stable: unstable_poles.is_empty(),
            unstable_poles,
            damping_ratios: poles.iter().map(|p| -p.re / p.norm()).collect(),
            natural_frequencies: poles.iter().map(|p| p.arg().abs()).collect(),
        }
    }
}

/// Stability analysis report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StabilityReport {
    pub poles: Vec<Complex64>,
    pub spectral_radius: f64,
    pub is_stable: bool,
    pub unstable_poles: Vec<Complex64>,
    pub damping_ratios: Vec<f64>,
    pub natural_frequencies: Vec<f64>,
}
