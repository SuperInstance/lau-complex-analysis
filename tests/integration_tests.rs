//! Tests for lau-complex-analysis.

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use num_complex::Complex64;
    use lau_complex_analysis::complex::ComplexExt;
    use lau_complex_analysis::holomorphic::HolomorphicCheck;
    use lau_complex_analysis::integration::{Contour, ContourIntegrator};
    use lau_complex_analysis::cauchy::CauchyIntegral;
    use lau_complex_analysis::series::{TaylorSeries, LaurentSeries};
    use lau_complex_analysis::residue::{Residue, ResidueTheorem, SingularityType};
    use lau_complex_analysis::argument::{ArgumentPrinciple, RoucheTheorem};
    use lau_complex_analysis::conformal::MobiusTransformation;
    use lau_complex_analysis::agent_frequency::{AgentFrequencyAnalysis, ZTransform};

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    // === Complex Arithmetic Tests ===

    #[test]
    fn test_polar_roundtrip() {
        let z = c(3.0, 4.0);
        let (r, theta) = z.to_polar();
        assert_relative_eq!(r, 5.0, epsilon = 1e-10);
        let z2 = Complex64::from_polar(r, theta);
        assert_relative_eq!(z2.re, 3.0, epsilon = 1e-10);
        assert_relative_eq!(z2.im, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_exp() {
        let z = c(0.0, std::f64::consts::PI);
        let result = z.exp();
        assert_relative_eq!(result.re, -1.0, epsilon = 1e-10);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_ln() {
        let z = c(-1.0, 0.0);
        let result = z.ln();
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(result.im, std::f64::consts::PI, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_power() {
        // i^i = e^(i*ln(i)) = e^(i*i*pi/2) = e^(-pi/2) ≈ 0.2079
        let i = c(0.0, 1.0);
        let result = i.powc(c(0.0, 1.0));
        assert_relative_eq!(result.re, std::f64::consts::E.powf(-std::f64::consts::FRAC_PI_2), epsilon = 1e-10);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_nth_root() {
        let z = c(1.0, 0.0);
        let roots = z.all_nth_roots(4);
        assert_eq!(roots.len(), 4);
        // Each root should satisfy root^4 ≈ 1
        for root in &roots {
            let r4 = root.powc(c(4.0, 0.0));
            assert_relative_eq!(r4.re, 1.0, epsilon = 1e-10);
            assert_relative_eq!(r4.im, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_complex_sin() {
        let z = c(0.0, 0.0);
        let result = z.csin();
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_cos() {
        let z = c(0.0, 0.0);
        let result = z.ccos();
        assert_relative_eq!(result.re, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_conjugate() {
        let z = c(3.0, -4.0);
        let conj = z.conj();
        assert_relative_eq!(conj.re, 3.0, epsilon = 1e-10);
        assert_relative_eq!(conj.im, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_from_polar() {
        let z = Complex64::from_polar(1.0, std::f64::consts::FRAC_PI_2);
        assert_relative_eq!(z.re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(z.im, 1.0, epsilon = 1e-10);
    }

    // === Holomorphic Function Tests ===

    #[test]
    fn test_cauchy_riemann_exp() {
        let check = HolomorphicCheck::new(1e-4);
        // f(z) = e^z is entire
        let result = check.check_at(|z| z.exp(), 1.0, 1.0);
        assert!(result.is_holomorphic);
    }

    #[test]
    fn test_cauchy_riemann_z_squared() {
        let check = HolomorphicCheck::new(1e-4);
        // f(z) = z^2 is entire
        let result = check.check_at(|z| z * z, 2.0, -1.0);
        assert!(result.is_holomorphic);
    }

    #[test]
    fn test_not_holomorphic_conjugate() {
        let check = HolomorphicCheck::new(1e-4);
        // f(z) = conj(z) is NOT holomorphic
        let result = check.check_at(|z| z.conj(), 1.0, 1.0);
        assert!(!result.is_holomorphic);
    }

    #[test]
    fn test_is_entire_exp() {
        let check = HolomorphicCheck::new(1e-3);
        assert!(check.is_entire(|z| z.exp(), (-2.0, 2.0), 5));
    }

    #[test]
    fn test_derivative_numerical() {
        let f = |z: Complex64| z * z;
        let deriv = HolomorphicCheck::derivative(&f, c(2.0, 0.0));
        assert_relative_eq!(deriv.re, 4.0, epsilon = 1e-5);
        assert_relative_eq!(deriv.im, 0.0, epsilon = 1e-5);
    }

    // === Integration Tests ===

    #[test]
    fn test_contour_line_integral() {
        let integrator = ContourIntegrator::new(10000);
        let contour = Contour::LineSegment { start: (0.0, 0.0), end: (1.0, 0.0) };
        // ∫_0^1 z dz = z^2/2 |_0^1 = 1/2
        let result = integrator.integrate(|z| z, &contour);
        assert_relative_eq!(result.re, 0.5, epsilon = 1e-4);
    }

    #[test]
    fn test_circle_integral_z() {
        let integrator = ContourIntegrator::new(50000);
        let contour = Contour::unit_circle();
        // ∮_|z|=1 z dz = 0 (holomorphic inside)
        let result = integrator.integrate(|z| z, &contour);
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-3);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_circle_integral_1_over_z() {
        let integrator = ContourIntegrator::new(50000);
        let contour = Contour::unit_circle();
        // ∮_|z|=1 1/z dz = 2πi
        let result = integrator.integrate(|z| c(1.0, 0.0) / z, &contour);
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-2);
        assert_relative_eq!(result.im, 2.0 * std::f64::consts::PI, epsilon = 1e-2);
    }

    #[test]
    fn test_polygon_contour() {
        let contour = Contour::Polygon {
            vertices: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        };
        let integrator = ContourIntegrator::new(50000);
        // ∮ f(z)=1 around a square = 0
        let result = integrator.integrate(|_| c(1.0, 0.0), &contour);
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-3);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-3);
    }

    // === Cauchy Integral Tests ===

    #[test]
    fn test_cauchy_formula_exp() {
        let cauchy = CauchyIntegral::new(50000);
        let contour = Contour::circle((0.0, 0.0), 2.0);
        // f(z) = e^z, f(1) = e
        let result = cauchy.cauchy_formula(|z| z.exp(), &contour, c(1.0, 0.0));
        assert_relative_eq!(result.re, std::f64::consts::E, epsilon = 1e-2);
    }

    #[test]
    fn test_cauchy_theorem_holomorphic() {
        let cauchy = CauchyIntegral::new(50000);
        let contour = Contour::circle((0.0, 0.0), 1.0);
        // ∮ z^2 dz = 0
        let result = cauchy.cauchy_theorem(|z| z * z, &contour);
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-3);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_cauchy_derivative() {
        let cauchy = CauchyIntegral::new(50000);
        let contour = Contour::circle((0.0, 0.0), 2.0);
        // f(z) = z^3, f'(0) = 0, f''(0) = 0, f'''(0) = 6
        let result = cauchy.cauchy_derivative(|z| z * z * z, &contour, c(0.0, 0.0), 3);
        assert_relative_eq!(result.re, 6.0, epsilon = 0.1);
    }

    // === Series Tests ===

    #[test]
    fn test_taylor_exp() {
        let ts = TaylorSeries::compute_numerical(|z| z.exp(), c(0.0, 0.0), 8);
        // e^z = 1 + z + z^2/2 + z^3/6 + ...
        assert_relative_eq!(ts.coefficients[0].re, 1.0, epsilon = 0.05);
        assert_relative_eq!(ts.coefficients[1].re, 1.0, epsilon = 0.05);
        assert_relative_eq!(ts.coefficients[2].re, 0.5, epsilon = 0.05);
        assert_relative_eq!(ts.coefficients[3].re, 1.0 / 6.0, epsilon = 0.05);
    }

    #[test]
    fn test_taylor_evaluate() {
        let ts = TaylorSeries::compute_numerical(|z| z.exp(), c(0.0, 0.0), 10);
        let approx = ts.evaluate(c(0.5, 0.0));
        assert_relative_eq!(approx.re, (0.5f64).exp(), epsilon = 0.01);
    }

    #[test]
    fn test_laurent_simple_pole() {
        // f(z) = 1/z has Laurent series with only a_{-1} = 1
        let ls = LaurentSeries::compute(
            |z| c(1.0, 0.0) / z,
            c(0.0, 0.0),
            -3, 3,
            0.1, 2.0,
        );
        // Coefficient a_{-1} should be ~1
        let a_neg1 = ls.coefficient(-1);
        assert!(a_neg1.is_some());
        assert_relative_eq!(a_neg1.unwrap().re, 1.0, epsilon = 0.1);
        // Other coefficients should be ~0
        let a_0 = ls.coefficient(0);
        assert!(a_0.is_some());
        assert_relative_eq!(a_0.unwrap().norm(), 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_laurent_principal_part() {
        let ls = LaurentSeries::compute(
            |z| c(1.0, 0.0) / (z * z),
            c(0.0, 0.0),
            -3, 3,
            0.1, 2.0,
        );
        let pp = ls.principal_part();
        // Should have entries for n=-3,-2,-1, but only n=-2 is non-zero
        assert!(!pp.is_empty());
    }

    // === Residue Tests ===

    #[test]
    fn test_residue_simple_pole() {
        // f(z) = 1/z, Res at 0 = 1
        let res = Residue::compute(|z| c(1.0, 0.0) / z, c(0.0, 0.0), 0.5);
        assert_relative_eq!(res.re, 1.0, epsilon = 0.05);
    }

    #[test]
    fn test_residue_double_pole() {
        // f(z) = 1/z^2, Res at 0 = 0
        let res = Residue::compute(|z| c(1.0, 0.0) / (z * z), c(0.0, 0.0), 0.5);
        assert_relative_eq!(res.norm(), 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_residue_rational() {
        // f(z) = 1/(z-1), Res at 1 = 1
        let res = Residue::compute(|z| c(1.0, 0.0) / (z - c(1.0, 0.0)), c(1.0, 0.0), 0.5);
        assert_relative_eq!(res.re, 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_residue_theorem() {
        // ∮ 1/(z^2+1) around |z|=2 = 2πi(Res(i) + Res(-i)) = 2πi(1/(2i) + 1/(-2i)) = 0
        let f = |z: Complex64| c(1.0, 0.0) / (z * z + c(1.0, 0.0));
        let contour = Contour::circle((0.0, 0.0), 2.0);
        let result = ResidueTheorem::apply(f, &contour, &[c(0.0, 1.0), c(0.0, -1.0)], 0.5);
        assert_relative_eq!(result.re, 0.0, epsilon = 0.5);
        assert_relative_eq!(result.im, 0.0, epsilon = 0.5);
    }

    #[test]
    fn test_classify_simple_pole() {
        let sing = Residue::classify(|z| c(1.0, 0.0) / z, c(0.0, 0.0));
        assert_eq!(sing.singularity_type, SingularityType::SimplePole);
    }

    #[test]
    fn test_classify_removable() {
        // f(z) = sin(z)/z has a removable singularity at 0
        let sing = Residue::classify(|z| z.csin() / z, c(0.0, 0.0));
        assert_eq!(sing.singularity_type, SingularityType::Removable);
    }

    // === Argument Principle Tests ===

    #[test]
    fn test_argument_principle_zeros() {
        // f(z) = z has one zero at origin
        let contour = Contour::circle((0.0, 0.0), 0.5);
        let count = ArgumentPrinciple::count_zeros(
            |z| z,
            |_| c(1.0, 0.0),
            &contour,
        );
        assert_eq!(count, 1);
    }

    #[test]
    fn test_argument_principle_z_squared() {
        // f(z) = z^2 has a zero of order 2 at origin
        let contour = Contour::circle((0.0, 0.0), 0.5);
        let count = ArgumentPrinciple::count_zeros(
            |z| z * z,
            |z| 2.0 * z,
            &contour,
        );
        assert_eq!(count, 2);
    }

    #[test]
    fn test_rouche_condition() {
        // f(z) = z^2 + 2, g(z) = z^2 on |z| = 1
        // |f-g| = 2, |f| >= |z^2+2| >= 1 on |z|=1
        // Actually: |f-g| = 2, min |f| on |z|=1: |z^2+2| min is |1+2|=3 (at z=1)
        // Wait: let's use f(z) = z^3 + 100, g(z) = z^3 on |z|=1
        // |f-g| = 100, min |f| = min |z^3+100| = 99 (at z=-1)... still not <
        //
        // Classic example: f(z)=5z^3+z, g(z)=5z^3 on |z|=1
        // |f-g| = |z| = 1, |f| >= |5z^3| - |z| = 5-1 = 4 > 1
        let contour = Contour::circle((0.0, 0.0), 1.0);
        let result = RoucheTheorem::check_condition(
            |z| 5.0 * z * z * z + z,
            |z| 5.0 * z * z * z,
            &contour,
            1000,
        );
        assert!(result.condition_holds);
    }

    // === Möbius Transformation Tests ===

    #[test]
    fn test_mobius_identity() {
        let t = MobiusTransformation::identity();
        let z = c(3.0, -2.0);
        let result = t.apply(z);
        assert_relative_eq!(result.re, 3.0, epsilon = 1e-10);
        assert_relative_eq!(result.im, -2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mobius_inverse() {
        let t = MobiusTransformation::new(c(2.0, 0.0), c(1.0, 0.0), c(0.0, 1.0), c(1.0, 0.0));
        let t_inv = t.inverse();
        let z = c(1.0, 2.0);
        let roundtrip = t_inv.apply(t.apply(z));
        assert_relative_eq!(roundtrip.re, z.re, epsilon = 1e-10);
        assert_relative_eq!(roundtrip.im, z.im, epsilon = 1e-10);
    }

    #[test]
    fn test_mobius_compose() {
        let t1 = MobiusTransformation::translation(c(1.0, 0.0));
        let t2 = MobiusTransformation::scaling(c(2.0, 0.0));
        let comp = t2.compose(&t1);
        let z = c(3.0, 0.0);
        let result = comp.apply(z);
        // t1(z) = z+1 = 4, t2(4) = 8
        assert_relative_eq!(result.re, 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mobius_three_points() {
        let t = MobiusTransformation::from_three_points(
            c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0),
            c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0),
        );
        // Should be close to identity
        let z = c(0.5, 0.5);
        let result = t.apply(z);
        assert_relative_eq!(result.re, 0.5, epsilon = 1e-6);
        assert_relative_eq!(result.im, 0.5, epsilon = 1e-6);
    }

    #[test]
    fn test_mobius_inversion() {
        let t = MobiusTransformation::inversion();
        let z = c(2.0, 0.0);
        let result = t.apply(z);
        assert_relative_eq!(result.re, 0.5, epsilon = 1e-10);
        assert_relative_eq!(result.im, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mobius_rotation() {
        let t = MobiusTransformation::rotation(std::f64::consts::FRAC_PI_2);
        let z = c(1.0, 0.0);
        let result = t.apply(z);
        assert_relative_eq!(result.re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(result.im, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mobius_fixed_points() {
        let t = MobiusTransformation::identity();
        let (fp1, fp2) = t.fixed_points();
        // Identity has all points as fixed points (degenerate)
        // Just check it doesn't panic
        assert!(fp1.is_some() || fp2.is_some() || (fp1.is_none() && fp2.is_none()));
    }

    #[test]
    fn test_mobius_determinant() {
        let t = MobiusTransformation::new(c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let det = t.determinant();
        assert_relative_eq!(det.re, 1.0, epsilon = 1e-10);
    }

    // === Agent Frequency Analysis Tests ===

    #[test]
    fn test_z_transform_impulse() {
        let signal = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let z = c(1.0, 0.0);
        let result = ZTransform::compute(&signal, z);
        assert_relative_eq!(result.re, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_z_transform_dc() {
        let signal = vec![1.0, 1.0, 1.0, 1.0];
        let z = c(1.0, 0.0);
        let result = ZTransform::compute(&signal, z);
        assert_relative_eq!(result.re, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_frequency_analysis_stable() {
        let analysis = AgentFrequencyAnalysis::new(128);
        // Damped oscillation
        let signal: Vec<f64> = (0..64).map(|n| (0.9f64).powi(n) * (2.0 * std::f64::consts::PI * 0.1 * n as f64).cos()).collect();
        let result = analysis.analyze(&signal);
        assert!(result.spectral_radius < 1.0);
        assert!(result.is_stable);
    }

    #[test]
    fn test_frequency_analysis_unstable() {
        let analysis = AgentFrequencyAnalysis::new(128);
        // Explicitly unstable poles
        let report = analysis.stability_analysis(&[c(0.5, 0.0), c(1.2, 0.0)]);
        assert!(!report.is_stable);
    }

    #[test]
    fn test_stability_analysis() {
        let analysis = AgentFrequencyAnalysis::new(64);
        let poles = vec![c(0.5, 0.0), c(0.3, 0.4)];
        let report = analysis.stability_analysis(&poles);
        assert!(report.is_stable);
        assert!(report.spectral_radius < 1.0);
    }

    #[test]
    fn test_stability_analysis_unstable() {
        let analysis = AgentFrequencyAnalysis::new(64);
        let poles = vec![c(0.5, 0.0), c(1.5, 0.0)];
        let report = analysis.stability_analysis(&poles);
        assert!(!report.is_stable);
        assert_eq!(report.unstable_poles.len(), 1);
    }

    #[test]
    fn test_periodicity_detection() {
        let analysis = AgentFrequencyAnalysis::new(256);
        // Pure sine wave with period 10
        let signal: Vec<f64> = (0..100).map(|n| (2.0 * std::f64::consts::PI * n as f64 / 10.0).sin()).collect();
        let period = analysis.detect_periodicity(&signal);
        assert!(period.is_some());
    }

    #[test]
    fn test_magnitude_spectrum() {
        let signal: Vec<f64> = (0..32).map(|n| (2.0 * std::f64::consts::PI * n as f64 / 8.0).sin()).collect();
        let spectrum = ZTransform::magnitude_spectrum(&signal, 64);
        assert_eq!(spectrum.len(), 64);
        // Should have a peak around omega = 2*pi/8
        let peak_omega = 2.0 * std::f64::consts::PI / 8.0;
        let max_idx = spectrum.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
        assert_relative_eq!(max_idx.0, peak_omega, epsilon = 0.2);
    }

    // === Additional Integration Tests ===

    #[test]
    fn test_arc_contour() {
        let contour = Contour::Arc {
            center: (0.0, 0.0),
            radius: 1.0,
            theta0: 0.0,
            theta1: std::f64::consts::PI,
        };
        let integrator = ContourIntegrator::new(10000);
        // ∫_arc 1 dz along semicircle from 1 to -1
        let result = integrator.integrate(|_| c(1.0, 0.0), &contour);
        // z(t) = e^(i*pi*t), dz = i*pi*e^(i*pi*t) dt
        // ∫ = e^(i*pi) - 1 = -2
        assert_relative_eq!(result.re, -2.0, epsilon = 0.01);
        assert_relative_eq!(result.im, 0.0, epsilon = 0.01);
    }

    #[test]
    fn test_contour_reverse() {
        let contour = Contour::LineSegment { start: (0.0, 0.0), end: (1.0, 0.0) };
        let rev = contour.reverse();
        let (z, dz) = rev.evaluate(0.0);
        assert_relative_eq!(z.re, 1.0, epsilon = 1e-10);
    }

    // === Additional Cauchy Tests ===

    #[test]
    fn test_cauchy_formula_z_cubed() {
        let cauchy = CauchyIntegral::new(50000);
        let contour = Contour::circle((0.0, 0.0), 2.0);
        // f(z) = z^3, f(1) = 1
        let result = cauchy.cauchy_formula(|z| z * z * z, &contour, c(1.0, 0.0));
        assert_relative_eq!(result.re, 1.0, epsilon = 0.05);
    }

    // === Combined tests ===

    #[test]
    fn test_residue_integral_consistency() {
        // f(z) = e^z / z has Res = 1 at z=0
        let res = Residue::compute(|z| z.exp() / z, c(0.0, 0.0), 0.5);
        assert_relative_eq!(res.re, 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_mobius_batch() {
        let t = MobiusTransformation::scaling(c(2.0, 0.0));
        let points = vec![c(1.0, 0.0), c(0.0, 1.0), c(1.0, 1.0)];
        let results = t.apply_batch(&points);
        assert_relative_eq!(results[0].re, 2.0, epsilon = 1e-10);
        assert_relative_eq!(results[1].im, 2.0, epsilon = 1e-10);
    }
}

