//! # lau-complex-analysis
//!
//! Complex analysis library: holomorphic functions, contour integration,
//! residue theory, conformal mapping, and agent frequency analysis.

pub mod complex;
pub mod holomorphic;
pub mod integration;
pub mod cauchy;
pub mod series;
pub mod residue;
pub mod argument;
pub mod conformal;
pub mod agent_frequency;

pub use complex::ComplexExt;
pub use holomorphic::HolomorphicCheck;
pub use integration::{Contour, ContourIntegrator};
pub use cauchy::CauchyIntegral;
pub use series::{TaylorSeries, LaurentSeries};
pub use residue::{Residue, ResidueTheorem};
pub use argument::{ArgumentPrinciple, RoucheTheorem};
pub use conformal::MobiusTransformation;
pub use agent_frequency::AgentFrequencyAnalysis;
