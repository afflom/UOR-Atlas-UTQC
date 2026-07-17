//! Discrete-phase weaving over the finite generator phase set.
//!
//! The exact `Q(ζ₂₄)` certificate refutes single-handle gate-set density (finite projective
//! Clifford image, order 24), so no dense epsilon-net exists and Solovay–Kitaev
//! approximation of arbitrary rotations is mathematically unavailable here. The weaver
//! therefore admits a rotation only when an **exact discrete generator phase** lies within
//! the caller's tolerance; anything else is an error, never an approximation.

use crate::BraidGen;

/// A net node: a short generator word and the exact discrete phase it realizes.
#[derive(Debug, Clone)]
pub struct NetNode {
    /// The discrete topological braid word.
    pub word: Vec<BraidGen>,
    /// The exact phase angle realized by this word.
    pub phase_angle: f64,
}

/// The discrete-phase weaver over the finite generator phase set.
pub struct SkWeaver {
    phase_net: Vec<NetNode>,
}

impl SkWeaver {
    /// Initializes the finite net of exact discrete generator phases, derived
    /// parametrically from the generator orders (σ: `2π/scope`, τ: `2π/context`, μ: `π`).
    #[must_use]
    pub fn new(p: &tqc_core::UseCaseParams) -> Self {
        let mut net = Vec::new();
        let pi = core::f64::consts::PI;

        // The identity word: the exact phase 0.
        net.push(NetNode {
            word: vec![],
            phase_angle: 0.0,
        });

        let sigma_phase = 2.0 * pi / (p.scope as f64);
        let tau_phase = 2.0 * pi / (p.context as f64);
        let mu_phase = pi;

        net.push(NetNode {
            word: vec![BraidGen::Sigma],
            phase_angle: sigma_phase,
        });
        net.push(NetNode {
            word: vec![BraidGen::Tau],
            phase_angle: tau_phase,
        });
        net.push(NetNode {
            word: vec![BraidGen::Sigma, BraidGen::Sigma],
            phase_angle: sigma_phase * 2.0,
        });
        net.push(NetNode {
            word: vec![BraidGen::Tau, BraidGen::Tau],
            phase_angle: tau_phase * 2.0,
        });
        net.push(NetNode {
            word: vec![BraidGen::Sigma, BraidGen::Tau, BraidGen::Mu],
            phase_angle: sigma_phase + tau_phase + mu_phase,
        });
        net.push(NetNode {
            word: vec![BraidGen::Mu, BraidGen::Sigma, BraidGen::Tau],
            phase_angle: mu_phase + sigma_phase + tau_phase,
        });

        Self { phase_net: net }
    }

    /// Admits a rotation `theta` only if an exact discrete generator phase lies within
    /// `epsilon` of it (on the circle), returning that word. There is no approximation
    /// path: gate-set density is exactly refuted, so Solovay–Kitaev synthesis of arbitrary
    /// rotations does not exist for this generator set.
    ///
    /// # Errors
    /// If no exact discrete phase lies within `epsilon` of `theta`.
    pub fn synthesize_rotation(&self, theta: f64, epsilon: f64) -> Result<Vec<BraidGen>, String> {
        let two_pi = 2.0 * core::f64::consts::PI;
        let mut target = theta % two_pi;
        if target < 0.0 {
            target += two_pi;
        }

        for node in &self.phase_net {
            let diff = (node.phase_angle - target).abs();
            let gap = diff.min(two_pi - diff);
            if gap < epsilon {
                return Ok(node.word.clone());
            }
        }
        Err(format!(
            "no exact discrete generator phase lies within {epsilon} of {theta}; single-handle \
             gate-set density is exactly refuted over Q(zeta_24) (finite projective Clifford \
             image, order 24), so arbitrary rotations cannot be approximated"
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tqc_core::UseCaseParams;

    #[test]
    fn exact_discrete_phases_are_admitted() {
        let p = UseCaseParams::new(4, 3, 8);
        let w = SkWeaver::new(&p);
        // 2π/4 is the exact σ phase.
        let word = w
            .synthesize_rotation(core::f64::consts::PI / 2.0, 1e-9)
            .unwrap();
        assert_eq!(word, vec![BraidGen::Sigma]);
        // The zero rotation is the identity word.
        assert!(w.synthesize_rotation(0.0, 1e-9).unwrap().is_empty());
    }

    #[test]
    fn arbitrary_rotations_are_rejected_not_approximated() {
        let p = UseCaseParams::new(4, 3, 8);
        let w = SkWeaver::new(&p);
        // 1.0 rad is no exact discrete phase; with a tight tolerance it must be an error.
        assert!(w.synthesize_rotation(1.0, 1e-9).is_err());
    }
}
