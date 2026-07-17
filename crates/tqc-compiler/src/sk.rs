//! Phase weaving over the Atlas generator alphabet.
//!
//! Two regimes, matching the exact density decision:
//!
//! - **Clifford (single-handle) carrier.** Gate-set density is exactly refuted over
//!   `Q(ζ₂₄)` (finite projective Clifford image, order 24), so no dense epsilon-net exists;
//!   the weaver admits a rotation only when an exact discrete generator phase lies within
//!   the caller's tolerance, and errors otherwise (never approximates).
//! - **Density-certified carrier** (`Certified22` / `Certified576`). The projective closure
//!   is dense (PU(22) / PU(576)), so arbitrary rotations are realizable. The weaver
//!   synthesizes a rotation by a **deterministic bounded search** over spectral-flow powers:
//!   the relative phase of the flow between two eigenspaces of integer eigenvalue gap `d` is
//!   `k·d` for integer `k`, and `{k·d \bmod 2\pi}` is dense (as `d/2\pi` is irrational), so
//!   some bounded `k` lands within `ε` of the target. The emitted word is `Φ^k`.

use crate::BraidGen;

/// The carrier the weaver targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Carrier {
    /// The single-handle Clifford carrier (gate-set density refuted).
    Clifford,
    /// The 22-dimensional density-certified carrier (PU(22)).
    Certified22,
    /// The 576-dimensional two-handle density-certified carrier (PU(576)).
    Certified576,
}

/// A net node: a short generator word and the exact discrete phase it realizes.
#[derive(Debug, Clone)]
pub struct NetNode {
    /// The discrete topological braid word.
    pub word: Vec<BraidGen>,
    /// The exact phase angle realized by this word.
    pub phase_angle: f64,
}

/// The phase weaver over the Atlas generator alphabet.
pub struct SkWeaver {
    phase_net: Vec<NetNode>,
    carrier: Carrier,
    /// The integer eigenvalue gap of the spectral-flow generator on a certified carrier;
    /// the relative flow phase over a step is `flow_gap` radians.
    flow_gap: i64,
    /// The bounded search depth for certified-carrier synthesis.
    flow_search_bound: i64,
}

impl SkWeaver {
    /// Initializes the finite net of exact discrete generator phases (Clifford carrier),
    /// derived parametrically from the generator orders (σ: `2π/scope`, τ: `2π/context`,
    /// μ: `π`).
    #[must_use]
    pub fn new(p: &tqc_core::UseCaseParams) -> Self {
        let mut w = Self {
            phase_net: Vec::new(),
            carrier: Carrier::Clifford,
            flow_gap: 0,
            flow_search_bound: 0,
        };
        w.build_discrete_net(p);
        w
    }

    /// Initializes the weaver on a density-certified carrier. The spectral-flow gap is a
    /// nonzero eigenvalue gap of the balanced operator `M` (eigenvalues `{10, 7, 2, -1}`);
    /// the gap `2 - (-1) = 3` is used. `{k·3 \bmod 2\pi}` is dense, so bounded search
    /// realizes any target phase.
    #[must_use]
    pub fn for_carrier(p: &tqc_core::UseCaseParams, carrier: Carrier) -> Self {
        let mut w = Self {
            phase_net: Vec::new(),
            carrier,
            flow_gap: 3,
            flow_search_bound: 1_000_000,
        };
        w.build_discrete_net(p);
        w
    }

    fn build_discrete_net(&mut self, p: &tqc_core::UseCaseParams) {
        let pi = core::f64::consts::PI;
        let sigma_phase = 2.0 * pi / (p.scope as f64);
        let tau_phase = 2.0 * pi / (p.context as f64);
        let mu_phase = pi;
        self.phase_net = vec![
            NetNode {
                word: vec![],
                phase_angle: 0.0,
            },
            NetNode {
                word: vec![BraidGen::Sigma],
                phase_angle: sigma_phase,
            },
            NetNode {
                word: vec![BraidGen::Tau],
                phase_angle: tau_phase,
            },
            NetNode {
                word: vec![BraidGen::Sigma, BraidGen::Sigma],
                phase_angle: sigma_phase * 2.0,
            },
            NetNode {
                word: vec![BraidGen::Tau, BraidGen::Tau],
                phase_angle: tau_phase * 2.0,
            },
            NetNode {
                word: vec![BraidGen::Sigma, BraidGen::Tau, BraidGen::Mu],
                phase_angle: sigma_phase + tau_phase + mu_phase,
            },
            NetNode {
                word: vec![BraidGen::Mu, BraidGen::Sigma, BraidGen::Tau],
                phase_angle: mu_phase + sigma_phase + tau_phase,
            },
        ];
    }

    /// Circle distance between two angles, in `[0, π]`.
    fn circle_dist(a: f64, b: f64) -> f64 {
        let two_pi = 2.0 * core::f64::consts::PI;
        let mut d = (a - b).rem_euclid(two_pi);
        if d > core::f64::consts::PI {
            d = two_pi - d;
        }
        d
    }

    /// Synthesizes a rotation `theta` to precision `epsilon`.
    ///
    /// On the Clifford carrier: returns the word of an exact discrete phase within `epsilon`,
    /// or an error (no approximation path exists).
    ///
    /// On a certified carrier: returns `Φ^k` (a bounded spectral-flow power) whose relative
    /// phase `k·flow_gap` lies within `epsilon` of `theta` — a deterministic search.
    ///
    /// # Errors
    /// If no admissible word within `epsilon` is found.
    pub fn synthesize_rotation(&self, theta: f64, epsilon: f64) -> Result<Vec<BraidGen>, String> {
        let two_pi = 2.0 * core::f64::consts::PI;
        let target = theta.rem_euclid(two_pi);

        match self.carrier {
            Carrier::Clifford => {
                for node in &self.phase_net {
                    if Self::circle_dist(node.phase_angle, target) < epsilon {
                        return Ok(node.word.clone());
                    }
                }
                Err(format!(
                    "no exact discrete generator phase lies within {epsilon} of {theta}; \
                     single-handle gate-set density is exactly refuted over Q(zeta_24) \
                     (finite projective Clifford image, order 24), so arbitrary rotations \
                     cannot be approximated"
                ))
            }
            Carrier::Certified22 | Carrier::Certified576 => {
                // Deterministic bounded search over flow powers k = 0..bound.
                let d = self.flow_gap as f64;
                let mut best_k = 0i64;
                let mut best = Self::circle_dist(0.0, target);
                for k in 0..=self.flow_search_bound {
                    let phase = (k as f64 * d).rem_euclid(two_pi);
                    let dist = Self::circle_dist(phase, target);
                    if dist < best {
                        best = dist;
                        best_k = k;
                        if best < epsilon {
                            break;
                        }
                    }
                }
                if best < epsilon {
                    Ok(vec![BraidGen::Flow; best_k as usize])
                } else {
                    Err(format!(
                        "certified-carrier spectral-flow search did not reach epsilon {epsilon} \
                         for theta {theta} within {} steps (nearest {best})",
                        self.flow_search_bound
                    ))
                }
            }
        }
    }

    /// The exact relative phase realized by `k` spectral-flow steps: `k·flow_gap` radians.
    /// Deterministic; used to verify a synthesized certified-carrier word.
    #[must_use]
    pub fn flow_phase(&self, k: usize) -> f64 {
        (k as f64 * self.flow_gap as f64).rem_euclid(2.0 * core::f64::consts::PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tqc_core::UseCaseParams;

    #[test]
    fn clifford_admits_exact_and_rejects_arbitrary() {
        let p = UseCaseParams::new(4, 3, 8);
        let w = SkWeaver::new(&p);
        assert_eq!(
            w.synthesize_rotation(core::f64::consts::PI / 2.0, 1e-9)
                .unwrap(),
            vec![BraidGen::Sigma]
        );
        assert!(w.synthesize_rotation(1.0, 1e-9).is_err());
    }

    #[test]
    fn certified_synthesizes_arbitrary_rotation_within_epsilon() {
        let p = UseCaseParams::new(4, 3, 8);
        let w = SkWeaver::for_carrier(&p, Carrier::Certified22);
        // An arbitrary irrational-looking target that has no exact discrete phase.
        for &theta in &[0.3f64, 1.0, 2.5, core::f64::consts::PI / 8.0] {
            let word = w.synthesize_rotation(theta, 0.05).unwrap();
            let k = word.len();
            assert!(word.iter().all(|g| *g == BraidGen::Flow));
            let phase = w.flow_phase(k);
            assert!(
                SkWeaver::circle_dist(phase, theta.rem_euclid(2.0 * core::f64::consts::PI)) < 0.05
            );
        }
    }
}
