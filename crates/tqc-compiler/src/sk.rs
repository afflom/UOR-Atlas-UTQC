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

/// A per-rotation synthesis residual certificate, shipped inside the compiled word.
///
/// For a rotation gate the certified-carrier synthesis records the Euler `Z`-flow powers
/// `(k_a, k_b, k_c)` of the `ZXZ` decomposition and the exact `f64` operator distance
/// `abs_error` between the realized `SU(2)` action (with `Φ = R_z(flow_gap)` and the Clifford
/// Hadamard template acting as the encoded Hadamard, whose exactness is the separate
/// `encoded-qubit-universality` certificate) and the target, together with the tolerance
/// `epsilon` the synthesis met (`abs_error ≤ epsilon`). This travels with the word so the
/// compilation claim and its residual are inseparable.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisResidual {
    /// Rotation axis (`x`, `y`, or `z`).
    pub axis: char,
    /// Target rotation angle (radians).
    pub target_angle: f64,
    /// The `ZXZ` Euler `Z`-flow powers `(k_a, k_b, k_c)`.
    pub euler_powers: [usize; 3],
    /// Exact `f64` operator distance between the realized and target `SU(2)` actions.
    pub abs_error: f64,
    /// The tolerance met.
    pub epsilon: f64,
}

/// A 2x2 complex matrix over `f64`, as `[[ (re, im); 2]; 2]`.
type M2 = [[(f64, f64); 2]; 2];

fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn m2mul(a: &M2, b: &M2) -> M2 {
    let mut out = [[(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                out[i][j] = cadd(out[i][j], cmul(a[i][k], b[k][j]));
            }
        }
    }
    out
}
fn expi(theta: f64) -> (f64, f64) {
    (theta.cos(), theta.sin())
}
/// `R_z(θ) = diag(e^{-iθ/2}, e^{iθ/2})`.
fn rz(theta: f64) -> M2 {
    [
        [expi(-theta / 2.0), (0.0, 0.0)],
        [(0.0, 0.0), expi(theta / 2.0)],
    ]
}
/// The Hadamard `H = (1/√2)[[1,1],[1,-1]]`.
fn hadamard() -> M2 {
    let s = 1.0 / 2.0f64.sqrt();
    [[(s, 0.0), (s, 0.0)], [(s, 0.0), (-s, 0.0)]]
}
/// The spectral flow's actual projective action: a relative phase `diag(1, e^{iφ})`.
/// (Equal to `R_z(φ)` only up to the global phase `e^{iφ/2}`, so realized words must be
/// compared to targets projectively.)
fn flow_op(phi: f64) -> M2 {
    [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), expi(phi)]]
}
/// An upper bound on the projective operator distance
/// `min_γ max_{ij} |A_{ij} - e^{iγ}B_{ij}|`, evaluated at the Frobenius-optimal global phase
/// `γ = arg(⟨B, A⟩_F)` (the exact minimizer of the Frobenius norm, a valid and conservative
/// gauge for the max-entry norm). Being an upper bound, `abs_error ≤ ε` implies the true
/// projective distance is `≤ ε`.
fn proj_dist(a: &M2, b: &M2) -> f64 {
    let mut ip = (0.0f64, 0.0f64);
    for i in 0..2 {
        for j in 0..2 {
            // conj(b) * a
            let cb = (b[i][j].0, -b[i][j].1);
            ip = cadd(ip, cmul(cb, a[i][j]));
        }
    }
    let gamma = ip.1.atan2(ip.0);
    let phase = expi(gamma);
    let mut bp = [[(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            bp[i][j] = cmul(phase, b[i][j]);
        }
    }
    m2_dist(a, &bp)
}
/// Operator distance `max_{ij} |A_{ij} - B_{ij}|` (a norm-equivalent, deterministic residual).
fn m2_dist(a: &M2, b: &M2) -> f64 {
    let mut d = 0.0f64;
    for i in 0..2 {
        for j in 0..2 {
            let dr = a[i][j].0 - b[i][j].0;
            let di = a[i][j].1 - b[i][j].1;
            d = d.max((dr * dr + di * di).sqrt());
        }
    }
    d
}

impl SkWeaver {
    /// Whether the weaver targets a density-certified carrier (arbitrary-axis synthesis
    /// available) rather than the Clifford carrier.
    #[must_use]
    pub fn is_certified(&self) -> bool {
        self.carrier != Carrier::Clifford
    }

    /// The least flow power `k ≤ bound` whose phase `k·flow_gap` is within `eps` of `phi`
    /// (deterministic bounded search); `None` if the bound is exhausted.
    fn best_flow_power(&self, phi: f64, eps: f64) -> Option<usize> {
        let two_pi = 2.0 * core::f64::consts::PI;
        let target = phi.rem_euclid(two_pi);
        let d = self.flow_gap as f64;
        let mut best_k = 0i64;
        let mut best = Self::circle_dist(0.0, target);
        for k in 0..=self.flow_search_bound {
            let dist = Self::circle_dist((k as f64 * d).rem_euclid(two_pi), target);
            if dist < best {
                best = dist;
                best_k = k;
                if best < eps {
                    break;
                }
            }
        }
        (best < eps).then_some(best_k as usize)
    }

    /// Synthesize an **arbitrary-axis** single-qubit rotation `R_axis(theta)` on a
    /// density-certified carrier via `ZXZ` Euler decomposition, returning the braid word
    /// (`R_x` uses the Clifford Hadamard sandwich `H·Φ^k·H`) and its residual certificate.
    ///
    /// Each Euler `Z`-factor is synthesized to `epsilon/3` as a heuristic budget; the
    /// composite residual is then measured directly---the realized `SU(2)` is formed and its
    /// projective distance to the target computed in `f64`---and the call errors if that
    /// measured residual exceeds `epsilon`. Nothing relies on a norm sub-additivity bound.
    ///
    /// # Errors
    /// On the Clifford carrier (no dense synthesis), or if any Euler factor misses `epsilon/3`.
    pub fn synthesize_axis_rotation(
        &self,
        axis: char,
        theta: f64,
        epsilon: f64,
    ) -> Result<(Vec<BraidGen>, SynthesisResidual), String> {
        if self.carrier == Carrier::Clifford {
            return Err("arbitrary-axis synthesis requires a density-certified carrier".into());
        }
        // Target SU(2) for the requested axis.
        let target = match axis {
            'z' => rz(theta),
            'x' => {
                let h = hadamard();
                m2mul(&m2mul(&h, &rz(theta)), &h)
            }
            'y' => {
                // R_y(θ) = S H R_z(θ) H S†, S = R_z(π/2) (up to global phase). Equivalently
                // build the ZXZ target directly; we form R_y from the standard identity.
                let h = hadamard();
                let s = rz(core::f64::consts::PI / 2.0);
                let sdg = rz(-core::f64::consts::PI / 2.0);
                m2mul(&m2mul(&m2mul(&m2mul(&s, &h), &rz(theta)), &h), &sdg)
            }
            other => return Err(format!("unknown rotation axis `{other}`")),
        };

        // ZXZ Euler angles of `target`: U = R_z(a) R_x(b) R_z(c) with
        // |u00| = cos(b/2), arg(u00) = -(a+c)/2, arg(u01) = -π/2 - (a-c)/2.
        let (u00, u01) = (target[0][0], target[0][1]);
        let mag00 = (u00.0 * u00.0 + u00.1 * u00.1).sqrt().clamp(0.0, 1.0);
        let mag01 = (u01.0 * u01.0 + u01.1 * u01.1).sqrt();
        let b = 2.0 * mag00.acos();
        // Degeneracy is decided by the off-diagonal magnitude (robust to |u00| ≈ 1 rounding):
        // when u01 ≈ 0 the target is a pure Z-rotation and a−c is free (set c = 0); when
        // u00 ≈ 0 (b ≈ π) it is a pure X-rotation and a+c is free (set c = 0).
        let (a, c) = if mag01 < 1e-9 {
            (-2.0 * u00.1.atan2(u00.0), 0.0)
        } else if mag00 < 1e-9 {
            let amc = -2.0 * u01.1.atan2(u01.0) - core::f64::consts::PI;
            (amc, 0.0)
        } else {
            let apc = -2.0 * u00.1.atan2(u00.0);
            let amc = -2.0 * u01.1.atan2(u01.0) - core::f64::consts::PI;
            ((apc + amc) / 2.0, (apc - amc) / 2.0)
        };

        // Per-factor synthesis budget (heuristic); the composite residual is measured and
        // gated directly below, so correctness does not depend on this split.
        let eps3 = epsilon / 3.0;
        let ka = self
            .best_flow_power(a, eps3)
            .ok_or_else(|| format!("Euler R_z(a={a}) missed epsilon/3"))?;
        let kb = self
            .best_flow_power(b, eps3)
            .ok_or_else(|| format!("Euler R_x(b={b}) missed epsilon/3"))?;
        let kc = self
            .best_flow_power(c, eps3)
            .ok_or_else(|| format!("Euler R_z(c={c}) missed epsilon/3"))?;

        // Word: Φ^{ka} · (H Φ^{kb} H) · Φ^{kc}, with H = the Clifford template στσ.
        let had = [BraidGen::Sigma, BraidGen::Tau, BraidGen::Sigma];
        let mut word = Vec::new();
        word.extend(std::iter::repeat_n(BraidGen::Flow, ka));
        word.extend(had.iter().cloned());
        word.extend(std::iter::repeat_n(BraidGen::Flow, kb));
        word.extend(had.iter().cloned());
        word.extend(std::iter::repeat_n(BraidGen::Flow, kc));

        // Realized projective action = flow(ka·d) H flow(kb·d) H flow(kc·d), compared to
        // the target up to a global phase (the carrier acts projectively).
        let realized = m2mul(
            &m2mul(
                &m2mul(
                    &m2mul(&flow_op(self.flow_phase(ka)), &hadamard()),
                    &flow_op(self.flow_phase(kb)),
                ),
                &hadamard(),
            ),
            &flow_op(self.flow_phase(kc)),
        );
        let abs_error = proj_dist(&realized, &target);
        if abs_error > epsilon {
            return Err(format!(
                "arbitrary-axis synthesis residual {abs_error} exceeds epsilon {epsilon}"
            ));
        }
        Ok((
            word,
            SynthesisResidual {
                axis,
                target_angle: theta,
                euler_powers: [ka, kb, kc],
                abs_error,
                epsilon,
            },
        ))
    }

    /// Recompute the realized `SU(2)` operator distance for a residual's Euler powers, so a
    /// shipped residual can be independently verified against the word it travels with.
    #[must_use]
    pub fn verify_residual(&self, r: &SynthesisResidual) -> f64 {
        let target = match r.axis {
            'x' => {
                let h = hadamard();
                m2mul(&m2mul(&h, &rz(r.target_angle)), &h)
            }
            'y' => {
                let h = hadamard();
                let s = rz(core::f64::consts::PI / 2.0);
                let sdg = rz(-core::f64::consts::PI / 2.0);
                m2mul(
                    &m2mul(&m2mul(&m2mul(&s, &h), &rz(r.target_angle)), &h),
                    &sdg,
                )
            }
            _ => rz(r.target_angle),
        };
        let [ka, kb, kc] = r.euler_powers;
        let realized = m2mul(
            &m2mul(
                &m2mul(
                    &m2mul(&flow_op(self.flow_phase(ka)), &hadamard()),
                    &flow_op(self.flow_phase(kb)),
                ),
                &hadamard(),
            ),
            &flow_op(self.flow_phase(kc)),
        );
        proj_dist(&realized, &target)
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
    fn certified_synthesizes_arbitrary_axis_within_epsilon() {
        let p = UseCaseParams::new(4, 3, 8);
        let w = SkWeaver::for_carrier(&p, Carrier::Certified22);
        for &axis in &['x', 'y', 'z'] {
            for &theta in &[0.3f64, 1.0, 2.5, core::f64::consts::PI / 4.0] {
                let (word, res) = w.synthesize_axis_rotation(axis, theta, 0.1).unwrap();
                assert!(!word.is_empty());
                assert!(
                    res.abs_error <= 0.1,
                    "residual {} over epsilon",
                    res.abs_error
                );
                // The shipped residual recomputes exactly (deterministic).
                assert!((w.verify_residual(&res) - res.abs_error).abs() < 1e-12);
            }
        }
        // The Clifford carrier has no arbitrary-axis synthesis.
        let cw = SkWeaver::new(&p);
        assert!(cw.synthesize_axis_rotation('x', 1.0, 0.1).is_err());
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
