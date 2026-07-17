//! The Atlas-native modular tensor category construction.
//!
//! This module constructs the pointed modular category `C(A, q)` of the abelian group
//! `A = Z_modality × Z_2^3` with the quadratic form `q(m, c) = ζ_modality^{m²} · i^{|c|}`
//! (for the Atlas instance: `A = Z_3 × Z_2^3`, 24 simple objects — three semion factors
//! times a `Z_3` anyon).
//!
//! # Construction (build-level, axiom-validated)
//!
//! - **Class-space quotient.** The Atlas class count is a `Z_scope`-graded extension of the
//!   carrier dimension; `construct_atlas_native` checks `class_count / scope = carrier_dim`
//!   and rejects any parameter tuple where the quotient does not close.
//! - **Fusion ring.** Group-law fusion on `A` — the non-negative associative quotient of
//!   the signed octonion composition (associativity of the absolute quotient is proven
//!   separately by `absolute_quotient_is_associative` in `tqc-core`).
//! - **Braiding and associator.** The `R`- and `F`-symbols are the Eilenberg–MacLane
//!   abelian 3-cocycle of the quadratic form: per semion factor
//!   `ω(a,b,c) = (−1)^{abc}`, `R(a,b) = i^{ab}`; the `Z_modality` factor is bilinear
//!   (`R(m,m') = ζ^{mm'}`, trivial cocycle). The full data passes the phase-exact
//!   pentagon, hexagon, balancing, Verlinde, and monodromy–S checks in
//!   [`crate::verifier::verify_mtc_axioms`].
//!
//! The category has nonzero central charge (for the Atlas instance `c ≡ 5 (mod 8)`), so
//! `(ST)³ = p⁺S²` with a non-trivial anomaly `p⁺ = e^{2πi c/8}`; the verifier derives `p⁺`
//! from the Gauss sum rather than assuming `(ST)³ = S²`.
//!
//! This is a `build`-level construction: it is validated against the universal MTC axioms;
//! it is *not* claimed as an F1-sourced fact.

use crate::verifier::ModularData;
use tqc_core::params::UseCaseParams;

/// Represents the failure to construct an Atlas-native MTC.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstructionObstruction {
    /// Dimension mismatch between generated class space and parameter bounds.
    DimensionMismatch(u64, u64),
    /// The context is not a power of two, so the `Z_2^k` composition quotient does not exist.
    ContextNotTwoGroup(u32),
}

impl core::fmt::Display for ConstructionObstruction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DimensionMismatch(a, b) => {
                write!(f, "class space {} does not quotient to carrier {}", a, b)
            }
            Self::ContextNotTwoGroup(o) => {
                write!(f, "context {o} is not a power of two: no Z_2^k quotient")
            }
        }
    }
}

impl std::error::Error for ConstructionObstruction {}

/// Attempt to construct an Atlas-native MTC from parameters.
use crate::{Matrix, C};

/// The true, explicit Atlas-native MTC, constructed from the structural absolute quotient
/// of the g2 composition ring (Context 8, Modality 3).
#[derive(Clone, Debug)]
#[allow(clippy::needless_range_loop)]
pub struct AtlasNative {
    /// The condensed carrier dimension (24).
    pub carrier_dim: usize,
    /// The use-case modality parameter.
    pub modality: usize,
    /// The use-case context parameter.
    pub context: usize,
}

impl ModularData for AtlasNative {
    fn dim(&self) -> usize {
        self.carrier_dim
    }

    #[allow(clippy::needless_range_loop)]
    fn s_matrix(&self) -> Matrix {
        let n = self.carrier_dim;
        let mut s = vec![vec![C::new(0.0, 0.0); n]; n];
        let root24 = (n as f64).sqrt();
        for x in 0..n {
            let m1 = x / self.context;
            let c1 = x % self.context;
            for y in 0..n {
                let m2 = y / self.context;
                let c2 = y % self.context;

                // Modality Z_{modality}: S carries the INVERSE monodromy bicharacter
                // χ(a,b)⁻¹ = e^{−2πi·k·m1m2/n} at the same level k as the twist
                // (k = 2 for odd n, k = 1 for even n). Using a fixed exponent +2·m1m2/n
                // is value-identical only for n ∈ {1,2,3} and breaks the monodromy–S and
                // anomaly relations for every larger modality.
                let k_mod = if self.modality % 2 == 0 { 1.0 } else { 2.0 };
                let theta = -2.0 * core::f64::consts::PI * k_mod * (m1 * m2) as f64
                    / (self.modality as f64);
                let phase3 = C::phase(theta);

                // Context Z_2^k (derived from the associative absolute quotient of octonions):
                // the semion bicharacter is real, hence self-inverse.
                let dot = (c1 & c2).count_ones();
                let phase2 = if dot % 2 == 1 { -1.0 } else { 1.0 };

                s[x][y] = phase3.scale(phase2 / root24);
            }
        }
        s
    }

    #[allow(clippy::needless_range_loop)]
    fn t_diag(&self) -> Vec<C> {
        let n = self.carrier_dim;
        let mut t = vec![C::new(0.0, 0.0); n];
        for x in 0..n {
            let m = x / self.context;
            let c = x % self.context;

            // Modality Z_{modality} pseudo-metric
            let k_mod = if self.modality % 2 == 0 { 1.0 } else { 2.0 };
            let theta = core::f64::consts::PI * k_mod * (m * m) as f64 / (self.modality as f64);
            let phase3 = C::phase(theta);

            // Context Z_2^3 pseudo-metric: q(c) = i^{c_0 + c_1 + c_2}
            let sum = c.count_ones();
            let phase2 = match sum % 4 {
                0 => C::new(1.0, 0.0),
                1 => C::new(0.0, 1.0),
                2 => C::new(-1.0, 0.0),
                3 => C::new(0.0, -1.0),
                _ => unreachable!(),
            };
            t[x] = phase3.times(phase2);
        }
        t
    }

    #[allow(clippy::needless_range_loop)]
    fn charge_conjugation(&self) -> Matrix {
        let n = self.carrier_dim;
        let mut c_mat = vec![vec![C::new(0.0, 0.0); n]; n];
        for x in 0..n {
            let m = x / self.context;
            let c = x % self.context;
            let m_inv = (self.modality - m % self.modality) % self.modality;
            // c is its own inverse in Z_2^k (context)
            let inv = m_inv * self.context + c;
            c_mat[x][inv] = C::new(1.0, 0.0);
        }
        c_mat
    }

    fn n_ijk(&self, i: usize, j: usize, k: usize) -> f64 {
        let m1 = i / self.context;
        let c1 = i % self.context;
        let m2 = j / self.context;
        let c2 = j % self.context;
        let m3 = k / self.context;
        let c3 = k % self.context;

        let m_add = (m1 + m2) % self.modality;
        let c_add = c1 ^ c2;

        if m_add == m3 && c_add == c3 {
            1.0
        } else {
            0.0
        }
    }

    fn f_symbol(&self, i: usize, j: usize, k: usize, l: usize, m: usize, n: usize) -> C {
        let m1 = i / self.context;
        let c1 = i % self.context;
        let m2 = j / self.context;
        let c2 = j % self.context;
        let m3 = k / self.context;
        let c3 = k % self.context;
        let m_m = m / self.context;
        let c_m = m % self.context;
        let m_n = n / self.context;
        let c_n = n % self.context;
        let m_l = l / self.context;
        let c_l = l % self.context;

        if (m1 + m2) % self.modality == m_m
            && c1 ^ c2 == c_m
            && (m_m + m3) % self.modality == m_l
            && c_m ^ c3 == c_l
            && (m2 + m3) % self.modality == m_n
            && c2 ^ c3 == c_n
            && (m1 + m_n) % self.modality == m_l
            && c1 ^ c_n == c_l
        {
            // Z_2^k 3-cocycle associator (per semion factor): F(c1, c2, c3) = (-1)^{\sum (c1_i * c2_i * c3_i)}
            let dot3 = (c1 & c2 & c3).count_ones();
            let mut phase = if dot3 % 2 == 1 { -1.0 } else { 1.0 };
            // Z_modality Eilenberg–MacLane cocycle for q(m) = e^{πi·k·m²/n}: trivial for odd
            // n (k = 2, bilinear form), and ω(m1,m2,m3) = (−1)^{m1·⌊(m2+m3)/n⌋} for even n
            // (k = 1, the level-2n quadratic form is not bilinear).
            if self.modality % 2 == 0 && (m1 * ((m2 + m3) / self.modality)) % 2 == 1 {
                phase = -phase;
            }
            C::new(phase, 0.0)
        } else {
            C::new(0.0, 0.0)
        }
    }

    fn r_symbol(&self, x: usize, y: usize, k: usize) -> C {
        let m1 = x / self.context;
        let c1 = x % self.context;
        let m2 = y / self.context;
        let c2 = y % self.context;
        let m3 = k / self.context;
        let c3 = k % self.context;

        if (m1 + m2) % self.modality == m3 && c1 ^ c2 == c3 {
            // Z_{modality} R-matrix phase e^{πi·k·m1m2/n}, at the same level k as the twist
            // q(m) = e^{πi·k·m²/n} (k = 2 for odd n, k = 1 for even n) — the braiding half
            // of the Eilenberg–MacLane abelian 3-cocycle of q.
            let k_mod = if self.modality % 2 == 0 { 1.0 } else { 2.0 };
            let theta = core::f64::consts::PI * k_mod * (m1 * m2) as f64 / (self.modality as f64);
            let phase3 = C::phase(theta);

            // Z_2^3 R-matrix phase: i^{c1 . c2}
            let dot = (c1 & c2).count_ones();
            let phase2 = match dot % 4 {
                0 => C::new(1.0, 0.0),
                1 => C::new(0.0, 1.0),
                2 => C::new(-1.0, 0.0),
                3 => C::new(0.0, -1.0),
                _ => unreachable!(),
            };

            phase3.times(phase2)
        } else {
            C::new(0.0, 0.0)
        }
    }
}

/// Construct the pointed abelian MTC stand-in from parameters.
pub fn construct_atlas_native(
    p: &UseCaseParams,
) -> Result<Box<dyn ModularData>, ConstructionObstruction> {
    // 1. Z_q Equivariant Gauging: The Atlas class count (96) is a Z_q extension
    // of the base carrier dimension (24). We quotient by the scope parameter q.
    let base_dim = p.class_count() / (p.scope as u64);
    if base_dim != p.carrier_dim() {
        return Err(ConstructionObstruction::DimensionMismatch(
            p.class_count(),
            p.carrier_dim(),
        ));
    }
    if !p.context.is_power_of_two() {
        return Err(ConstructionObstruction::ContextNotTwoGroup(p.context));
    }

    // 2. Structural Absolute Quotient: The non-negative fusion quotient of the
    // signed octonion algebra is strictly associative (proven by `absolute_quotient_is_associative(8)`),
    // mapping the pseudo-unitary algebra cleanly into a unitary MTC.

    // Return the pointed abelian quotient construction.
    Ok(Box::new(AtlasNative {
        carrier_dim: p.carrier_dim() as usize,
        modality: p.modality as usize,
        context: p.context as usize,
    }))
}
