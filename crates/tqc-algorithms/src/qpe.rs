//! Quantum Phase Estimation readout, computed exactly.
//!
//! An exact-arithmetic *reference evaluation* of the QPE readout: for a rational phase
//! `θ = a/b` and `M = 2^m` counting states, the measurement distribution of textbook QPE
//! is the Fejér-type kernel `P(k) = |Σ_j e^{2πij(θ − k/M)}|²/M²`, whose peak over the grid
//! is at the `k` minimizing `|θ − k/M|` (a property of the kernel; see e.g. Nielsen–Chuang
//! §5.2). The peak is located here by **exact integer minimization** of `|a·M − k·b|` —
//! no floats, no rounding of rationals through `f64`. This computes the analytic readout
//! of the algorithm; it does not simulate a circuit and claims no speedup.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive, Zero};

/// A QPE readout solver over `2^m` counting states.
pub struct QpeSolver {
    /// The number of counting qubits (precision of the phase estimation).
    pub counting_qubits: usize,
    /// The number of state qubits modeling the eigenvector register. The analytic readout
    /// does not depend on it; it is retained to describe the modeled instance size.
    pub state_qubits: usize,
}

/// Exact readout report for QPE.
#[derive(Debug)]
pub struct ExactQpeReport {
    /// The number of precision bits.
    pub counting_qubits: usize,
    /// The exact underlying phase eigenvalue.
    pub exact_phase: BigRational,
    /// The estimated phase `k/M` at the exact kernel peak.
    pub estimated_phase: BigRational,
    /// The peak measurement integer `k`.
    pub measured_integer: usize,
    /// The exact estimation error `|θ − k/M|` on the circle, guaranteed `≤ 1/(2M)`.
    pub estimation_error: BigRational,
}

impl QpeSolver {
    /// Initializes the QPE solver.
    #[must_use]
    pub fn new(counting_qubits: usize, state_qubits: usize) -> Self {
        Self {
            counting_qubits,
            state_qubits,
        }
    }

    /// Locates the exact peak of the QPE measurement kernel for the given rational phase
    /// by integer minimization of `|a·M − k·b|` over `k ∈ [0, M)` (phase taken mod 1, and
    /// distance taken on the circle), then verifies the QPE guarantee `|θ − k/M| ≤ 1/(2M)`
    /// exactly.
    ///
    /// # Errors
    /// If the phase normalization or the exact guarantee fails (an internal defect).
    pub fn execute_exact_witness(
        &self,
        true_phase: &BigRational,
    ) -> Result<ExactQpeReport, String> {
        let m_states = BigInt::from(1u8) << self.counting_qubits;

        // Normalize θ into [0, 1): phases are defined mod 1.
        let theta = true_phase - true_phase.floor();
        if theta < BigRational::zero() || theta >= BigRational::from_integer(BigInt::from(1)) {
            return Err("phase normalization failed".into());
        }

        // Exact integer minimization of the circle distance |a·M − k·b| over k ∈ [0, M):
        // the minimizer is floor(θM) or its successor mod M.
        let a = theta.numer().clone();
        let b = theta.denom().clone();
        let am = &a * &m_states;
        let k_floor = &am / &b; // exact floor: a, b, M are non-negative
        let k_next = (&k_floor + 1) % &m_states;
        let err_at = |k: &BigInt| -> BigInt {
            let d = (&am - k * &b).abs();
            let wrap = (&b * &m_states) - &d;
            d.min(wrap)
        };
        let e_floor = err_at(&k_floor);
        let e_next = err_at(&k_next);
        let (k_best, err_best) = if e_next < e_floor {
            (k_next, e_next)
        } else {
            (k_floor, e_floor)
        };

        // The QPE guarantee, verified exactly: |θ − k/M| ≤ 1/(2M)  ⇔  2·|aM − kb| ≤ b.
        if BigInt::from(2) * &err_best > b {
            return Err(format!(
                "exact QPE guarantee violated: 2·|aM − kb| = {} > b = {}",
                BigInt::from(2) * &err_best,
                b
            ));
        }

        let measured_integer = k_best
            .to_usize()
            .ok_or_else(|| "peak index does not fit in usize".to_owned())?;
        let estimated_phase = BigRational::new(k_best, m_states.clone());
        let estimation_error = BigRational::new(err_best, b * m_states);

        Ok(ExactQpeReport {
            counting_qubits: self.counting_qubits,
            exact_phase: true_phase.clone(),
            estimated_phase,
            measured_integer,
            estimation_error,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qpe_exact_witness() {
        let solver = QpeSolver::new(3, 1);
        let true_phase = BigRational::new(BigInt::from(1), BigInt::from(8));
        let report = solver.execute_exact_witness(&true_phase).unwrap();
        assert_eq!(report.measured_integer, 1);
        assert_eq!(report.estimated_phase, true_phase);
        assert!(report.estimation_error.is_zero());
    }

    #[test]
    fn test_qpe_non_representable_phase_meets_guarantee() {
        // θ = 1/3 with 3 counting bits: not representable on the grid; the peak must be
        // the nearest grid point with exact circle error ≤ 1/16.
        let solver = QpeSolver::new(3, 1);
        let true_phase = BigRational::new(BigInt::from(1), BigInt::from(3));
        let report = solver.execute_exact_witness(&true_phase).unwrap();
        assert_eq!(report.measured_integer, 3); // 3/8 is the nearest grid point to 1/3
        let bound = BigRational::new(BigInt::from(1), BigInt::from(16));
        assert!(report.estimation_error <= bound);
    }

    #[test]
    fn test_qpe_negative_phase_normalizes() {
        // θ = −1/8 ≡ 7/8 (mod 1): must not silently map to 0.
        let solver = QpeSolver::new(3, 1);
        let true_phase = BigRational::new(BigInt::from(-1), BigInt::from(8));
        let report = solver.execute_exact_witness(&true_phase).unwrap();
        assert_eq!(report.measured_integer, 7);
    }

    #[test]
    fn test_qpe_wraparound_phase_peaks_at_zero() {
        // θ = 31/32 with 3 counting bits: the circle-nearest grid point is 0/8 (distance
        // 1/32 across the wrap), not 7/8 (distance 3/32).
        let solver = QpeSolver::new(3, 1);
        let true_phase = BigRational::new(BigInt::from(31), BigInt::from(32));
        let report = solver.execute_exact_witness(&true_phase).unwrap();
        assert_eq!(report.measured_integer, 0);
        assert_eq!(
            report.estimation_error,
            BigRational::new(BigInt::from(1), BigInt::from(32))
        );
    }
}
