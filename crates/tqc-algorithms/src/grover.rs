//! Grover amplitude evolution, evaluated exactly.
//!
//! An exact-arithmetic *reference evaluation* of Grover's amplitude recurrence at a fixed
//! instance: the two amplitude coefficients (target / non-target) evolve under the exact
//! rational reflection recurrence, with the optimal iteration count from the exact integer
//! square root. This validates the algorithm mathematics; it is a classical exact
//! computation and claims no speedup.

use num_bigint::BigInt;
use num_rational::BigRational;

/// A Grover's Search solver mapped to the topological space.
pub struct GroverSolver {
    /// The number of virtual qubits in the search space.
    pub num_qubits: usize,
}

/// Exact certified execution report for Grover's algorithm.
#[derive(Debug)]
pub struct ExactGroverReport {
    /// Number of virtual qubits.
    pub num_qubits: usize,
    /// Target state to search for.
    pub target_state: usize,
    /// Exact optimal number of iterations executed.
    pub iterations: usize,
    /// Target state amplitude coefficient (multiplied by 1/sqrt(N)).
    pub target_amplitude_coeff: BigRational,
    /// Non-target state amplitude coefficient (multiplied by 1/sqrt(N)).
    pub non_target_amplitude_coeff: BigRational,
}

impl GroverSolver {
    /// Initializes the solver for a given number of virtual qubits.
    #[must_use]
    pub fn new(num_qubits: usize) -> Self {
        Self { num_qubits }
    }

    /// Executes the certified exact Grover witness natively.
    ///
    /// Evaluates the sequence of Oracle and Diffuser reflections exactly on the
    /// amplitude basis, at the fixed reference instance; no #P-hardness claim is attached.
    pub fn execute_exact_witness(&self, target_state: usize) -> Result<ExactGroverReport, String> {
        let n_states = 1 << self.num_qubits;
        if target_state >= n_states {
            return Err("Target state out of bounds".into());
        }

        // Calculate optimal iterations: floor(pi/4 * sqrt(N))
        // using exact scaled integer arithmetic without f64
        let n_states_u128 = n_states as u128;
        // pi^2 / 16 scaled by 10^16
        let pi_sq_over_16_scaled = 6168502750680849u128;
        let iters_scaled_sq = n_states_u128 * pi_sq_over_16_scaled;

        let mut x = iters_scaled_sq;
        let mut y = x.div_ceil(2);
        while y < x {
            x = y;
            y = (x + iters_scaled_sq / x) / 2;
        }
        let iterations = (x / 100_000_000) as usize;

        // Exact amplitude tracking on the two-coefficient symmetric subspace (the
        // Grover state stays in span{|target>, uniform-rest} throughout)
        let mut a_t = BigRational::new(BigInt::from(1), BigInt::from(1));
        let mut a_u = BigRational::new(BigInt::from(1), BigInt::from(1));
        let n_states_bi = BigRational::new(BigInt::from(n_states), BigInt::from(1));
        let n_minus_1_bi = BigRational::new(BigInt::from(n_states - 1), BigInt::from(1));
        let two = BigRational::new(BigInt::from(2), BigInt::from(1));

        for _ in 0..iterations {
            // 1. Oracle: Invert phase of the target state
            a_t = -a_t;

            // 2. Diffuser: Inversion about the mean
            let mean = (&a_t + &n_minus_1_bi * &a_u) / &n_states_bi;
            a_t = &two * &mean - a_t;
            a_u = &two * &mean - a_u;
        }

        Ok(ExactGroverReport {
            num_qubits: self.num_qubits,
            target_state,
            iterations,
            target_amplitude_coeff: a_t,
            non_target_amplitude_coeff: a_u,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grover_exact_witness() {
        let solver = GroverSolver::new(3);
        let report = solver.execute_exact_witness(5).unwrap();
        // For 3 qubits, 2 iterations is optimal.
        // Probability = (11/4)^2 / 8 = 121 / 128
        let p_t = &report.target_amplitude_coeff * &report.target_amplitude_coeff
            / BigRational::new(BigInt::from(8), BigInt::from(1));
        let expected_p_t = BigRational::new(BigInt::from(121), BigInt::from(128));
        assert_eq!(p_t, expected_p_t);
    }
}
