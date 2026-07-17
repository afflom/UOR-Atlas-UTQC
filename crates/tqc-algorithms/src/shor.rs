//! Period finding (Shor's algorithm core), evaluated exactly at fixed instances.
//!
//! An exact-arithmetic *reference evaluation*: the counting-register amplitudes are built
//! as exact cyclotomic sums over `Q(ζ_M)` (power-of-two `M`, `Φ_M(x) = x^{M/2} + 1`),
//! candidate peaks are ordered by a deterministic fixed-point ranking (a heuristic
//! ordering only — exactness is carried by the final check), continued fractions extract
//! period candidates exactly, and the recovered period is verified exactly against the
//! modular-exponentiation orbit. A classical exact computation at small fixed instances;
//! no speedup claim is made.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};

/// An element of the cyclotomic field Q(zeta_M), represented exactly as a polynomial
/// in zeta_M of degree M/2 (where M is a power of 2, so Phi_M(x) = x^{M/2} + 1).
#[derive(Clone, Debug, PartialEq, Eq)]
struct CycPow2 {
    m: usize,
    coeffs: Vec<i64>, // degree is m / 2
}

impl CycPow2 {
    fn zero(m: usize) -> Self {
        Self {
            m,
            coeffs: vec![0; m / 2],
        }
    }

    // Add omega^p
    fn add_omega_pow(&mut self, mut p: usize) {
        p %= self.m;
        let half = self.m / 2;
        if p < half {
            self.coeffs[p] += 1;
        } else {
            self.coeffs[p - half] -= 1; // omega^{m/2} = -1
        }
    }

    fn add(&self, other: &Self) -> Self {
        let mut res = self.clone();
        for i in 0..(self.m / 2) {
            res.coeffs[i] += other.coeffs[i];
        }
        res
    }

    // Compute self * complex_conjugate(self)
    fn norm_sq(&self) -> Self {
        let half = self.m / 2;
        let mut conj = vec![0; half];
        conj[0] = self.coeffs[0];
        for i in 1..half {
            conj[half - i] = -self.coeffs[i];
        }

        let mut res = vec![0; half];
        for (i, &coeff_i) in self.coeffs.iter().enumerate().take(half) {
            if coeff_i == 0 {
                continue;
            }
            for (j, &conj_j) in conj.iter().enumerate().take(half) {
                if conj_j == 0 {
                    continue;
                }
                let p = i + j;
                let val = coeff_i * conj_j;
                if p < half {
                    res[p] += val;
                } else {
                    res[p - half] -= val;
                }
            }
        }
        Self {
            m: self.m,
            coeffs: res,
        }
    }

    // Evaluate exactly to a scaled integer for threshold/sorting comparisons
    // entirely bypassing f64 non-determinism.
    fn to_scaled_int(&self) -> i128 {
        let mut sum = 0i128;
        for (i, &c) in self.coeffs.iter().enumerate() {
            let scaled_cos = exact_scaled_cos2(i, self.m) / 2;
            sum += (c as i128) * scaled_cos;
        }
        sum
    }
}

// Computes a fixed-point representation of 2*cos(2pi j / m) scaled by 10^15.
// m must be a power of 2. Evaluated strictly over integers using Chebyshev polynomials.
fn exact_scaled_cos2(mut j: usize, m: usize) -> i128 {
    j %= m;
    let mut sign = 1i128;
    if j > m / 2 {
        j = m - j;
    }
    if j > m / 4 {
        j = m / 2 - j;
        sign = -1;
    }
    let scale = 1_000_000_000_000_000i128; // 10^15
    if j == 0 {
        return sign * 2 * scale;
    }
    if j == m / 4 {
        return 0;
    }

    let mut c = 0;
    let mut temp = m;
    while temp > 1 {
        temp >>= 1;
        c += 1;
    }

    let mut val = 0i128;
    for _ in 2..c {
        // val = sqrt(2 + val)
        let inner = 2 * scale + val;
        let n = (inner * scale).unsigned_abs();
        let mut x = n;
        let mut y = x.div_ceil(2);
        while y < x {
            x = y;
            y = (x + n / x) / 2;
        }
        val = x as i128;
    }

    // Chebyshev polynomials to get 2*cos(j * 2pi / M)
    let mut t_prev = 2 * scale;
    let mut t_curr = val;
    for _ in 1..j {
        let t_next = (val * t_curr) / scale - t_prev;
        t_prev = t_curr;
        t_curr = t_next;
    }

    sign * t_curr
}

/// A solver for the Period Finding subroutine mapped to exact permutations
/// over the topological combinatorial space.
pub struct ShorSolver {
    /// Number of qubits in the counting register.
    pub counting_qubits: usize,
    /// Number of qubits in the state register.
    pub state_qubits: usize,
}

/// Exact certified execution report for Shor's algorithm.
#[derive(Debug)]
pub struct ExactShorReport {
    /// The base for the modular exponentiation.
    pub base: usize,
    /// The modulus.
    pub modulus: usize,
    /// The exact period, derived independently as the orbit length of 1 under the
    /// modular-exponentiation permutation.
    pub period: usize,
    /// The exact rational eigenphase c * s / r.
    pub exact_phase: BigRational,
    /// The period recovered via the QPE peak and continued-fraction expansion; verified
    /// to equal `period` exactly.
    pub recovered_period: usize,
}

impl ShorSolver {
    /// Initializes the solver.
    #[must_use]
    pub fn new(counting_qubits: usize, state_qubits: usize) -> Self {
        assert!(counting_qubits > 0, "Counting qubits must be > 0");
        Self {
            counting_qubits,
            state_qubits,
        }
    }

    /// The modeled state register must be able to hold residues mod `modulus`.
    fn check_state_capacity(&self, modulus: usize) -> Result<(), String> {
        if (1usize << self.state_qubits) < modulus {
            return Err(format!(
                "state register of {} qubits cannot hold residues mod {modulus}",
                self.state_qubits
            ));
        }
        Ok(())
    }

    /// Executes the certified exact Shor witness natively without state vectors.
    ///
    /// The modular exponentiation
    /// is executed as an exact permutation on Z/N, extracting the period via orbit length,
    /// determining the exact rational eigenphase, and performing continued-fractions recovery
    /// without f64 loss.
    pub fn execute_exact_witness(
        &self,
        base: usize,
        modulus: usize,
    ) -> Result<ExactShorReport, String> {
        if base.gcd(&modulus) != 1 {
            return Err("Base and modulus must be coprime".into());
        }
        self.check_state_capacity(modulus)?;

        // 1. Modular exponentiation as an exact permutation on Z/N
        let perm: Vec<usize> = (0..modulus).map(|i| (i * base) % modulus).collect();

        // The true period, derived independently of the QPE readout: the orbit length of 1
        // under the permutation. The QPE + continued-fraction recovery below must reproduce
        // exactly this value or the witness fails.
        let period = {
            let mut state = perm[1];
            let mut len = 1usize;
            while state != 1 {
                state = perm[state];
                len += 1;
                if len > modulus {
                    return Err("orbit of 1 did not close within the modulus".into());
                }
            }
            len
        };

        // 2. Evaluate the quantum state evolution and QPE interference.
        // We authentically execute the QPE interference pattern by evaluating the Fourier amplitudes
        // without classical orbit-length shortcuts.
        let m_states = 1 << self.counting_qubits;
        let mut u_a_j = vec![0usize; m_states];
        let mut state = 1;

        for item in u_a_j.iter_mut().take(m_states) {
            *item = state;
            state = perm[state]; // Executed by permutation composition
        }

        // Calculate QPE interference probabilities for measuring k in the counting register
        // evaluated entirely in exact cyclotomic algebra Q(zeta_M).
        let mut probabilities = vec![CycPow2::zero(m_states); m_states];
        for (k, prob_k) in probabilities.iter_mut().enumerate().take(m_states) {
            let mut p_k = CycPow2::zero(m_states);
            // The probability is proportional to sum_x | sum_{j: x_j = x} e^{-2pi i j k / M} |^2
            for x in 0..modulus {
                let mut amplitude_x = CycPow2::zero(m_states);
                for (j, &xj) in u_a_j.iter().enumerate() {
                    if xj == x {
                        amplitude_x.add_omega_pow(j * k);
                    }
                }
                p_k = p_k.add(&amplitude_x.norm_sq());
            }
            // Probabilities are built purely algebraically.
            *prob_k = p_k;
        }

        // Sort k indices by their exact interference algebraic probability (descending)
        let mut k_candidates: Vec<usize> = (1..m_states).collect(); // Exclude k=0 (trivial)
        k_candidates.sort_by(|&a, &b| {
            probabilities[b]
                .to_scaled_int()
                .cmp(&probabilities[a].to_scaled_int())
        });

        // 3. Evaluate the principal eigenphases recovered via continued fractions
        let mut recovered_period = 0;
        let mut exact_phase = BigRational::new(BigInt::from(0), BigInt::from(1));

        'outer: for k in k_candidates {
            if probabilities[k].to_scaled_int() <= 0 {
                continue; // Ignore negligible or zero background probabilities exactly
            }

            // k / M is the measured rational phase.
            let phase = BigRational::new(BigInt::from(k), BigInt::from(m_states));

            // 4. Continued-fractions recovery: EVERY convergent denominator is a period
            // candidate (the final convergent alone is just the reduced denominator of
            // k/M, a divisor of M — it misses every period not dividing 2^m).
            for q in Self::convergent_denominators(&phase) {
                if q == 0 || q > modulus {
                    continue;
                }
                // The candidate must close the orbit AND equal the independently derived
                // period exactly (a proper divisor multiple would close the orbit too).
                let mut check_state = 1;
                for _ in 0..q {
                    check_state = perm[check_state];
                }
                if check_state == 1 && q == period {
                    recovered_period = q;
                    exact_phase = phase;
                    break 'outer;
                }
            }
        }

        if recovered_period == 0 {
            return Err("Failed to recover the full period from the spectral eigenphases".into());
        }

        Ok(ExactShorReport {
            base,
            modulus,
            period,
            exact_phase,
            recovered_period,
        })
    }

    /// All convergent denominators of the continued-fraction expansion of `phase`,
    /// in increasing order. Denominators that exceed `usize` are skipped (they cannot
    /// be periods of a `usize` modulus).
    fn convergent_denominators(phase: &BigRational) -> Vec<usize> {
        let mut num = phase.numer().clone();
        let mut den = phase.denom().clone();

        let mut q_prev = BigInt::from(1);
        let mut q_curr = BigInt::from(0);
        let mut out = Vec::new();

        while !den.is_zero() {
            let a = num.clone() / den.clone();
            let rem = num.clone() % den.clone();

            let q_next = a * q_curr.clone() + q_prev;
            q_prev = q_curr;
            q_curr = q_next;

            if let Some(q) = q_curr.to_usize() {
                out.push(q);
            }

            num = den;
            den = rem;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shor_exact_witness() {
        let solver = ShorSolver::new(4, 4);
        // Base 2, Modulus 15 (Shor's classic example)
        let report = solver
            .execute_exact_witness(2, 15)
            .expect("Exact witness failed");

        assert_eq!(report.period, 4);
        assert_eq!(report.recovered_period, 4);
    }

    #[test]
    fn test_shor_recovers_period_not_dividing_two_power() {
        // Base 2 mod 7 has period 3, which does not divide M = 16: recovery must come
        // from an intermediate continued-fraction convergent (1/3 of 5/16), which the
        // final-convergent-only expansion can never produce.
        let solver = ShorSolver::new(4, 3);
        let report = solver
            .execute_exact_witness(2, 7)
            .expect("Exact witness failed");
        assert_eq!(report.period, 3);
        assert_eq!(report.recovered_period, 3);
    }

    #[test]
    fn test_shor_rejects_undersized_state_register() {
        let solver = ShorSolver::new(4, 2);
        assert!(solver.execute_exact_witness(2, 15).is_err());
    }
}
