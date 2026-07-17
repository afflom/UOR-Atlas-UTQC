//! Callable check bodies for the algorithm dictionary rows.
//!
//! One implementation per row (DRY): both the cucumber steps in `tqc-conformance` and the
//! `xtask report` ledger call these. Each is an exact-arithmetic reference evaluation at a
//! fixed instance — see the module docs of the corresponding solver for the precise scope.

use crate::{grover::GroverSolver, qft::QftSolver, qpe::QpeSolver, shor::ShorSolver};
use num_bigint::BigInt;
use num_rational::BigRational;
use tqc_compiler::{BraidGen, Compiler};
use tqc_core::generators::{Generators, Permutation};
use tqc_core::UseCaseParams;

/// Grover at the fixed instance n = 3 qubits, target 5: the exact rational amplitude
/// recurrence must yield success probability exactly 121/128 after the optimal 2 iterations.
///
/// # Errors
/// If the exact amplitude evolution deviates from the closed-form value.
pub fn grover_check() -> Result<(), String> {
    let solver = GroverSolver::new(3);
    let report = solver.execute_exact_witness(5)?;
    let p_t = &report.target_amplitude_coeff * &report.target_amplitude_coeff
        / BigRational::new(BigInt::from(8), BigInt::from(1));
    let expected = BigRational::new(BigInt::from(121), BigInt::from(128));
    if p_t != expected {
        return Err(format!(
            "Grover exact success probability {p_t} != {expected}"
        ));
    }
    Ok(())
}

/// QPE at the fixed instance m = 3 counting bits, θ = 1/8: the exact integer peak must be
/// k = 1 with zero estimation error; and at the non-representable θ = 1/3 the exact
/// guarantee |θ − k/M| ≤ 1/(2M) must hold.
///
/// # Errors
/// If the exact readout deviates.
pub fn qpe_check() -> Result<(), String> {
    let solver = QpeSolver::new(3, 1);
    let exact = BigRational::new(BigInt::from(1), BigInt::from(8));
    let report = solver.execute_exact_witness(&exact)?;
    if report.measured_integer != 1 || report.estimated_phase != exact {
        return Err("QPE exact peak at θ=1/8 deviates".into());
    }
    let third = BigRational::new(BigInt::from(1), BigInt::from(3));
    let report = solver.execute_exact_witness(&third)?;
    if report.measured_integer != 3 {
        return Err("QPE peak at θ=1/3 is not the nearest grid point 3/8".into());
    }
    Ok(())
}

/// Shor period finding at the fixed instance a = 2, N = 15: the exact cyclotomic
/// amplitude construction plus continued fractions must recover period 4, verified
/// exactly against the modular-exponentiation orbit.
///
/// # Errors
/// If the recovered period deviates.
pub fn shor_check() -> Result<(), String> {
    let solver = ShorSolver::new(4, 4);
    let report = solver.execute_exact_witness(2, 15)?;
    if report.period != 4 || report.recovered_period != 4 {
        return Err(format!(
            "Shor period {} / recovered {} != 4",
            report.period, report.recovered_period
        ));
    }
    Ok(())
}

/// QFT scheduling at the fixed instance of 4 virtual qubits: the circuit must schedule
/// onto a braid word (rotations admitted only as exact discrete phases, ε = 0.5 — at this
/// tolerance the small controlled-phase half-rotations resolve to the exact identity word,
/// so the scheduled content is the H/CNOT template words), the word must stay within the
/// documented bound, and executing it as a permutation of the class space must be total.
/// Returns the final integer state for the caller to content-address.
///
/// # Errors
/// If scheduling fails or the word exceeds its bound.
pub fn qft_word_check(p: &UseCaseParams) -> Result<Vec<i64>, String> {
    let solver = QftSolver::new(4);
    let circuit = solver.build_circuit();
    let compiler = Compiler::new(p);
    let word = compiler.compile(&circuit, 0.5)?;
    if word.sequence.is_empty() {
        return Err("QFT scheduled to an empty braid word".into());
    }
    if word.sequence.len() >= 2000 {
        return Err(format!(
            "QFT braid word length {} exceeds the documented bound",
            word.sequence.len()
        ));
    }
    let g = Generators::new(p);
    let mut perm = Permutation::identity(p.class_count());
    for op in &word.sequence {
        let p_op = match op {
            BraidGen::Sigma => &g.sigma,
            BraidGen::Tau => &g.tau,
            BraidGen::Mu => &g.mu,
        };
        perm = perm.then(p_op);
    }
    let n = p.class_count() as i64;
    let base: Vec<i64> = (0..n).collect();
    Ok(perm.permute_amplitudes(&base))
}
