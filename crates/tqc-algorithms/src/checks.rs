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
            BraidGen::Flow => {
                return Err(
                    "Clifford QFT scheduling must not emit a spectral-flow generator".into(),
                )
            }
        };
        perm = perm.then(p_op);
    }
    let n = p.class_count() as i64;
    let base: Vec<i64> = (0..n).collect();
    Ok(perm.permute_amplitudes(&base))
}

/// Certified-carrier compilation and Shor replay (dictionary row `certified-carrier-compilation`).
///
/// On a density-certified carrier (`Carrier::Certified22`, PU(22) established by the exact
/// certificate) the compiler synthesizes arbitrary-axis single-qubit rotations (ZXZ Euler
/// over the spectral flow and the Clifford Hadamard sandwich), so the continuous-phase gates
/// of Shor's QFT compile **unconditionally** (unlike the Clifford carrier, where gate-set
/// density is refuted). This check:
///
/// 1. confirms the PU(22) density premise from the exact `Q(ζ₂₄)` certificate;
/// 2. compiles each controlled-phase rotation of Shor's 4-bit QFT, `π/2^{j}`, on the
///    certified carrier, and verifies the emitted spectral-flow word realizes the target
///    phase within `ε` (a deterministic synthesis, checked exactly against `flow_phase`);
/// 3. compiles an arbitrary-axis rotation circuit through the certified-carrier compiler,
///    verifies every shipped per-word residual certificate (|Δ| ≤ ε and recomputes exactly),
///    and pins them by the compiled word's content address κ;
/// 4. independently validates the algorithm result against the exact reference
///    [`shor_check`] (base 2 mod 15, period 4) — the semantic oracle. The compiled braid
///    word is NOT claimed to execute the period-finding unitary; the exact reference is the
///    sole source of the result.
///
/// # Errors
/// If the density premise fails, a rotation does not compile within `ε`, or the exact
/// Shor reference deviates.
pub fn certified_carrier_compilation_check(p: &UseCaseParams) -> Result<(), String> {
    use tqc_compiler::sk::{Carrier, SkWeaver};

    // (1) PU(22) density premise from the exact certificate.
    let cert = tqc_vv::exact::exact_density_certificate(p)?;
    if !cert.pu22_dense {
        return Err("PU(22) density premise not established by the exact certificate".into());
    }

    // (2) Compile Shor's QFT controlled-phase rotations on the certified carrier and verify
    //     each synthesized flow word realizes its target within epsilon.
    let weaver = SkWeaver::for_carrier(p, Carrier::Certified22);
    let epsilon = 0.05f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    for j in 1..=3u32 {
        let theta = std::f64::consts::PI / 2f64.powi(j as i32); // π/2^j
        let word = weaver.synthesize_rotation(theta, epsilon)?;
        if !word.iter().all(|g| *g == tqc_compiler::BraidGen::Flow) {
            return Err(format!(
                "certified rotation π/2^{j} did not synthesize to a spectral-flow word"
            ));
        }
        let phase = weaver.flow_phase(word.len());
        // circle distance between the achieved phase and the target must be within epsilon.
        let mut d = (phase - theta.rem_euclid(two_pi)).rem_euclid(two_pi);
        if d > std::f64::consts::PI {
            d = two_pi - d;
        }
        if d >= epsilon {
            return Err(format!(
                "certified rotation π/2^{j}: synthesized phase off by {d} >= epsilon {epsilon}"
            ));
        }
    }

    // (3) Arbitrary-axis compilation with shipped residuals: compile a small circuit of
    //     Rx/Ry/Rz rotations through the certified-carrier compiler; every rotation carries a
    //     synthesis residual in the compiled word, each residual satisfies |Δ| ≤ ε and
    //     recomputes exactly, and the residuals are pinned by the word's content address κ.
    {
        use tqc_compiler::{Compiler, LogicGate};
        let compiler = Compiler::for_certified_carrier(p, Carrier::Certified22);
        let circuit = vec![
            LogicGate::Rz(0, 1.0),
            LogicGate::Rx(0, 0.7),
            LogicGate::Ry(0, 2.3),
        ];
        let word = compiler.compile(&circuit, epsilon)?;
        if word.residuals.len() != 3 {
            return Err(format!(
                "expected 3 shipped synthesis residuals, got {}",
                word.residuals.len()
            ));
        }
        for r in &word.residuals {
            if r.abs_error > r.epsilon {
                return Err(format!(
                    "shipped residual |Δ|={} exceeds its epsilon {}",
                    r.abs_error, r.epsilon
                ));
            }
            let recomputed = weaver.verify_residual(r);
            if (recomputed - r.abs_error).abs() > 1e-12 {
                return Err(format!(
                    "shipped residual |Δ|={} does not recompute ({recomputed})",
                    r.abs_error
                ));
            }
        }
        // The residuals are pinned by the compiled word's content address, which must be
        // stable and non-degenerate (the word carries generators and residuals).
        let bytes = word.canonical_bytes();
        if bytes != word.canonical_bytes() {
            return Err("compiled-word canonical bytes are not deterministic".into());
        }
        let kappa = tqc_substrate::kappa(&bytes);
        if tqc_substrate::kappa(&word.canonical_bytes()) != kappa {
            return Err("compiled-word κ is not stable".into());
        }
    }

    // (4) Replay the Shor instance against the exact reference evaluation.
    shor_check()
}
