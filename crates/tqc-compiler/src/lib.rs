//! Classical-gate front-end: schedules logical gate sequences (H, X, CNOT, T, rotations)
//! onto fixed generator-word templates over the Atlas class space.
//!
//! **Scope, precisely:** this crate emits braid words; it does **not** establish unitary
//! equivalence between a logical gate and its emitted word, and no such claim is made
//! anywhere in the V&V dictionary. Rotations are admitted only as *exact discrete phases*
//! of the generator set: the exact `Q(ζ₂₄)` certificate in `tqc-vv` refutes gate-set
//! density for the single-handle group (finite projective Clifford image, order 24), so
//! Solovay–Kitaev approximation of arbitrary angles is mathematically unavailable and the
//! weaver is constructed non-dense, matching the theorem.

pub mod qasm;
pub mod sk;

use sk::{Carrier, SkWeaver};
use tqc_core::UseCaseParams;

/// A topological Braid Word composed of native Atlas generators.
///
/// `Sigma`, `Tau`, `Mu` are the single-handle Clifford sub-alphabet (permutations of the
/// class space). `Flow` is the certified-carrier spectral-flow generator: it is available
/// only on a density-certified carrier and carries a continuous phase, so it never appears
/// in a Clifford (permutation) word.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BraidGen {
    /// The \(\sigma\) generator (fusion twist).
    Sigma,
    /// The \(\tau\) generator (S4 transposition).
    Tau,
    /// The \(\mu\) generator (conjugation mirror).
    Mu,
    /// The spectral-flow generator on a density-certified carrier (continuous phase).
    Flow,
}

impl BraidGen {
    /// Formats the generator to a character.
    #[must_use]
    pub fn as_char(&self) -> char {
        match self {
            Self::Sigma => 'σ',
            Self::Tau => 'τ',
            Self::Mu => 'μ',
            Self::Flow => 'Φ',
        }
    }

    /// Whether this generator is a Clifford (permutation) generator.
    #[must_use]
    pub fn is_clifford(&self) -> bool {
        matches!(self, Self::Sigma | Self::Tau | Self::Mu)
    }
}

/// The synthesized topological operator block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BraidWord {
    /// The linear sequence of generators.
    pub sequence: Vec<BraidGen>,
}

impl BraidWord {
    /// Creates an empty word.
    #[must_use]
    pub fn new() -> Self {
        Self {
            sequence: Vec::new(),
        }
    }

    /// Appends a generator.
    pub fn push(&mut self, gen: BraidGen) {
        self.sequence.push(gen);
    }
}

impl Default for BraidWord {
    fn default() -> Self {
        Self::new()
    }
}

/// A standard quantum logic gate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LogicGate {
    /// Hadamard gate on a target qubit.
    Hadamard(usize),
    /// Pauli-X gate.
    PauliX(usize),
    /// T gate (pi/8 phase).
    TGate(usize),
    /// Controlled-NOT.
    CNot(usize, usize),
    /// X-axis rotation parameterized by angle.
    Rx(usize, f64),
    /// Y-axis rotation parameterized by angle.
    Ry(usize, f64),
    /// Z-axis rotation parameterized by angle.
    Rz(usize, f64),
}

/// The Compiler synthezises logical circuits into topological braids.
pub struct Compiler<'a> {
    /// The target topological parameters.
    pub params: &'a UseCaseParams,
    weaver: SkWeaver,
}

impl<'a> Compiler<'a> {
    /// Initializes the single-handle **Clifford** compiler: rotations are admitted only at
    /// exact discrete phases, since single-handle gate-set density is exactly refuted (finite
    /// projective Clifford image). See `tqc-vv::exact::exact_density_certificate`.
    #[must_use]
    pub fn new(params: &'a UseCaseParams) -> Self {
        Self {
            params,
            weaver: SkWeaver::new(params),
        }
    }

    /// Initializes the compiler on a **density-certified** carrier (`Carrier::Certified22`
    /// or `Carrier::Certified576`), where the projective closure is dense (Theorem: PU(22)
    /// / PU(576) density) and arbitrary rotations are synthesized by deterministic
    /// spectral-flow search to precision `epsilon`.
    #[must_use]
    pub fn for_certified_carrier(params: &'a UseCaseParams, carrier: Carrier) -> Self {
        Self {
            params,
            weaver: SkWeaver::for_carrier(params, carrier),
        }
    }

    /// Schedules a logical circuit onto a contiguous braid word. `epsilon` is the phase
    /// tolerance for rotation gates: a rotation is admitted only if an exact discrete
    /// generator phase lies within `epsilon` of it (no dense approximation exists).
    ///
    /// # Errors
    /// Returns an error if a rotation has no exact discrete phase within `epsilon`.
    pub fn compile(&self, circuit: &[LogicGate], epsilon: f64) -> Result<BraidWord, String> {
        let mut word = BraidWord::new();

        for gate in circuit {
            match gate {
                LogicGate::Hadamard(_) => {
                    word.push(BraidGen::Sigma);
                    word.push(BraidGen::Tau);
                    word.push(BraidGen::Sigma);
                }
                LogicGate::PauliX(_) => {
                    word.push(BraidGen::Mu);
                }
                LogicGate::TGate(_) => {
                    for _ in 0..5 {
                        word.push(BraidGen::Sigma);
                        word.push(BraidGen::Tau);
                        word.push(BraidGen::Tau);
                        word.push(BraidGen::Sigma);
                    }
                }
                LogicGate::CNot(_, _) => {
                    word.push(BraidGen::Sigma);
                    word.push(BraidGen::Sigma);
                    word.push(BraidGen::Mu);
                    word.push(BraidGen::Tau);
                    word.push(BraidGen::Sigma);
                    word.push(BraidGen::Tau);
                }
                LogicGate::Rx(_, theta) | LogicGate::Ry(_, theta) | LogicGate::Rz(_, theta) => {
                    let seq = self.weaver.synthesize_rotation(*theta, epsilon)?;
                    for gen in seq {
                        word.push(gen);
                    }
                }
            }
        }
        Ok(word)
    }
}
