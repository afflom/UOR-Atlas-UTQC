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

use sk::SkWeaver;
use tqc_core::UseCaseParams;

/// A topological Braid Word composed of native Atlas generators.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BraidGen {
    /// The \(\sigma\) generator (fusion twist).
    Sigma,
    /// The \(\tau\) generator (S4 transposition).
    Tau,
    /// The \(\mu\) generator (conjugation mirror).
    Mu,
}

impl BraidGen {
    /// Formats the generator to a character.
    #[must_use]
    pub fn as_char(&self) -> char {
        match self {
            Self::Sigma => 'σ',
            Self::Tau => 'τ',
            Self::Mu => 'μ',
        }
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
    /// Initializes the compiler against the target topological `UseCaseParams`.
    #[must_use]
    pub fn new(params: &'a UseCaseParams) -> Self {
        Self {
            params,
            // The weaver admits only exact discrete phases: single-handle gate-set density
            // is exactly refuted (finite projective Clifford image), so no approximation
            // path exists. See `tqc-vv::exact::exact_density_certificate`.
            weaver: SkWeaver::new(params),
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
