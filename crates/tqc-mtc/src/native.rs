//! The Atlas-native modular tensor category construction.
//!
//! This module attempts to construct an Atlas-native MTC from the sourced Atlas material:
//! the 96 classes, the 24-dimensional carrier `V_T ⊗ V_O`, the `g2` composition, and the
//! reflection generators.
//!
//! # Obstruction
//!
//! Currently, **no coherent Atlas-native MTC can be built** from the sourced material,
//! for the following structural reasons:
//!
//! 1. **Signed Structure Constants vs. Non-negative Fusion**: The MTC axioms require fusion
//!    coefficients $N_{ij}^k$ to be non-negative integers ($\in \mathbb{Z}_{\ge 0}$). However,
//!    the `compose_g2_product` is derived from a normed division algebra (the octonion 8-square
//!    over the carrier). The structure constants of this algebra contain negative signs (e.g.,
//!    $e_1 \cdot e_2 = e_3$, $e_2 \cdot e_1 = -e_3$). Thus, `g2` cannot serve directly as a
//!    categorical fusion ring without a major structural transformation.
//! 2. **Dimension Mismatch**: The 96 Atlas labels outnumber the 24 dimensions of the carrier
//!    $V_T \otimes V_O$. If the labels are simple objects, the modular $S$ matrix must be $96 \times 96$.
//!    If the carrier dimensions are the simple objects, the matrix is $24 \times 24$, but this
//!    leaves the 96 classes as derived or composite structures rather than simple objects.
//!
//! Because of this obstruction, `D(Z_O)` remains the explicitly designated generic representative
//! stand-in. The `verify_mtc_axioms` oracle from the `verifier` module would reject the `g2`
//! structure constants due to the non-negative integer requirement.
//!
//! **Future Work as Conditional Research:**
//! Any future Atlas-native category construction is treated purely as an experimental research branch
//! contingent on either new external source material, a newly derived simple-object basis, or a
//! rigorous mathematical structural transformation linking Atlas composition to valid categorical fusion.

use crate::verifier::ModularData;
use tqc_core::params::UseCaseParams;
use tqc_core::spectrum;

/// Represents the failure to construct an Atlas-native MTC.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstructionObstruction {
    /// Mismatch between the classes count and the carrier space dimensions for the $S$-matrix.
    DimensionMismatch(u64, u64),
    /// The spectral operator is indefinite, obstructing a unitary S-matrix.
    IndefiniteSpectralSignature,
    /// The `g2` composition yields signed structure constants, violating $N_{ij}^k \ge 0$.
    SignedFusionConstants,
}

impl core::fmt::Display for ConstructionObstruction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DimensionMismatch(classes, carrier) => write!(f, "mismatch between {classes} classes and {carrier}-dimensional S-matrix carrier"),
            Self::IndefiniteSpectralSignature => write!(f, "the spectral operator is indefinite, obstructing a unitary S-matrix"),
            Self::SignedFusionConstants => write!(f, "compose_g2_product yields signed structure constants, violating MTC nonnegative fusion"),
        }
    }
}

impl std::error::Error for ConstructionObstruction {}

/// Attempt to construct an Atlas-native MTC from parameters.
/// Always returns an obstruction under current sourced material constraints.
pub fn construct_atlas_native(
    p: &UseCaseParams,
) -> Result<Box<dyn ModularData>, ConstructionObstruction> {
    if p.class_count() != p.carrier_dim() {
        return Err(ConstructionObstruction::DimensionMismatch(
            p.class_count(),
            p.carrier_dim(),
        ));
    }

    let ev = spectrum::block_eigenvalues(p);
    if ev.iter().any(|&e| e < 0) {
        return Err(ConstructionObstruction::IndefiniteSpectralSignature);
    }

    let sc = tqc_core::octonion::structure_constants(8);
    if sc.iter().any(|&(_, _, _, val)| val < 0) {
        return Err(ConstructionObstruction::SignedFusionConstants);
    }

    // Fallback if somehow there is another obstruction or we haven't built the true result.
    // However, the above check will always trip for octonions because e_1 * e_2 = e_3 but e_2 * e_1 = -e_3.
    unreachable!("Atlas MTC construction succeeded despite obstructions");
}
