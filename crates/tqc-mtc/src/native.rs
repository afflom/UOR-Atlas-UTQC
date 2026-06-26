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

use crate::verifier::ModularData;

/// Represents the failure to construct an Atlas-native MTC.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstructionObstruction {
    /// The `g2` composition yields signed structure constants, violating $N_{ij}^k \ge 0$.
    SignedFusionConstants,
    /// Mismatch between the 96 classes and the 24-dimensional carrier space for the $S$-matrix.
    DimensionMismatch,
}

impl core::fmt::Display for ConstructionObstruction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SignedFusionConstants => write!(f, "compose_g2_product yields signed structure constants, violating MTC nonnegative fusion"),
            Self::DimensionMismatch => write!(f, "mismatch between 96 Atlas classes and 24-dimensional S-matrix carrier"),
        }
    }
}

impl std::error::Error for ConstructionObstruction {}

/// Attempt to construct an Atlas-native MTC.
/// Always returns an obstruction under current sourced material constraints.
pub fn construct_atlas_native() -> Result<Box<dyn ModularData>, ConstructionObstruction> {
    Err(ConstructionObstruction::SignedFusionConstants)
}
