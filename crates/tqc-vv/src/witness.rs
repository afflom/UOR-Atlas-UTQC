//! The V&V witnesses. Each binds **one** parametric computation (from `tqc-core`) to **one**
//! authoritative oracle value (from F1), returning a readable error on mismatch.
//!
//! These functions are the single implementation of each check; both the `#[test]`s below and
//! the cucumber step definitions in `tqc-conformance` call them (DRY).

use crate::oracle::F1Constants;
use tqc_core::generators::Generators;
use tqc_core::inner::{euclidean_norm_sq, preserves_norm};
use tqc_core::{coxeter, labels, modular, spectrum, UseCaseParams};
use tqc_model::Model;

/// Outcome of a witness.
pub type Witness = Result<(), String>;

fn check(cond: bool, msg: impl Into<String>) -> Witness {
    if cond {
        Ok(())
    } else {
        Err(msg.into())
    }
}

/// VV-00 — the committed F1 artifact matches its recorded pin, and was extracted at the
/// canonical parameters.
///
/// # Errors
/// If the sha256 differs from `model/oracles.toml`, or the parameters disagree.
pub fn oracle_provenance(model: &Model, f1: &F1Constants) -> Witness {
    let oracle = model
        .oracle("f1-atlas")
        .ok_or_else(|| "model is missing the `f1-atlas` oracle".to_owned())?;
    check(
        F1Constants::sha256() == oracle.sha256,
        format!(
            "F1 artifact sha256 {} != manifest {}",
            F1Constants::sha256(),
            oracle.sha256
        ),
    )?;
    let uc = model
        .canonical_usecase()
        .ok_or_else(|| "no canonical use-case".to_owned())?;
    check(
        (f1.params.scope_q, f1.params.modality_t, f1.params.context_o)
            == (uc.scope, uc.modality, uc.context),
        "F1 oracle parameters disagree with the canonical use-case",
    )
}

/// VV — objects / anyon labels: class count, stride, and the `classIndex` bijection.
///
/// # Errors
/// On any mismatch with F1 or a non-bijective index.
pub fn objects_labels(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(p.class_count() == f1.classes.count, "class count != F1")?;
    check(p.stride() == f1.classes.stride, "stride != F1")?;
    check(
        labels::class_index_is_bijection(p),
        "classIndex is not a bijection onto [0,count)",
    )
}

/// VV — label / state-space index: the belt extent and its factorizations.
///
/// # Errors
/// On any mismatch with F1 or a non-contiguous belt.
pub fn label_space_belt(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(
        p.belt_extent() == f1.classes.belt_extent,
        "belt extent != F1",
    )?;
    let mut got: Vec<(u64, u64)> = p.belt_factorizations();
    let mut want: Vec<(u64, u64)> = f1
        .classes
        .belt_factorizations
        .iter()
        .map(|v| (v[0], v[1]))
        .collect();
    got.sort_unstable();
    want.sort_unstable();
    check(got == want, "belt factorizations != F1")?;
    let (pages, page) = (want[0].0.min(want[0].1), want[0].0.max(want[0].1));
    let addrs = labels::belt_addresses(pages, page);
    check(
        addrs.len() as u64 == p.belt_extent() && addrs.last() == Some(&(p.belt_extent() - 1)),
        "belt addresses are not contiguous over the extent",
    )
}

/// VV — the inner product is the definite Euclidean companion `Σxᵢ²`.
///
/// Validates definiteness (not merely PSD), keeping it distinct from the signed prime form.
///
/// # Errors
/// If the form is not a definite sum of squares.
pub fn inner_product(p: &UseCaseParams) -> Witness {
    let n = p.carrier_dim() as usize;
    let zero = vec![0i64; n];
    check(
        euclidean_norm_sq(&zero) == 0,
        "norm of the zero vector must be 0",
    )?;
    let mut v = zero.clone();
    if let Some(first) = v.first_mut() {
        *first = 1;
    }
    check(
        euclidean_norm_sq(&v) > 0,
        "norm of a non-zero vector must be positive (definite)",
    )?;
    check(
        euclidean_norm_sq(&[3, 4]) == 25,
        "norm must be the plain sum of squares",
    )
}

/// VV — the reflection generators: orders match F1, and each is orthogonal on `Σxᵢ²`.
///
/// This is "the unlock": orthogonality (= unitarity) holds with no positivity assumption.
///
/// # Errors
/// On an order mismatch or a generator that fails to preserve the norm.
pub fn reflection_generators(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let g = Generators::new(p);
    check(
        g.sigma.order() == u64::from(f1.generators.sigma_order),
        "sigma order != F1",
    )?;
    check(
        g.tau.order() == u64::from(f1.generators.tau_order),
        "tau order != F1",
    )?;
    check(
        g.mu.order() == u64::from(f1.generators.mu_order),
        "mu order != F1",
    )?;
    // Orthogonality over the whole class space.
    let n = p.class_count() as usize;
    let v: Vec<i64> = (0..n as i64).map(|i| i % 11 - 5).collect();
    for (name, perm) in [("sigma", &g.sigma), ("tau", &g.tau), ("mu", &g.mu)] {
        check(
            preserves_norm(perm, &v),
            format!("generator {name} does not preserve the norm"),
        )?;
    }
    Ok(())
}

/// VV — the spectrum: parametric block eigenvalues plus F1 multiplicities reconcile to the
/// F1 signature and trace.
///
/// # Errors
/// On any disagreement with F1.
pub fn spectrum(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(
        spectrum::block_eigenvalues(p).as_slice() == f1.spectrum.eigenvalues.as_slice(),
        "block eigenvalues != F1",
    )?;
    let sig = spectrum::reconcile(p, &f1.spectrum.eigenvalues, &f1.spectrum.multiplicities)
        .map_err(|e| format!("spectrum reconciliation failed: {e:?}"))?;
    check(
        f1.spectrum.signature == vec![sig.positive, sig.negative],
        format!(
            "signature ({},{}) != F1 {:?}",
            sig.positive, sig.negative, f1.spectrum.signature
        ),
    )?;
    check(
        i64::try_from(p.carrier_dim()) == Ok(f1.spectrum.trace),
        "trace != carrier dim",
    )?;
    check(p.carrier_dim() == f1.spectrum.dim, "dim != F1")?;
    check(
        f1.spectrum.indefinite,
        "F1 records the operator as indefinite",
    )
}

/// VV — Coxeter / Weyl: `rank = φ(h) = context`, and the exponents are consistent.
///
/// # Errors
/// On any disagreement with F1.
pub fn coxeter_weyl(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let rank = coxeter::euler_phi(f1.coxeter.number_h);
    check(rank == f1.coxeter.rank, "phi(h) != F1 rank")?;
    check(rank == p.context, "rank != context (O)")?;
    check(
        f1.coxeter.exponents.len() as u32 == rank,
        "exponent count != rank",
    )?;
    check(
        f1.coxeter.exponents.iter().sum::<u32>() == f1.coxeter.exponent_sum,
        "exponent sum != F1",
    )
}

/// VV — the definite anchor: the E8 Gram is `4 × Cartan` (diag 8, edges -4) and is
/// positive-definite, matching the F1 `e8_seed`.
///
/// # Errors
/// On any disagreement with F1 or a non-PD Gram.
pub fn definite_anchor_e8(f1: &F1Constants) -> Witness {
    let scale = f1.e8_seed.gram_scale;
    let cartan = tqc_atlas::e8_cartan();
    for (i, row) in cartan.iter().enumerate() {
        check(row[i] == f1.e8_seed.cartan_diag, "E8 Cartan diagonal != F1")?;
    }
    let gram = tqc_atlas::e8_gram(scale);
    for i in 0..8 {
        check(
            gram[i][i] == f1.e8_seed.gram_diag,
            "E8 Gram diagonal != F1 gram_diag",
        )?;
        for j in 0..8 {
            check(
                gram[i][j] == scale * cartan[i][j],
                "E8 Gram != scale*Cartan",
            )?;
            if i != j && gram[i][j] != 0 {
                check(
                    gram[i][j] == f1.e8_seed.gram_edge,
                    "E8 Gram edge != F1 gram_edge",
                )?;
            }
        }
    }
    check(
        tqc_core::anchor::is_positive_definite(&gram) == f1.e8_seed.psd,
        "E8 positive-definiteness != F1",
    )?;
    check(f1.e8_seed.psd, "F1 records the E8 seed as PSD")
}

/// VV — the generic definite anchor: the use-case's Euclidean companion is positive-definite.
///
/// # Errors
/// If the companion is not positive-definite.
pub fn definite_anchor(p: &UseCaseParams) -> Witness {
    let gram = tqc_core::anchor::euclidean_companion(p.carrier_dim() as usize);
    check(
        tqc_core::anchor::is_positive_definite(&gram),
        "the use-case Euclidean companion must be positive-definite",
    )
}

/// VV — the modular identity `E4³ = E6² + 1728·Δ`, plus the weight `T·O/2`.
///
/// # Errors
/// If the identity fails on the F1 coefficients or the weight is inconsistent.
pub fn modular_identities(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let e4: Vec<i128> = f1.modular.e4.iter().map(|&x| i128::from(x)).collect();
    let e6: Vec<i128> = f1.modular.e6.iter().map(|&x| i128::from(x)).collect();
    let delta: Vec<i128> = f1.modular.delta.iter().map(|&x| i128::from(x)).collect();
    check(
        modular::identity_holds(&e4, &e6, &delta, i128::from(f1.modular.constant)),
        "E4^3 = E6^2 + 1728*Delta failed on the F1 coefficients",
    )?;
    check(
        u64::from(f1.modular.weight) * 2 == p.carrier_dim(),
        "weight*2 != carrier dim (T*O)",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use tqc_atlas::canonical;

    fn atlas() -> (Model, F1Constants, UseCaseParams) {
        let model = Model::load().unwrap();
        let f1 = F1Constants::load().unwrap();
        let p = canonical(&model).unwrap();
        (model, f1, p)
    }

    #[test]
    fn vv_oracle_provenance() {
        let (m, f1, _) = atlas();
        oracle_provenance(&m, &f1).unwrap();
    }

    #[test]
    fn vv_all_some_true_suite_witnesses_pass_on_the_atlas() {
        let (_, f1, p) = atlas();
        objects_labels(&p, &f1).unwrap();
        label_space_belt(&p, &f1).unwrap();
        inner_product(&p).unwrap();
        reflection_generators(&p, &f1).unwrap();
        spectrum(&p, &f1).unwrap();
        coxeter_weyl(&p, &f1).unwrap();
        modular_identities(&p, &f1).unwrap();
        definite_anchor_e8(&f1).unwrap();
        definite_anchor(&p).unwrap();
    }

    #[test]
    fn definite_anchor_holds_at_an_arbitrary_use_case() {
        let p = UseCaseParams::new(2, 2, 4);
        definite_anchor(&p).unwrap();
    }
}
