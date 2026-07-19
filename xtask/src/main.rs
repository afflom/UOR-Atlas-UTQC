//! Workspace automation.
//!
//! - `oracle-verify`    — every committed oracle artifact matches its `model/oracles.toml` sha256 (offline).
//! - `atlas-pin-check`  — the F1 pin exists upstream and the artifact digest matches (online).
//! - `report`           — run the suite witnesses and emit a conformance ledger.

#![forbid(unsafe_code)]

use anyhow::{anyhow, bail, Context, Result};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};
use tqc_atlas::canonical;
use tqc_model::{Model, Tier};
use tqc_vv::{witness, F1Constants};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn sha256_file(path: &Path) -> Result<String> {
    let bytes = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let mut h = Sha256::new();
    h.update(&bytes);
    Ok(hex::encode(h.finalize()))
}

fn main() -> Result<()> {
    let cmd = std::env::args().nth(1).unwrap_or_default();
    match cmd.as_str() {
        "oracle-verify" => oracle_verify(),
        "atlas-pin-check" => atlas_pin_check(),
        "report" => report(),
        "paper-check" => paper_check(),
        other => bail!(
            "unknown command `{other}`; use: oracle-verify | atlas-pin-check | report | paper-check"
        ),
    }
}

/// The maintained paper-claim ⇄ dictionary-row crosswalk (`docs/paper/claims_crosswalk.toml`).
#[derive(serde::Deserialize)]
struct Crosswalk {
    #[serde(default)]
    assert: Vec<AssertClaim>,
    #[serde(default)]
    forbid: Vec<ForbidClaim>,
    #[serde(default)]
    global: GlobalForbid,
}

#[derive(serde::Deserialize)]
struct AssertClaim {
    claim: String,
    row: String,
}

#[derive(serde::Deserialize)]
struct ForbidClaim {
    row: String,
    phrases: Vec<String>,
}

#[derive(serde::Deserialize, Default)]
struct GlobalForbid {
    #[serde(default)]
    forbidden: Vec<String>,
}

/// The honesty gate, extended to the paper. Mechanical and conservative:
///
/// 1. Every `[[assert]]` entry names a real dictionary row whose status is **assertable**
///    (`some-true` / `build`, i.e. gating). A row demoted to `open` under an asserted claim
///    is a status regression and fails the gate.
/// 2. Every `open` model row must appear in a `[[forbid]]` block, so a newly-open row cannot
///    silently escape the scan.
/// 3. No `[[forbid]]` phrase appears in the paper as an **affirmative** assertion — a
///    negation window (`not`, `never`, `no`, `without`, `rather than`, `cannot`) before the
///    phrase marks a permitted disclaimer.
/// 4. No `[global].forbidden` phrase appears anywhere (phrases with no legitimate use).
fn paper_check() -> Result<()> {
    let model = Model::load().map_err(|e| anyhow!(e.to_string()))?;
    let root = root();
    let cross_path = root.join("docs/paper/claims_crosswalk.toml");
    let cross_text =
        fs::read_to_string(&cross_path).context("reading docs/paper/claims_crosswalk.toml")?;
    let cross: Crosswalk = toml::from_str(&cross_text).context("parsing claims_crosswalk.toml")?;

    // Concatenate the paper section sources.
    let sections_dir = root.join("docs/paper/sections");
    let mut paper = String::new();
    let mut files: Vec<_> = fs::read_dir(&sections_dir)
        .context("reading docs/paper/sections")?
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("tex"))
        .collect();
    files.sort();
    for f in &files {
        paper.push_str(&fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?);
        paper.push('\n');
    }
    let paper_lc = paper.to_lowercase();

    // (1) Asserted claims must name assertable rows.
    for a in &cross.assert {
        let row = model
            .rows
            .iter()
            .find(|r| r.id == a.row)
            .ok_or_else(|| anyhow!("crosswalk asserts unknown row `{}`", a.row))?;
        let status = model
            .status(&row.status)
            .ok_or_else(|| anyhow!("row `{}` has unknown status `{}`", row.id, row.status))?;
        if !status.gating {
            bail!(
                "paper asserts claim \"{}\" backed by row `{}`, but that row is now non-gating \
                 (status `{}`): status regression",
                a.claim,
                a.row,
                row.status
            );
        }
    }

    // (2) Every open (non-gating) model row must have a forbid block.
    for row in &model.rows {
        let status = model
            .status(&row.status)
            .ok_or_else(|| anyhow!("row `{}` has unknown status `{}`", row.id, row.status))?;
        if !status.gating && !cross.forbid.iter().any(|f| f.row == row.id) {
            bail!(
                "open row `{}` has no [[forbid]] block in the crosswalk: its claim could be \
                 asserted in the paper unchecked",
                row.id
            );
        }
    }

    // (3) No forbidden open-row phrase appears as an affirmative assertion.
    for f in &cross.forbid {
        for phrase in &f.phrases {
            for pos in find_all(&paper_lc, &phrase.to_lowercase()) {
                if !negated_before(&paper_lc, pos) {
                    bail!(
                        "paper affirmatively asserts open-row `{}` claim phrase \"{}\" (no \
                         negation in the preceding window); disclaim it or demote the claim",
                        f.row,
                        phrase
                    );
                }
            }
        }
    }

    // (4) Globally-forbidden phrases must not appear at all.
    for phrase in &cross.global.forbidden {
        if paper_lc.contains(&phrase.to_lowercase()) {
            bail!("paper contains globally-forbidden phrase \"{phrase}\" (removed at a159aed)");
        }
    }

    println!(
        "paper-check: {} asserted claims backed by gating rows; {} open row(s) guarded; \
         {} globally-forbidden phrase(s) absent",
        cross.assert.len(),
        cross.forbid.len(),
        cross.global.forbidden.len()
    );
    Ok(())
}

/// All byte offsets at which `needle` occurs in `hay`.
fn find_all(hay: &str, needle: &str) -> Vec<usize> {
    let mut out = Vec::new();
    let mut from = 0;
    while let Some(rel) = hay[from..].find(needle) {
        let at = from + rel;
        out.push(at);
        from = at + needle.len();
    }
    out
}

/// Whether a negation token appears within the ~40-char window before `pos`.
fn negated_before(hay: &str, pos: usize) -> bool {
    let start = pos.saturating_sub(40);
    let window = &hay[start..pos];
    [
        "not ",
        "never ",
        "no ",
        "without ",
        "rather than ",
        "cannot ",
        "nor ",
    ]
    .iter()
    .any(|tok| window.contains(tok))
}

/// Verify every committed oracle artifact against its recorded sha256, and that the
/// human-readable provenance table (`oracles/PROVENANCE.md`) agrees with the
/// machine-readable manifest: every oracle id must appear in the table, and every
/// artifact-bearing oracle's sha256 prefix must appear beside it.
fn oracle_verify() -> Result<()> {
    let model = Model::load().map_err(|e| anyhow!(e.to_string()))?;
    let root = root();
    let provenance = fs::read_to_string(root.join("oracles/PROVENANCE.md"))
        .context("reading oracles/PROVENANCE.md")?;
    let mut checked = 0u32;
    for o in &model.oracles {
        if !provenance.contains(&format!("`{}`", o.id)) {
            bail!(
                "oracle `{}` is missing from oracles/PROVENANCE.md; the provenance table \
                 must agree with model/oracles.toml",
                o.id
            );
        }
        if o.artifact.is_empty() {
            continue;
        }
        let got = sha256_file(&root.join(&o.artifact))?;
        if got != o.sha256 {
            bail!(
                "oracle `{}` sha256 mismatch: {got} != {} (manifest)",
                o.id,
                o.sha256
            );
        }
        if !provenance.contains(&got[..8]) {
            bail!(
                "oracle `{}`: sha256 prefix {} missing from oracles/PROVENANCE.md",
                o.id,
                &got[..8]
            );
        }
        println!("ok   {:<22} {}  {}…", o.id, o.artifact, &got[..12]);
        checked += 1;
    }
    println!(
        "oracle-verify: {checked} artifact(s) verified; PROVENANCE.md agrees with \
         model/oracles.toml on all {} oracles",
        model.oracles.len()
    );
    Ok(())
}

/// Confirm the pinned F1 commit exists upstream and the artifact digest matches.
fn atlas_pin_check() -> Result<()> {
    let model = Model::load().map_err(|e| anyhow!(e.to_string()))?;
    let o = model
        .oracle("f1-atlas")
        .context("model is missing the f1-atlas oracle")?;
    let got = sha256_file(&root().join(&o.artifact))?;
    if got != o.sha256 {
        bail!("F1 artifact sha256 mismatch: {got} != {}", o.sha256);
    }
    let out = std::process::Command::new("git")
        .args(["ls-remote", o.source.as_str()])
        .output()
        .with_context(|| format!("git ls-remote {}", o.source))?;
    if !out.status.success() {
        bail!(
            "git ls-remote {} failed: {}",
            o.source,
            String::from_utf8_lossy(&out.stderr)
        );
    }
    // `git ls-remote` lists ref tips as `<sha>\t<ref>`. Require the pin to be a live tip
    // (line-start match), so the pin must be an immutable tag/branch head — not merely some
    // historical commit that ls-remote cannot see.
    let listing = String::from_utf8_lossy(&out.stdout);
    let tip = listing.lines().find(|l| l.starts_with(o.pin.as_str()));
    match tip {
        Some(line) => {
            let reference = line.split('\t').nth(1).unwrap_or("<unknown ref>");
            println!(
                "atlas-pin-check: F1 pin {} is a live upstream tip ({}); artifact digest ok",
                o.pin, reference
            );
            println!("note: authoritative numeric re-derivation (lake build + Extract.lean) is the deferred CI Lean job");
            Ok(())
        }
        None => bail!(
            "pinned F1 commit {} is not a live ref tip at {} (pin a release tag, not a transient commit)",
            o.pin,
            o.source
        ),
    }
}

/// Run the suite witnesses and emit a conformance ledger.
fn report() -> Result<()> {
    let model = Model::load().map_err(|e| anyhow!(e.to_string()))?;
    let f1 = F1Constants::load().map_err(|e| anyhow!(e))?;
    let p = canonical(&model)?;

    let mut lines = Vec::new();
    lines.push("# Conformance ledger\n".to_owned());
    lines.push(format!(
        "F1 oracle digest: `{}…`\n",
        &F1Constants::sha256()[..16]
    ));
    lines.push("| row | status | tier | stage | result |".to_owned());
    lines.push("|---|---|---|---|---|".to_owned());

    let mut suites = 0u32;
    let mut passed = 0u32;
    for row in &model.rows {
        let result = match row.tier {
            Tier::Suite => {
                suites += 1;
                match run_suite_witness(&row.id, &model, &p, &f1) {
                    Some(Ok(())) => {
                        passed += 1;
                        "PASS".to_owned()
                    }
                    Some(Err(e)) => format!("FAIL: {e}"),
                    None => "NO WITNESS".to_owned(),
                }
            }
            Tier::Target => match row.id.as_str() {
                // Open rows are measured and reported, never asserted.
                "advantage" => match witness::advantage_probe(&p) {
                    Ok(m) => format!(
                        "MEASURED (non-gating): degeneracy {:.3}x over {} paths / {} distinct κ; \
                         compute savings {:.2}%; storage compression {:.2}x ({} → {} bytes)",
                        m.topological_degeneracy,
                        m.total_paths,
                        m.distinct_states,
                        m.compute_savings_pct,
                        m.memory_compression_ratio,
                        m.total_state_bytes,
                        m.unique_state_bytes
                    ),
                    Err(e) => format!("MEASUREMENT FAILED (non-gating): {e}"),
                },
                _ => "target (expected-RED, non-gating)".to_owned(),
            },
        };
        lines.push(format!(
            "| {} | {} | {:?} | {} | {} |",
            row.id, row.status, row.tier, row.stage, result
        ));
    }
    lines.push(format!("\nsuites: {passed}/{suites} passed"));

    let cert_str = match tqc_vv::exact::exact_density_certificate(&p) {
        Ok(m) => m.description,
        Err(e) => format!("certificate error: {e}"),
    };
    lines.push(format!("\nexact certificate: {cert_str}\n"));

    let report = lines.join("\n");
    let out_dir = root().join("target/conformance");
    fs::create_dir_all(&out_dir).context("creating target/conformance")?;
    fs::write(out_dir.join("ledger.md"), &report).context("writing ledger.md")?;
    println!("{report}");
    if passed != suites {
        bail!("conformance: {passed}/{suites} suites passed");
    }
    Ok(())
}

fn run_suite_witness(
    id: &str,
    model: &Model,
    p: &tqc_core::UseCaseParams,
    f1: &F1Constants,
) -> Option<Result<(), String>> {
    Some(match id {
        "objects-labels" => witness::objects_labels(p, f1),
        "label-space-belt" => witness::label_space_belt(p, f1),
        "inner-product" => witness::inner_product(p),
        "reflection-generators" => witness::reflection_generators(p, f1),
        "spectrum" => witness::spectrum(p, f1),
        "coxeter-weyl" => witness::coxeter_weyl(p, f1),
        "modular-identities" => witness::modular_identities(p, f1),
        "definite-anchor-e8" => witness::definite_anchor_e8(f1),
        "fusion-g2" => witness::fusion_g2(p),
        "dual-f4" => witness::dual_f4(p),
        "categorical-structure" => witness::categorical_structure(p),
        "ground-space-protection" => witness::ground_space_protection(p),
        "complex-amplitude-encoding" => witness::complex_amplitude_encoding(p, f1),
        "modular-s-t" => witness::modular_s_t(p),
        "braiding-r-matrix" => witness::braiding_r_matrix(p),
        "holospace-cycle" => witness::holospace_cycle(p),
        "atlas-native-mtc" => witness::atlas_native_mtc(p),
        "generative-closure" => witness::generative_closure_probe(p),
        "utqc-proven" => witness::utqc_proven_probe(model, f1, p),
        "quantum-realization" => witness::quantum_realization(p),
        "universality" => witness::equivalency_universality_probe(p),
        "two-qubit-universality" => witness::two_qubit_universality_probe(p).and_then(|m| {
            if m.is_entangling && m.is_coherent {
                Ok(())
            } else {
                Err("no coherent native entangler decided".into())
            }
        }),
        "solovay-kitaev" => witness::solovay_kitaev_decision_witness(p),
        "archimedean-continuity" => witness::archimedean_continuity_witness(p),
        "pair-carrier-structure" => witness::pair_carrier_witness(p),
        "encoded-qubit-universality" => witness::encoded_qubit_universality_witness(p),
        "certified-carrier-compilation" => {
            tqc_algorithms::checks::certified_carrier_compilation_check(p)
        }
        "eps-free-decision-path" => witness::eps_free_decision_path_witness(p),
        "canonical-kappa-form" => witness::canonical_kappa_witness(p),
        "reduction-crux" => witness::reduction_crux_witness(p),
        "contraction-class-elimination" => witness::contraction_class_elimination_witness(p),
        "finite-closure" => witness::finite_closure_probe(p).map(|_| ()),
        "s4-modal-logic" => witness::s4_frame_witness(p),
        "mac-lane-coherence" => witness::mac_lane_coherence(p),
        "fault-tolerance" => witness::deterministic_replay_witness(p),
        "complexity-bound" => witness::complexity_bound_witness(p),
        "reconstructability" => witness::reconstruction_witness(p),
        "tensor-contraction-bypass" => witness::isotopy_collision_witness(p),
        "topological-entanglement" => witness::topological_entanglement_probe(p).and_then(|m| {
            if m.is_logarithmic_scaling && m.entropy_bound > 0.0 {
                Ok(())
            } else {
                Err(format!(
                    "measured entanglement profile fails the derived verdict: {:?}",
                    m.depth_profile
                ))
            }
        }),
        "grover-search" => tqc_algorithms::checks::grover_check(),
        "qft-algorithm" => tqc_algorithms::checks::qft_word_check(p).map(|_| ()),
        "qpe-algorithm" => tqc_algorithms::checks::qpe_check(),
        "shor-algorithm" => tqc_algorithms::checks::shor_check(),
        "whitepaper-formatting" => whitepaper_formatting_check(),
        _ => return None,
    })
}

/// The whitepaper-formatting row: docs/paper/main.tex conforms to the RevTeX 4-2 spec
/// (document class) and carries the tikz visual-aid dependency — the same checks the BDD
/// steps make.
fn whitepaper_formatting_check() -> Result<(), String> {
    let path = root().join("docs/paper/main.tex");
    let text = fs::read_to_string(&path).map_err(|e| format!("reading main.tex: {e}"))?;
    if !text.contains("{revtex4-2}") {
        return Err("main.tex does not use the revtex4-2 document class".into());
    }
    if !text.contains("\\usepackage{tikz}") {
        return Err("main.tex does not include tikz".into());
    }
    Ok(())
}
