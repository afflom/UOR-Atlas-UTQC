# TQC task runner. Mirrors the holospaces substrate's Justfile conventions and adds
# the docs-as-code / BDD / V&V targets. `just` with no args lists everything.

default:
    @just --list

# --- Quality gates (each is also a CI gate) ---

# Format check (CI uses --check; locally `just fmt-fix` rewrites).
fmt:
    cargo fmt --all --check
fmt-fix:
    cargo fmt --all

# Deny-heavy clippy across all targets/features. Warnings are errors.
lint:
    cargo clippy --workspace --all-targets --all-features -- -D warnings

# Unit + integration tests across the workspace.
test:
    cargo test --workspace --all-features

# Rustdoc must build clean (docs-as-code).
doc:
    RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --all-features

# --- BDD / conceptual-model / V&V ---

# Run the Gherkin BDD suite (cucumber). The runner itself fails on any skipped/undefined step.
bdd:
    cargo test -p tqc-conformance --test bdd

# The conceptual-model honesty gate: every dictionary row has a feature + witness;
# status discipline holds; no open/none claim is asserted as established.
honesty:
    cargo test -p tqc-conformance --test honesty -- --nocapture

# Verify every external oracle artifact against its committed checksum + provenance.
oracles:
    cargo run -p xtask -- oracle-verify

# Confirm the pinned F1 commit exists upstream and the artifact digest matches (online).
atlas-pin-check:
    cargo run -p xtask -- atlas-pin-check

# Emit the conformance ledger (suite witnesses against the F1 oracle).
report:
    cargo run -p xtask -- report

# Build + test the substrate facade (excluded crate; real holospaces/hologram/uor-addr).
substrate:
    cargo test --manifest-path crates/tqc-substrate/Cargo.toml

# Supply-chain / license gate.
deny:
    cargo deny check

# Build the academic whitepaper strictly with warnings-as-errors and run chktex linting
paper:
    cd docs/paper && latexmk -pdf -Werror -interaction=nonstopmode main.tex
    cd docs/paper && chktex main.tex -q -n 22 -n 30 -n 46

# The full local gate (what CI runs).
vv: fmt lint test doc bdd honesty oracles paper
    @echo "V&V: all gates green."
