# Oracle provenance

Every external validation artifact this repository checks against, with its authority,
pin, license, and checksum. The machine-readable twin is
[`../model/oracles.toml`](../model/oracles.toml); `cargo run -p xtask -- oracle-verify`
(and CI) assert the two agree and that each committed artifact matches its `sha256`.

| Oracle | Authority | Pin | License | Artifact | Verified by |
|---|---|---|---|---|---|
| `f1-atlas` | F1 — UOR Atlas (Lean 4, machine-checked, sorry-free) | tag `v0.21.0` / commit `b64df2a` | MIT | `oracles/f1/atlas-constants.json` (sha256 `fa08d368…`) | sha256 + `atlas-pin-check` (tip) + deferred lake re-derivation |
| `uor-addr-composition` | uor-addr — realized g2/f4/e6/e7/e8 composition | `0.2.0` (crates.io) | Apache-2.0 | _(linked code; via facade)_ | byte-for-byte κ reduction in witnesses |
| `holospaces-cc` | holospaces `vv` witnesses (CC-1/2/29/30) | commit `f241562…` | MIT | _(linked code; via facade)_ | facade against holospaces public API |
| `mtc-axioms` | MTC axioms (pentagon/hexagon, Yang–Baxter, SL(2,ℤ), Verlinde) | Bakalov–Kirillov; Turaev; Kitaev App. E; Rowell–Stong–Wang | — | _(executable predicates)_ | axiom predicates; unit-tested on Fibonacci/Ising |

## F1 re-derivation

The `f1-atlas` artifact is a committed snapshot of `decide`-proved theorem values from F1.
To regenerate authoritatively (the devcontainer provides Lean/elan):

```sh
cargo run -p xtask -- atlas-extract     # clones F1 @ pin, lake build, emits JSON, diffs
```

CI runs this in a dedicated, cached job so a silent F1 drift is caught. Day-to-day CI only
re-verifies the committed `sha256` (fast).
