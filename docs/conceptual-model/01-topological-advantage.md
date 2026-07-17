# 01 Measured Content-Addressed Deduplication (Topological Degeneracy)

The substrate's measured efficiency is **content-addressed deduplication of topological degeneracy via UOR cache-collapse**, over the finite modular sector. It is an engineering measurement, reported via `just report`; it is never asserted as a proven quantum advantage.

When evaluating braids, isotopic generator words (distinct paths that compose to the same topological operator) evaluate to the same quantum state. Rather than tracking an exponentially growing tree of paths, Holospaces hashes each evaluated state to a `κ` address via its Universal Object Reference (UOR) addressing. Because multiple paths collapse to the identical `κ`, the physical hardware (L1/L2/L3 caches) naturally absorbs the degeneracy.

This is a real, measurable efficiency for computations that carry topological degeneracy.

## The Finite-Image Plateau vs Universal Computation

It is critical to bound this claim appropriately: **this is measured content-addressed deduplication (topological degeneracy) over the finite modular sector — not a proven quantum advantage.**

While a generic, dense $N$-qubit quantum circuit classically tracking an arbitrary $O(2^N)$ continuous state vector requires exponential memory, the workloads measured here operate over a *finite* braid orbit: the number of distinct operators plateaus (at 32 in the trace below), so the exponentially many braid words deduplicate onto a bounded set of κ addresses. The deduplication ratio is measured and reported; no claim is made that this beats classical simulation of universal quantum computation, and the finite (Clifford-like) closure of the measured orbit is precisely why the plateau exists.

## Empirical Benchmark: The Degeneracy Plateau

A scaling trace executed at deep braiding networks over the generators $\{\sigma, \tau, \mu\}$ measures exactly this cache-collapse deduplication when operating over a finite orbit.

| Braid Depth | Total Paths | Distinct Operators ($\kappa$) | Cache Hit Elision % |
|-------------|-------------|-------------------------------|---------------------|
| 2           | 9           | 9                             | 0.00%               |
| 4           | 81          | 27                            | 66.67%              |
| 6           | 729         | 32                            | 95.61%              |
| 8           | 6,561       | 32                            | 99.51%              |
| 10          | 59,049      | 32                            | 99.95%              |

*The plateau of unique `κ` addresses (at 32 distinct states in this specific finite orbit) demonstrates the closure of the group image. As the evaluation scales beyond $N=6$, the number of possible combinatorial paths grows exponentially, while the topological evaluations perfectly flat-line in memory due to content deduplication. This content-addressed elision provides significant performance multipliers specifically for workloads possessing deep topological repetition.*

## Isotopy κ-Collision (decision by equivalence)

Crucially, standard quantum complexity theory bounds classical emulation via the `#P-hard` tensor contraction problem—the exponentially expensive task of extracting continuous scalar probabilities by projecting the final global state against the genesis state.

The UOR Atlas sidesteps that particular cost model by framing algorithmic outputs as **Topological Decision Problems**. Because Holospaces' UOR cache-collapse guarantees that all distinct computational histories resulting in isotopic topologies collapse to identical cryptographic $\kappa$-addresses, the system never needs to perform a numeric tensor contraction over $\mathbb{C}$ *for these workloads*. The final invariants can be interrogated directly via their topological identity. Thus, if the inputs, operations, and outputs remain exclusively as $\kappa$-addressed $k$-forms, the cost of continuous amplitude extraction does not arise for finite topological compilation and evaluation — a measured deduplication over the finite modular sector, not a proven speedup for general quantum simulation (which would require escaping the finite closure this measurement relies on).
