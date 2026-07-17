@row:advantage @stage:S4 @status:open @oracle:holospaces-cc
Feature: Advantage metrics (measured, open)
  An open target: the content-addressed deduplication (topological degeneracy) of braid
  evaluation over the finite modular sector is measured and reported via `just report`.
  It is never asserted as a proven quantum advantage; this feature is non-gating and is
  not executed by the BDD suite runner.

  Scenario: the deduplication metrics are measured and reported
    Given the UOR Atlas use-case
    Then the cache-collapse deduplication metrics are measured and reported without any advantage assertion
