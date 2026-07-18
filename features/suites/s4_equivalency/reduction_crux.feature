@row:reduction-crux @stage:S4 @status:build @oracle:exact-algebra
Feature: Graded reduction crux (decided, two-sided)
  The two-sided evaluation boundary is decided by measurement plus exact certificates: the
  finite (diagonal) sector plateaus at its finite closure, while the universal sector exhibits
  exponentially many distinct operators per length (no kappa-collapse) with exact graded
  coefficient bit-size bounded linearly in word length. The W4b invariants (monotone
  distinct-kappa; per-length count exactly 2^L below L0; the diagonal positive control) hold.

  Scenario: the diagonal sector plateaus and the universal sector grows exponentially
    Given the UOR Atlas use-case
    Then the finite sector plateaus and the universal operator orbit grows exponentially with linear coefficient cost
