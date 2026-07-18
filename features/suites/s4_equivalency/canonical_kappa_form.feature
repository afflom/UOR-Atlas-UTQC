@row:canonical-kappa-form @stage:S4 @status:build @oracle:exact-algebra
Feature: Canonical-form kappa witnesses (crux prerequisite)
  Before any distinct-kappa count is interpretable, the canonical serialization must roundtrip
  byte-identically, one operator reached by two factorizations must land on the identical
  kappa, and the harness kappa producer must be the substrate content addresser.

  Scenario: the canonical-form kappa invariants hold
    Given the UOR Atlas use-case
    Then the canonical-form kappa serialization roundtrips and two factorization paths agree
