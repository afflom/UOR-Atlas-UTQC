@row:shor-algorithm @stage:S4 @status:build @oracle:exact-algebra
Feature: Shor's Algorithm (exact period finding at a fixed instance)
  Scenario: the period is recovered by exact cyclotomic evaluation
    Given the UOR Atlas use-case
    Then the Shor period is recovered by exact cyclotomic evaluation and verified against the orbit
