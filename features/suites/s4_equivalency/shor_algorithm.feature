@row:shor-algorithm @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Shor's Algorithm (Period Finding) Rollup
  Scenario: Shor's algorithm executes efficiently over the combinatorial manifold
    Given the UOR Atlas use-case
    Then a complex algorithmic rollup executes Shor's period finding with polynomial braid compilation and fully evaluates bypassing tensor contraction
