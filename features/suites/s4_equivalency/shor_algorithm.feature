@row:shor-algorithm @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Shor's Algorithm (Period Finding) Rollup
  Scenario: Shor's algorithm executes efficiently over the combinatorial manifold
    Given the UOR Atlas use-case
    Then the Shor's period finding algorithmic rollup is classically simulable due to the finite closure
