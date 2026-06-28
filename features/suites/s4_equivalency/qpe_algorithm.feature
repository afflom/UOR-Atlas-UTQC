@row:qpe-algorithm @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Quantum Phase Estimation Rollup
  Scenario: QPE executes efficiently over the combinatorial manifold
    Given the UOR Atlas use-case
    Then the QPE algorithmic rollup is classically simulable due to the finite closure
