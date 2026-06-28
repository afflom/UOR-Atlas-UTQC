@row:qft-algorithm @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Quantum Fourier Transform Rollup
  Scenario: QFT executes efficiently over the combinatorial manifold
    Given the UOR Atlas use-case
    Then the QFT algorithmic rollup is classically simulable due to the finite closure
