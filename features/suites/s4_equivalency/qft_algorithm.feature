@row:qft-algorithm @stage:S4 @status:build @oracle:exact-algebra
Feature: Quantum Fourier Transform (scheduling and class-space execution)
  Scenario: the QFT circuit schedules onto a bounded braid word
    Given the UOR Atlas use-case
    Then the QFT circuit schedules onto a bounded braid word and executes on the class space
