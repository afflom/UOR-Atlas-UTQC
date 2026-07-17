@row:fault-tolerance @stage:S4 @status:build @oracle:exact-algebra
Feature: Deterministic Discrete Execution (replay determinism)
  Scenario: the execution is classically deterministic under replay
    Given the UOR Atlas use-case
    Then the discrete execution replays identical braid words to identical states and kappa
