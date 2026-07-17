@row:reconstructability @stage:S4 @status:build @oracle:exact-algebra
Feature: Topological Reconstructability from Serialized Artifacts
  Scenario: a validator reconstructs the state and kappa from the serialized genesis and word
    Given the UOR Atlas use-case
    Then any validator can perfectly mathematically reconstruct the final state and identical kappa from the genesis configuration and braid word
