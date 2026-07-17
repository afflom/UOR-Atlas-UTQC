@row:s4-modal-logic @stage:S4 @status:build @oracle:exact-algebra
Feature: S4 Modal Frame (reflexivity and transitivity, verified by enumeration)
  Scenario: the generator-reachability frame satisfies the S4 axioms
    Given the UOR Atlas use-case
    Then the S4 modal logic frame satisfies reflexivity and transitivity
