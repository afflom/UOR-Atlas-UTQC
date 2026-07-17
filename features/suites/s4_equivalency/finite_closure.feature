@row:finite-closure @stage:S4 @status:build @oracle:exact-algebra
Feature: Finite Braid Closure
  Scenario: the modular representation image is derived to be finite
    Given the UOR Atlas use-case
    Then the generated braiding subgroup is measured as mathematically finite
