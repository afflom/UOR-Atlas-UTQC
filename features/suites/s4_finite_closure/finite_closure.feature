@row:finite-closure @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Finite-closure representation
  Scenario: the generated subgroup is mathematically finite
    Given the UOR Atlas use-case
    Then the generated subgroup is proven mathematically finite precluding density
