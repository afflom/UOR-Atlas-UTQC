@row:tensor-contraction-bypass @stage:S4 @status:build @oracle:holospaces-cc
Feature: Isotopy kappa-collision (decision by equivalence)
  Topological decision problems of the form "do these braid words realize the same operator"
  are answered by content-addressed collision of their kappa forms. No #P-hardness claim is
  attached to this mechanism.

  Scenario: isotopic braid operations collide on identical kappa forms
    Given the UOR Atlas use-case
    Then isomorphic topological braid operations naturally collide on identical kappa forms via classical equivalence
