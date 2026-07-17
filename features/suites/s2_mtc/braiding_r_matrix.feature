@row:braiding-r-matrix @stage:S2 @status:build @oracle:mtc-axioms
Feature: Braiding R-matrix
  The braiding of the Atlas-native pointed category, decided with exact integer arithmetic
  on root-of-unity exponents: the monodromy bicharacter is symmetric, bimultiplicative,
  twist-consistent, and non-degenerate (the modularity criterion for pointed categories),
  tied back to the runtime R-symbols, then validated through the phase-exact hexagon,
  balancing, and monodromy-S axiom checks. A build construction validated against the axioms.

  Scenario: the R-matrix satisfies the hexagon and Yang-Baxter for the Atlas
    Given the UOR Atlas use-case
    Then the braiding R satisfies the hexagon and Yang-Baxter

  Scenario Outline: the braiding holds for arbitrary use-cases
    Given an arbitrary use-case with scope <q> modality <T> context <O>
    Then the braiding R satisfies the hexagon and Yang-Baxter

    Examples:
      | q | T | O |
      | 4 | 3 | 8 |
      | 2 | 2 | 4 |
