@row:modular-s-t @stage:S2 @status:build @oracle:mtc-axioms
Feature: Modular S/T matrices
  Constructed as the Atlas-native pointed category C(Z_modality x Z_2^k, q) and validated
  phase-exactly against the full MTC axiom set: S symmetric and unitary, S^4 = I, S^2 = C
  (charge conjugation), twists of finite order, (ST)^3 = p+ S^2 with the Gauss-sum anomaly
  p+, non-negative integer fusion, exact Verlinde, phase-exact pentagon/hexagon/balancing,
  and the monodromy-S relation. A build construction validated against the axioms.

  Scenario: S and T satisfy the SL(2,Z) relations for the Atlas
    Given the UOR Atlas use-case
    Then the modular S and T satisfy the SL(2,Z) relations

  Scenario Outline: the modular data holds for arbitrary use-cases
    Given an arbitrary use-case with scope <q> modality <T> context <O>
    Then the modular S and T satisfy the SL(2,Z) relations

    Examples:
      | q | T | O |
      | 4 | 3 | 8 |
      | 2 | 2 | 4 |
