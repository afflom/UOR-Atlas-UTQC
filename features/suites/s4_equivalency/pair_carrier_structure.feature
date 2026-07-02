@row:pair-carrier-structure @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Pair-Carrier Density (PU(576))
  @row:pair-carrier-structure
  Scenario: the two-handle carrier structure is exactly decided
    Given the UOR Atlas use-case
    Then the pair carrier is irreducible and its closure is dense in PU(576)
