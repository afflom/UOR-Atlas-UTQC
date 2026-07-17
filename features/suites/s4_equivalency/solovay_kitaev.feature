@row:solovay-kitaev @stage:S4 @status:build @oracle:exact-algebra
Feature: Solovay-Kitaev Density (Decided)
  @row:solovay-kitaev
  Scenario: the density question is exactly decided over Q(zeta_24)
    Given the UOR Atlas use-case
    Then the Solovay-Kitaev density question is exactly decided
