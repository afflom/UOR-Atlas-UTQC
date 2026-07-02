@row:archimedean-continuity @stage:S4 @status:some-true @oracle:mtc-axioms
Feature: Archimedean Continuity (PU(22)-Dense)
  @row:archimedean-continuity
  Scenario: the coupled generators are dense in PU(22) on the 22-dim block
    Given the UOR Atlas use-case
    Then the archimedean continuity is exactly located on the 22-dim block and saturates PU(22)
