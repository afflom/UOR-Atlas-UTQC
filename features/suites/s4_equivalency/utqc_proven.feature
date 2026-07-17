@row:utqc-proven @stage:S4 @status:build @oracle:exact-algebra
Feature: UTQC roll-up conjunction
  The roll-up invokes every pillar witness explicitly and fails when any pillar fails.

  Scenario: Roll-up
    Given the UOR Atlas use-case
    Then every UTQC pillar witness passes in the roll-up conjunction
