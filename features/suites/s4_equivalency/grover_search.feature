@row:grover-search @stage:S4 @status:build @oracle:exact-algebra
Feature: Grover's Search Algorithm (exact reference evaluation)
  Scenario: Grover's amplitude recurrence is evaluated exactly at the fixed reference instance
    Given the UOR Atlas use-case
    Then the Grover amplitude recurrence is evaluated exactly at the fixed reference instance
