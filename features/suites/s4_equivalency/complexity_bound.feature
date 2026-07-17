@row:complexity-bound @stage:S4 @status:build @oracle:exact-algebra
Feature: Linear Execution Cost Bound (exact operation count)
  Scenario: execution cost is counted exactly and is linear in braid depth
    Given the UOR Atlas use-case
    Then the execution cost is exactly the operation count linear in braid depth with no exponential state
