@row:eps-free-decision-path @stage:S4 @status:build @oracle:exact-algebra
Feature: Epsilon-free decision path (witnessed)
  A source scan asserts the exact certifier's verdict path carries no floating-point value or
  external entropy outside explicitly delimited instrumentation spans, and that the mod-p
  projection PRNGs are deterministic fixed-seed. The determinism requirement of the
  verification architecture is thereby witnessed, not merely stated.

  Scenario: the certifier verdict path is scanned and found epsilon-free
    Given the UOR Atlas use-case
    Then the exact certifier verdict path carries no floating-point value or external entropy
