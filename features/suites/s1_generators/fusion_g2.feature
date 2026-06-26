@row:fusion-g2 @stage:S1 @status:some-true @oracle:uor-addr-composition
Feature: Atlas composition operation g2
  The Atlas composition operation g2 is the realized CS-G2 commutative product on every sigma-axis (sha256, blake3,
  sha3-256, keccak256, sha512), and the composition norm is multiplicative at the use-case's
  context level (the octonion eight-square at O=8).

  Scenario: fusion is commutative and norm-multiplicative for the Atlas
    Given the UOR Atlas use-case
    Then fusion reduces to the realized g2 product on every sigma-axis

  Scenario Outline: fusion holds for arbitrary use-cases
    Given an arbitrary use-case with scope <q> modality <T> context <O>
    Then fusion reduces to the realized g2 product on every sigma-axis

    Examples:
      | q | T | O |
      | 4 | 3 | 8 |
      | 2 | 2 | 4 |
