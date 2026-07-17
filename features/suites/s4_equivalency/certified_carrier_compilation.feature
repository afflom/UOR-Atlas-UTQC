@row:certified-carrier-compilation @stage:S4 @status:build @oracle:exact-algebra
Feature: Certified-carrier compilation (unconditional rotations)
  On the density-certified 22-dim carrier the compiler synthesizes arbitrary rotations by a
  deterministic bounded spectral-flow search, so Shor's continuous-phase QFT rotations
  compile unconditionally; the compiled Shor instance replays against the exact reference.

  Scenario: Shor's QFT rotations compile on the certified carrier and the instance replays
    Given the UOR Atlas use-case
    Then the certified carrier compiles Shor's rotations and the instance replays against the exact reference
