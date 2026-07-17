@row:qpe-algorithm @stage:S4 @status:build @oracle:exact-algebra
Feature: Quantum Phase Estimation (exact readout)
  Scenario: the QPE readout is computed by exact integer minimization
    Given the UOR Atlas use-case
    Then the QPE readout peak is computed by exact integer minimization and meets the exact guarantee
