@row:contraction-class-elimination @stage:S4 @status:build @oracle:exact-algebra
Feature: Contraction-Class Elimination (proved finite-sector polytime evaluation)
  On the framework's own finite modular sector -- the finite (collapsing) side of the
  evaluation boundary -- evaluating any braid word to its kappa is PROVED polynomial-time
  (exactly n*|W| operations, no state exponential in depth, O(n) canonicalization, a finite
  kappa-orbit), parametrically for all word lengths and all instances, by induction from the
  single machine-checked fact that each generator is a bijection of the n class-slots. The
  universal (PU(22)/PU(576)-dense) sector is excluded (its kappa-orbit is exponential); no #P
  or quantum-computational-advantage claim is attached.

  Scenario: finite-sector evaluation is proved polynomial-time with a finite kappa-orbit
    Given the UOR Atlas use-case
    Then the finite-sector contraction class evaluates in proved polynomial time with a finite kappa-orbit
