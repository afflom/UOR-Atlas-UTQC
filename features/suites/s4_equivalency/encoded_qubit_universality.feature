@row:encoded-qubit-universality @stage:S4 @status:build @oracle:exact-algebra
Feature: Encoded-qubit universality (corollary of PU(24^n) density)
  On the two-handle 576-dim carrier the k-qubit code subgroup is closed in SU(576), so the
  established PU(576) density yields dense encoded single- and two-qubit gates. The encoded
  gate set is verified exactly over Q(zeta_24) and the encoding is kappa-pinned.

  Scenario: the encoded qubit gate set embeds faithfully and the density premise holds
    Given the UOR Atlas use-case
    Then the encoded-qubit gate set embeds faithfully over the certified carrier and universality follows
