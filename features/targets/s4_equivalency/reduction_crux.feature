@row:reduction-crux @stage:S4 @status:open @oracle:exact-algebra
Feature: Graded reduction crux (measured, open)
  A measured, non-gating target. The exact graded diagonal-sector kappa-representation on the
  576-dim pair carrier is evolved under two-handle monodromy composition, and its distinct
  graded-kappa count and coefficient degree are measured against braid word length and
  reported via `just report`. It is never asserted as a reduction in either direction; a
  plateau is evidence of finite closure on the diagonal sector. This feature is non-gating
  and is not executed by the BDD suite runner.

  Scenario: the graded diagonal-sector kappa growth is measured and reported
    Given the UOR Atlas use-case
    Then the graded diagonal-sector kappa growth is measured and reported without any reduction assertion
