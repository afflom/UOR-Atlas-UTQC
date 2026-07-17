@row:generative-closure @stage:S4 @status:build @oracle:exact-algebra
Feature: Generative closure (orbit decomposition)
  The generator group closure is computed by BFS: the class space partitions into the
  derived mirror orbits, covering every label from one seed per mirror class.

  Scenario: Orbit decomposition
    Given the UOR Atlas use-case
    Then the generator closure partitions the class space into the derived mirror orbits covering every label
