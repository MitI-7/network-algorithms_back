# network-algorithms

## Maximum Flow

* ford fulkerson
    * O(nmU)
* edmonds karp
    * O(n^2 m)
* dinic
    * O(n^2 m)
* capacity scaling(dinic)
    * O(nm log U)

## Minimum Cost Flow

* cost scaling push relabel
* cycle canceling
    * O(nm^2 CU)
* out of kilter
    * O(nU * (m + n) log n)
* primal dual
    * O(min{nU, nC} n^2 m)
* successive shortest path
    * O(nU * (m + n) log n)
* primal network simplex
* dual network simplex
* parametric network simplex

## Generalized Maximum Flow
