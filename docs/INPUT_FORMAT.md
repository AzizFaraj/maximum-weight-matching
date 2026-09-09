# Input format and normalization

The benchmark driver reads Matrix Market coordinate files. Its normalization
is part of the experimental identity and must not be changed silently.

## General graphs

- The declared matrix must be square.
- Matrix indices are converted from one-based to zero-based vertex IDs.
- Each off-diagonal entry denotes one undirected edge.
- Self-loops are discarded.
- Endpoint order is canonicalized.
- Duplicate endpoint pairs are merged by retaining their maximum weight.
- Isolated vertices declared by the matrix dimensions are preserved.

If a symmetric file stores both `(i,j)` and `(j,i)`, those records therefore
become one undirected edge.

## Bipartite matrices

Bipartite inputs are supported by Edmonds and Gabow. Rows and columns are
different vertex sets:

```text
number of vertices = rows + columns
edge(i,j) = row vertex i -- column vertex (rows + j)
```

Thus equal row and column indices do not identify one vertex. This disjoint
namespace rule avoids the historical conversion error that produced different
graphs by using only `max(rows, columns)` vertices.

Scaling currently accepts only general-graph inputs through this driver.

## Weights

In `file` mode, the third Matrix Market field is used as the edge weight; a
missing field defaults to one. The benchmark contract requires finite,
nonnegative integral weights representable by the implementation. In
`uniform` mode, weights in `[1,W]` are generated deterministically from the
recorded seed.

Every run reports a digest over the normalized vertex count, edge set, and
exact integer weights. Reproduction scripts must verify this digest before
comparing solver outputs.
