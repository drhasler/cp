# Graph theory

## Requirements
you should be familiar with the terminology:
node, edge, graph, tree, cycle, directed

## Content
### components
disjoint set union, topological sort, articulation points, bridges, biconnected components, connected components, lowest common ancestor
### spanning trees
prim and kruskal MST, Kirchhoff's thm

## doesnt contain
- DP on trees (will be in dp folder)
- traversal (will be added)

## tips n triks
- rule of thumb: if there is a problem about a graph which is not
"traversal specific" (doesnt mean much), you should try
reducing it to a problem on a tree (often MST) / forest
to process information in a well defined order.
- for many tree problems, you should DFS post order and
merge information about "lower" subtrees
- BFS only to find short paths
