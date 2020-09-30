# Listing k-cliques

## cklist.jl

List all the $k$-cliques from an given undirected graph. For now, the output is just the number of $k$-cliques.

### Execution

> julia cklist.jl k edgelist.txt

with :

- *k* : size of the $k$-cliques to list
- *edgelist.txt* : file of the input undirected graph. The file is the list of edges, one edge in each line as two int nodes separated by a space