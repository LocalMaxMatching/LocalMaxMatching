# LocalMaxMatching

Efficient parallel and sequential implementations of the **Local Maximum Matching** algorithm for computing maximal and approximate maximum weight matchings in graphs.

The code was authored by Marcel Birn.

## Overview

A matching M of a graph G = (V, E) is a subset of edges such that no two elements of M share an endpoint. This repository implements the *local max* algorithm: an edge is called *locally maximal* if its weight is larger than the weight of any incident edge (with tie-breaking via random perturbations). The algorithm repeatedly adds locally maximal edges to the matching and removes their incident edges until no edges remain. The result is a maximal matching that provides a **1/2-approximation** of the maximum weight matching.

Key properties:
- **Linear expected work** for computing maximal matchings
- **O(log^2 n) expected time** on the CREW PRAM model
- Adapts to CRCW PRAM, MapReduce, external memory, and distributed memory (MPI) models
- In practice, approximately 75% of edges are removed per iteration
- Solution quality is close to the greedy algorithm and within a few percent of optimal

## Paper

If you use this software, please cite:

> Marcel Birn, Vitaly Osipov, Peter Sanders, Christian Schulz, Nodari Sitchinava. **Efficient Parallel and External Matching.** Euro-Par 2013. Lecture Notes in Computer Science, vol 8097. Springer.

```bibtex
@inproceedings{DBLP:conf/europar/BirnOSSS13,
  author    = {Marcel Birn and Vitaly Osipov and Peter Sanders and Christian Schulz and Nodari Sitchinava},
  title     = {Efficient Parallel and External Matching},
  booktitle = {Euro-Par 2013 Parallel Processing},
  series    = {Lecture Notes in Computer Science},
  volume    = {8097},
  pages     = {659--670},
  publisher = {Springer},
  year      = {2013},
  doi       = {10.1007/978-3-642-40047-6_66}
}
```

ArXiv preprint: [https://arxiv.org/abs/1302.4587](https://arxiv.org/abs/1302.4587)

## Requirements

- C++ compiler with C++11 support (g++, clang++, or icpc)
- [Boost](https://www.boost.org/) (program_options library)
- MPI (OpenMPI or MPICH) for parallel programs
- Optional: [LEMON](https://lemon.cs.elte.hu/) graph library for computing optimal matchings

## Building

### CMake (recommended)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Binaries are placed in `build/`. MPI targets are built automatically if MPI is found; otherwise they are skipped.

### Make (legacy)

```bash
# Edit the makefile to set OPTIMIZATION = -O3, then:
make bin/local_max_matching

# Parallel variants (requires MPI)
make bin/parallel_local_max_matching
```

You can set `BOOST_INC` and `BOOST_LIB` if Boost is installed in a non-standard location:

```bash
make bin/local_max_matching BOOST_INC=/path/to/boost/include BOOST_LIB=/path/to/boost/lib
```

### Build Targets

**Sequential matching programs:**

| Target | Description |
|--------|-------------|
| `local_max_matching` | Local maximum matching |
| `local_tree_matching` | Local tree maximum matching (tree DP refinement) |
| `karp_sipser_matching` | Karp-Sipser approximate matching |
| `mixed-karp_sipser-local_max` | Hybrid Karp-Sipser and local max |

**Parallel matching programs (MPI):**

| Target | Description |
|--------|-------------|
| `parallel_local_max_matching` | Distributed local maximum matching |
| `parallel_local_max_matching_more_info` | With communication statistics |

**Utility programs:**

| Target | Description |
|--------|-------------|
| `add_random_weights_to_metis_graph` | Add random edge weights to a METIS graph |
| `create_node_index_for_metis_graph` | Create byte-offset index for fast graph access |
| `arrange_nodes_along_a_curve` | Reorder nodes spatially (improves locality) |
| `matrix_market_to_bipartite_metis` | Convert Matrix Market format to bipartite METIS |
| `compute_edge_distribution` | Analyze edge distribution across MPI processes |

## Usage

### Sequential

```bash
# Compute a maximum weight matching on a graph in METIS format
./bin/local_max_matching --read_graph graph.metis --edge_rating weight

# Compute a maximal cardinality matching (unit weights)
./bin/local_max_matching --read_graph graph.metis --edge_rating const

# Use random edge weights with a specific seed
./bin/local_max_matching --read_graph graph.metis --edge_rating rand --seed 42

# Run on a complete graph K_100
./bin/local_max_matching --kn --vertices 100 --edge_rating rand --seed 1

# Run on a 2D grid graph
./bin/local_max_matching --grid --dim 2 --dim_length 100 --edge_rating rand --seed 1

# Run multiple repetitions
./bin/local_max_matching --read_graph graph.metis --edge_rating weight --repetitions 10
```

### Parallel (MPI)

```bash
# Run with 4 MPI processes
mpirun -np 4 ./bin/parallel_local_max_matching --read_graph graph.metis --edge_rating weight

# With extended communication statistics
mpirun -np 16 ./bin/parallel_local_max_matching_more_info --read_graph graph.metis --edge_rating weight
```

### Command-Line Options

| Option | Description |
|--------|-------------|
| `--help, -h` | Show help message |
| `--read_graph PATH` | Read graph from METIS/DIMACS format file |
| `--kn` | Generate a complete graph K_n |
| `--vertices N` | Number of vertices (for `--kn`) |
| `--grid` | Generate a grid graph |
| `--dim D` | Grid dimension (for `--grid`) |
| `--dim_length L` | Grid side length (for `--grid`) |
| `--edge_rating TYPE` | Edge rating function (see below) |
| `--seed SEED` | Random seed (for `rand` edge rating) |
| `--repetitions N` | Number of repetitions (default: 1) |

**Edge rating functions:**

| Rating | Description |
|--------|-------------|
| `const` | Unit weights (computes maximal cardinality matching) |
| `weight` | Use edge weights from the graph file |
| `expansionstar2` | Weight function based on node/edge properties |
| `rand` | Random weights (requires `--seed`) |

### Output Format

```
input info: <input description>
read time: 0.12345 sec
Nodes: 1000   Edges: 5000
Times of rounds: 0.001 0.0005 0.0002
Cardinalities of rounds: 400 80 15
Weights of rounds: 150.5 30.2 5.1
Elapsed time: 0.00170
Rounds: 3
Size of Matching: 495
Weight of Matching: 185.8
Is maximal Matching: true
```

## Graph Input Format

The primary input format is **METIS/DIMACS**:

```
<num_vertices> <num_edges> [<format>]
<adjacency list for vertex 1>
<adjacency list for vertex 2>
...
```

The optional `format` field is a bitmask: bit 0 = vertex weights, bit 1 = edge weights. For weighted graphs (format `10` in binary = `2` decimal or `11` in binary = `3` for both vertex and edge weights), each adjacency list entry is `neighbor weight` pairs. Vertices are 1-indexed in the file.

Example (4 vertices, 4 edges, with edge weights):

```
4 4 1
2 1.0 3 2.0
1 1.0 4 3.0
1 2.0 4 1.5
2 3.0 3 1.5
```

The utility `matrix_market_to_bipartite_metis` converts Matrix Market sparse matrices to bipartite METIS format.

## Algorithms

### Local Maximum Matching

Each iteration performs three passes over the remaining edges:

1. **Candidate selection:** For each node, find the incident edge with the highest weight
2. **Matching:** If both endpoints of an edge mutually select each other, add the edge to the matching
3. **Cleanup:** Remove all edges incident to matched nodes

In expectation, each iteration removes at least 50% of the remaining edges (approximately 75% in practice), yielding O(log n) iterations and O(m) total work.

### Local Tree Maximum Matching

Extends local maximum matching by building spanning trees and applying dynamic programming to find optimal matchings on these trees. This can improve matching quality at the cost of additional computation.

### Karp-Sipser

Prioritizes degree-1 vertices (which can always be optimally matched) and falls back to random selection otherwise.

### Mixed Karp-Sipser / Local Max

Combines Karp-Sipser's degree-1 vertex prioritization with the local maximum approach for remaining vertices.

## Project Structure

```
LocalMaxMatching/
  common/                 # I/O, graph generation, utilities
  graphs/                 # Graph data structures (sequential and parallel)
  matching_algorithms/    # Algorithm implementations (header-only)
  programs/               # Main programs (.cpp files)
  optimal_matchings/      # Optimal matching via LEMON (optional)
  bin/                    # Build output directory
  makefile                # Build system
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
