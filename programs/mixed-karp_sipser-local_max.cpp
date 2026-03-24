#include "matching.h"

#include <graphs/adjacency_array_graph_with_active_edges_vector.h>
#include <matching_algorithms/mixed-karp_sipser_and_local_max.h>


int main(int argc, char **argv)
{

	typedef dipl::AdjacencyArrayGraphWithActiveEdgesVector<double, unsigned int, unsigned int> Graph;

	int err = matching< dipl::MixedKarpSipserLocalMax<Graph>, Graph >(argc, argv);

	return err;
}
