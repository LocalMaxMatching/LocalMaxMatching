#include "parallel_matching.h"

#include <matching_algorithms/parallel_local_tree_maximum_matching.h>
#include <graphs/parallel_edge_graph.h>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int proc_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

//	if(proc_id!=0) std::cout.rdbuf(0);

	typedef dipl::ParallelEdgeGraph<> Graph;

	int err = parallel_matching< dipl::ParallelLocalTreeMaximumMatching<Graph>, Graph >(argc, argv);

	MPI_Finalize();

	return err;
}










