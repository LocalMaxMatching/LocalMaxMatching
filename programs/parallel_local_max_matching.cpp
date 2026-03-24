#include "parallel_matching.h"

#include <graphs/parallel_edge_graph.h>
#include <matching_algorithms/parallel_local_maximum_matching.h>


int main(int argc, char **argv)
{
//	{
//		int i = 1;
//		char hostname[256];
//		gethostname(hostname, sizeof(hostname));
//		printf("PID %d on %s ready for attach\n", getpid(), hostname);
//		fflush(stdout);
//		while (0 == i)
//			sleep(5);
//	}

	MPI_Init(&argc, &argv);

	int proc_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

//	if(proc_id!=0) std::cout.rdbuf(0);

	typedef dipl::ParallelEdgeGraph<> Graph;

	int err = parallel_matching< dipl::ParallelLocalMaximumMatching<Graph>, Graph >(argc, argv);

	MPI_Finalize();

	return err;
}










