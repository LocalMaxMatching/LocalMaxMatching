#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>


#include <graphs/parallel_edge_graph.h>
#include <common/hypercube.h>


int main(int argc, char **argv)
{

	int procs;

	typedef dipl::ParallelEdgeGraph<> Graph;

	unsigned int depth;

	int dim_length = 0;

	if(argc>1)
	{
		std::stringstream io;
		io << argv[1];
		io >> procs;

		io.clear();
		io << argv[2];
		io >> dim_length;

		dipl::mem_usg_for_hypercube_matching<Graph>(procs, dim_length);
	}
	else
	{
		std::cout << "First argument number of procs, second argument dim_length" << std::endl;
	}

	return 0;
}



