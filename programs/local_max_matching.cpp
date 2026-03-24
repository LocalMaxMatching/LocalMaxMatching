#ifdef MORE_INFO_LOCAL_MAX
#define MORE_INFORMATION
#endif

#include "matching.h"

#include <graphs/edge_graph.h>
#include <matching_algorithms/local_maximum_matching.h>


int main(int argc, char **argv)
{
	srand(time(NULL));

	typedef dipl::EdgeGraph<> Graph;

	int err = matching< dipl::LocalMaximumMatching<Graph>, Graph >(argc, argv);

	return err;
}










