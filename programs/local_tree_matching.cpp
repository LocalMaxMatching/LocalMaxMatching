#ifdef MORE_INFO_LOCAL_TREE
#define MORE_INFORMATION
#endif

#include "matching.h"

#include <graphs/edge_graph.h>
#include <matching_algorithms/local_tree_maximum_matching.h>

int main(int argc, char **argv)
{

	typedef dipl::EdgeGraph<> Graph;

	int err = matching< dipl::LocalTreeMaximumMatching<Graph>, Graph >(argc, argv);

	return err;
}
