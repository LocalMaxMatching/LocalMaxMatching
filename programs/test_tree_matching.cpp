#include <matching_algorithms/tree_matching.h>
#include <graphs/edge_graph.h>

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <list>



int main(int argc, char **argv)
{
	typedef dipl::EdgeGraph<> Graph;
	typedef Graph::Edge Edge;
	typedef Graph::EdgeIDType EdgeIDType;

	unsigned int node_count = 31;

	std::list<Edge> edges;

//	get_random_graph(node_count, node_count*sqrt(node_count), max_weight, edges);

	edges.push_back(Edge(0, 1, 1));
	edges.push_back(Edge(1, 2, 1));
	edges.push_back(Edge(1, 3, 1));
	edges.push_back(Edge(1, 4, 1));
	edges.push_back(Edge(2, 5, 1));
	edges.push_back(Edge(2, 10, 1));
	edges.push_back(Edge(2, 11, 1));
	edges.push_back(Edge(2, 12, 1));
	edges.push_back(Edge(3, 6, 1));
	edges.push_back(Edge(3, 7, 1));
	edges.push_back(Edge(3, 8, 1));
	edges.push_back(Edge(4, 9, 1));
	edges.push_back(Edge(6, 13, 1));
	edges.push_back(Edge(7, 14, 1));
	edges.push_back(Edge(7, 15, 1));
	edges.push_back(Edge(8, 16, 1));
	edges.push_back(Edge(8, 17, 1));
	edges.push_back(Edge(9, 18, 1));
	edges.push_back(Edge(10, 19, 1));
	edges.push_back(Edge(10, 20, 1));
	edges.push_back(Edge(17, 21, 1));
	edges.push_back(Edge(18, 23, 1));
	edges.push_back(Edge(18, 24, 1));
	edges.push_back(Edge(20, 22, 1));

	edges.push_back(Edge(25, 27, 1));
	edges.push_back(Edge(26, 27, 1));
	edges.push_back(Edge(27, 28, 1));
	edges.push_back(Edge(28, 29, 1));
	edges.push_back(Edge(28, 30, 1));

	std::cout << "Vertices of Graph: " << node_count << ",  Edges of Graph: " << edges.size() << std::endl;

	Graph g(node_count, edges);

	std::list<EdgeIDType> tree_edges;

	for(Graph::edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
	{
		tree_edges.push_back(g.get_edge_id(e_it));
	}

	std::cout << "correct  initial state: " << g.correct_state() << std::endl;

	std::list<EdgeIDType> matching;

	max_weighted_matching_of_edge_list(g, tree_edges, matching);

	for(std::list<EdgeIDType>::iterator it=matching.begin(); it!=matching.end(); it++)
	{
		std::cout << g.get_edge(*it) << std::endl;
	}

	std::cout << "Size of Matching: " << matching.size() << std::endl;

#ifdef MORE_INFORMATION
	active_edge_count.push_back(0);

	for(std::list<Graph::EdgeIDType>::iterator it=active_edge_count.begin(); it!=active_edge_count.end(); it++)
	{
		Graph::EdgeIDType current = *it;
		std::cout << current << "\t";
	}
	std::cout << std::endl;
#endif


	return 0;
}
