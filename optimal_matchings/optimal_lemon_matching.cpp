#include <iostream>

#include <time.h>
#include <cmath>
#include <vector>
#include <list>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <common/read_command_line_arguments.h>
#include <common/get_graph.h>
#include <common/input_info.h>
#include <common/edge.h>

#include <common/hash_functions.h>

#include <lemon/smart_graph.h>
#include <lemon/matching.h>


int main(int argc, char **argv)
{
	typedef unsigned int NodeIDType;
	typedef unsigned int EdgeIDType;
	typedef double WeightType;

	typedef dipl::Edge<WeightType, NodeIDType> DiplEdge;

	typedef lemon::SmartGraph 		Graph;
	typedef Graph::EdgeMap<double>	WeightMap;

	typedef lemon::MaxWeightedMatching<Graph, WeightMap> Matching;


	bool is_root = true;

	bool read_graph, create_kn, create_grid;
	std::string path, rating;
	unsigned int num_vertices_cmd_lin, dim, dim_length;
	int seed;
	unsigned int repetitions;

	dipl::read_command_line_args(argc, argv,
								 read_graph, path,
								 create_kn, num_vertices_cmd_lin,
								 create_grid, dim, dim_length,
								 rating, seed, repetitions,
								 is_root);

	NodeIDType num_vertices = num_vertices_cmd_lin;
	EdgeIDType num_edges;
	num_edges = 0;


	std::vector<DiplEdge> *edges = new std::vector<DiplEdge>();

	boost::posix_time::ptime start_file = boost::posix_time::microsec_clock::local_time();

	dipl::get_edge_graph<WeightType, NodeIDType, EdgeIDType, DiplEdge>(read_graph, path,
																   create_kn,
																   create_grid, dim, dim_length,
																   rating,
																   seed,
																   is_root,
																   num_vertices,
																   num_edges,
																   *edges);


	Graph g;
	WeightMap wm(g);

	g.reserveNode(num_vertices);

	for(EdgeIDType n=0; n<num_vertices; n++)
	{
		g.addNode();
	}

	g.reserveEdge(edges->size());

	for(EdgeIDType n=0; n<edges->size(); n++)
	{
		Graph::Edge e;

		e = g.addEdge(g.nodeFromId((*edges)[n].n1), g.nodeFromId((*edges)[n].n2));
		wm[e] = (*edges)[n].weight;
	}

	boost::posix_time::ptime end_file = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed_time_file(end_file-start_file);

	long duration_in_micros_file = elapsed_time_file.total_microseconds();

	std::cout << "input info: " << dipl::input_info( read_graph, path, create_kn, num_vertices, create_grid, dim, dim_length, rating) << std::endl;
	std::cout << "read time: " << (double)duration_in_micros_file/1000000. << " sec" << std::endl;

	// we no longer need the edges vector
	delete edges;


	Matching m(g,wm);

//	bool maximal_matching = false;
	unsigned long matching_size = 0;
	double global_weight = 0;


	std::cout << "Nodes: " << g.nodeNum() << "   Edges: " << g.edgeNum() << std::endl;

	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();

	m.run();

	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed_time(end-start);

	long duration_in_micros = elapsed_time.total_microseconds();
	double elapsed_time_in_sec = (double)duration_in_micros/1000000.;

	matching_size = m.matchingSize();
	global_weight = m.matchingWeight();

	std::cout << "Elapsed time: " << elapsed_time_in_sec << std::endl;

	std::cout << "Size of Matching: " << matching_size << std::endl;

	std::cout << "Weight of Matching: " << global_weight << std::endl;

//	std::cout << "Is maximal Matching: " << maximal_matching << std::endl;


	return 0;
}
