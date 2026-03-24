#ifndef MATCHING_H_
#define MATCHING_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <list>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <common/read_command_line_arguments.h>
#include <common/get_graph.h>
#include <common/input_info.h>

#include <common/hash_functions.h>


template <class Matching, class Graph>
int matching(int argc, char **argv)
{
	bool is_root = true;

	typedef typename Graph::EdgeWithoutID EdgeWithoutID;
	typedef typename Graph::Edge Edge;
	typedef typename Graph::EdgeIDType EdgeIDType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::WeightType WeightType;

	unsigned int rounds = 0;

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


	std::vector<EdgeWithoutID> *edges = new std::vector<EdgeWithoutID>();

	boost::posix_time::ptime start_file = boost::posix_time::microsec_clock::local_time();

	dipl::get_edge_graph<WeightType, NodeIDType, EdgeIDType, EdgeWithoutID>(read_graph, path,
																   create_kn,
																   create_grid, dim, dim_length,
																   rating,
																   seed,
																   is_root,
																   num_vertices,
																   num_edges,
																   *edges);

	boost::posix_time::ptime end_file = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed_time_file(end_file-start_file);

	long duration_in_micros_file = elapsed_time_file.total_microseconds();

	std::cout << "input info: " << dipl::input_info( read_graph, path, create_kn, num_vertices, create_grid, dim, dim_length, rating) << std::endl;
	std::cout << "read time: " << (double)duration_in_micros_file/1000000. << " sec" << std::endl;

	std::vector<double> elapsed_times(repetitions);
	std::vector<unsigned long> cardinalities(repetitions);
	std::vector<double> weights(repetitions);
	std::vector<unsigned int> r_hash(repetitions);
	std::vector<unsigned int> all_rounds(repetitions);

	bool maximal_matching = false;
	unsigned long matching_size = 0;
	double global_weight = 0;

	for(unsigned int i=0; i<repetitions; i++)
	{
		dipl::r = rand();
		r_hash[i] = dipl::r;

		Graph g(num_vertices, *edges);

		if(i==0)
		{
			std::cout << "Nodes: " << g.num_vertices() << "   Edges: " << g.num_edges() << std::endl;
		}

		std::list<Edge> matching;

		boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();

#ifdef MORE_INFO_LOCAL_MAX
		std::list<EdgeIDType> active_edge_count;
		Matching::compute_weighted_matching(g, matching, rounds, active_edge_count);
#elif defined MORE_INFO_LOCAL_TREE
		std::list<EdgeIDType> active_edge_count;
		std::list< std::list<NodeIDType> > depths_of_trees;
		std::list< std::list<NodeIDType> > sizes_of_trees;
		Matching::compute_weighted_matching(g, matching, rounds,
											active_edge_count, depths_of_trees, sizes_of_trees);
#else
		Matching::compute_weighted_matching(g, matching, rounds);
#endif

		boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration elapsed_time(end-start);

		long duration_in_micros = elapsed_time.total_microseconds();
		double elapsed_time_in_sec = (double)duration_in_micros/1000000.;

		elapsed_times[i] = elapsed_time_in_sec;

		NodeIDType unmatched_nodes;
		maximal_matching = Matching::is_maximal_matching(matching, g, unmatched_nodes);

		matching_size = matching.size();
		cardinalities[i] = matching_size;

		global_weight = Matching::get_weight(g, matching);
		weights[i] = global_weight;

		all_rounds[i] = rounds;
	}

	if(is_root)
	{
		std::cout << "Times of rounds: ";

		for(unsigned long i=0; i<elapsed_times.size(); i++)
		{
			std::cout << elapsed_times[i] << " ";
		}

		std::cout << std::endl;


		std::cout << "Cardinalities of rounds: ";

		for(unsigned long i=0; i<elapsed_times.size(); i++)
		{
			std::cout << cardinalities[i] << " ";
		}

		std::cout << std::endl;

		std::cout << "Weights of rounds: ";

		for(unsigned long i=0; i<elapsed_times.size(); i++)
		{
			std::cout << weights[i] << " ";
		}

		std::cout << std::endl;

		std::cout << "R-Hashes of rounds: ";

		for(unsigned long i=0; i<elapsed_times.size(); i++)
		{
			std::cout << r_hash[i] << " ";
		}

		std::cout << std::endl;

		std::cout << "All rounds: ";

		for(unsigned long i=0; i<elapsed_times.size(); i++)
		{
			std::cout << all_rounds[i] << " ";
		}

		std::cout << std::endl;

		std::sort(elapsed_times.begin(), elapsed_times.end());

		unsigned long first = (double)repetitions*0.25;
		unsigned long last  = (double)repetitions*0.75;
		last += (last<repetitions);

		double elapsed_time_in_sec = 0;

		for(unsigned long i=first; i<last; i++)
		{
			elapsed_time_in_sec += elapsed_times[i];
		}

		elapsed_time_in_sec /= (double)(last-first);

		std::cout << "Elapsed time: " << elapsed_time_in_sec << std::endl;

		std::cout << "Rounds: " << rounds << std::endl;

		std::cout << "Size of Matching: " << matching_size << std::endl;

		std::cout << "Weight of Matching: " << global_weight << std::endl;

		std::cout << "Is maximal Matching: " << maximal_matching << std::endl;

#ifdef MORE_INFO_LOCAL_MAX
		active_edge_count.push_back(0);
		std::cout << "edge count development: ";
		for(typename std::list<EdgeIDType>::iterator it=active_edge_count.begin();
				it!=active_edge_count.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
#elif defined MORE_INFO_LOCAL_TREE
	    active_edge_count.push_back(0);
		std::cout << "edge count development: ";
		for(typename std::list<EdgeIDType>::iterator it=active_edge_count.begin();
				it!=active_edge_count.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;

		std::cout << "depths of trees: ";
		for(typename std::list<std::list<NodeIDType> >::iterator it=depths_of_trees.begin();
				it!=depths_of_trees.end(); it++)
		{
			for(typename std::list<NodeIDType>::iterator d_it=it->begin();
					d_it!=it->end(); d_it++)
			{
				std::cout << *d_it << " ";
			}
		}
		std::cout << std::endl;


		std::cout << "sizes of trees: ";
		for(typename std::list<std::list<NodeIDType> >::iterator it=sizes_of_trees.begin();
				it!=sizes_of_trees.end(); it++)
		{
			for(typename std::list<NodeIDType>::iterator d_it=it->begin();
					d_it!=it->end(); d_it++)
			{
				std::cout << *d_it << " ";
			}
		}
		std::cout << std::endl;
#endif

	}

	return 0;
}

#endif

