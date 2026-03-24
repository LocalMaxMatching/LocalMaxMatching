#ifndef PARALLEL_MATCHING_H_
#define PARALLEL_MATCHING_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <list>

#include <fstream>

#include <mpi.h>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <common/read_command_line_arguments.h>
#include <common/get_graph.h>
#include <common/input_info.h>
#include <common/mpi_system_support.h>


template <class Matching, class Graph>
int parallel_matching(int argc, char **argv)
{
	int proc_id, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// warm up the network - necessary for the IC1
	dipl::warm_up_network(MPI_COMM_WORLD);


	bool is_root = proc_id==0;

	typedef typename Graph::Edge Edge;
	typedef typename Graph::EdgeIDType EdgeIDType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::WeightType WeightType;

	unsigned int depth;

	bool read_graph, create_kn, create_grid;
	std::string path, rating;
	std::string matching_output_file;
	unsigned int num_vertices_cmd_lin, dim, dim_length;
	int seed;
	unsigned int repetitions;

	dipl::read_command_line_args(argc, argv,
								 read_graph, path,
								 create_kn, num_vertices_cmd_lin,
								 create_grid, dim, dim_length,
								 rating, seed, repetitions,
								 is_root,
								 &matching_output_file);

	NodeIDType num_vertices = num_vertices_cmd_lin;
	EdgeIDType num_edges;
	num_edges = 0;

	NodeIDType first_global_vertex;
	NodeIDType end_global_vertex;

	std::vector<std::pair<NodeIDType, NodeIDType> > proc_ranges;

//	std::vector<Edge> *edges = new std::vector<Edge>();

	boost::posix_time::ptime start_file = boost::posix_time::microsec_clock::local_time();

	Graph g;

	dipl::get_edge_graph<WeightType, NodeIDType, EdgeIDType, Edge>(read_graph, path,
																   create_kn,
																   create_grid, dim, dim_length,
																   rating,
																   seed,
																   is_root,
																   num_vertices,
																   num_edges,
																   g,
																   first_global_vertex, end_global_vertex,
																   proc_ranges);

	MPI_Barrier(MPI_COMM_WORLD);

	boost::posix_time::ptime end_file = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed_time_file(end_file-start_file);

	long duration_in_micros_file = elapsed_time_file.total_microseconds();

	if(is_root)
	{
		std::cout << "input info: " << dipl::input_info( read_graph, path, create_kn, num_vertices, create_grid, dim, dim_length, rating) << std::endl;
		std::cout << "using " << procs << " processes" << std::endl;
		std::cout << "read time: " << (double)duration_in_micros_file/1000000. << " sec" << std::endl;
	}

	g.initialize(first_global_vertex, end_global_vertex, proc_ranges);

//	for(int p=0; p<procs; p++)
//	{
//		if(proc_id==p)
//		{
//			std::cout << "proc " << p << ": " << std::endl;
//			for(typename Graph::edge_iterator e_it=g.begin_local_edges(); e_it<g.end_local_edges(); e_it++)
//			{
//				std::cout << "(" << g.local_vertex_id_to_global_id(e_it->n1)
//						  << ", " << g.local_vertex_id_to_global_id(e_it->n2)
//						  << ",, " << e_it->weight << std::endl;
//			}
//
//			for(typename Graph::edge_iterator e_it=g.begin_local_cross_edges(); e_it<g.end_local_cross_edges(); e_it++)
//			{
//				std::cout << "(" << g.local_vertex_id_to_global_id(e_it->n1)
//						  << ", " << g.local_vertex_id_to_global_id(e_it->n2)
//						  << ",, " << e_it->weight << std::endl;
//			}
//		}
//		MPI_Barrier(MPI_COMM_WORLD);
//	}

//	bool delete_edge_vector=true;

//	Graph g(num_vertices, *edges, first_global_vertex, end_global_vertex, proc_ranges, delete_edge_vector);


//	delete edges; // we no longer need those edges

	unsigned long global_local_edge_count, global_cross_edge_count;
	unsigned long local_edge_count = g._num_local_edges;
	unsigned long local_cross_edge_count = g._num_local_cross_edges;

	MPI_Reduce(&local_edge_count, &global_local_edge_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_cross_edge_count, &global_cross_edge_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if(is_root)
	{
		std::cout << "vertices: " << num_vertices << "  local_edges: " << global_local_edge_count << "  cross edges: " << global_cross_edge_count << std::endl;
	}


	std::vector<Edge> matching;

	// synchronize before computing the matching, other process might still read the input graphs
	MPI_Barrier(MPI_COMM_WORLD);
	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();

#ifdef MORE_INFORMATION
	Matching::communication_factor = 1;
#endif

	Matching::compute_weighted_matching(g, matching, depth);

//	you might want to use the following function, when using 2D-Grids
//	dipl::ParallelLocalTreeMaximumMatching<Graph, dipl::Grid2D::gvid_to_nvid>::compute_weighted_matching(g, matching, depth);

	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed_time(end-start);

	long duration_in_micros = elapsed_time.total_microseconds();

	// get the maximal running time of all processes
	long max_duration_in_micros;
	MPI_Reduce(&duration_in_micros, &max_duration_in_micros, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

	double max_time_in_sec = (double)max_duration_in_micros/1000000.;

	bool maximal_matching = g.is_maximal_matching(matching);

	unsigned int global_matching_size;
	unsigned int local_matching_size = matching.size();

	MPI_Reduce(&local_matching_size, &global_matching_size, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	unsigned int global_depth;
	MPI_Reduce(&depth, &global_depth, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);

	double global_weight = Matching::get_weight(g, matching);

#ifdef MORE_INFORMATION

	std::stringstream filename;
	filename << "more_info_" << procs << "_" << proc_id;
	std::ofstream outfile(filename.str().c_str());

	outfile << "proc " << proc_id << ": " << Matching::sent_bytes << " Byte,  "
				  << Matching::number_of_sent_messages << " msgs,  "
				  << Matching::communication_time << " sec"
				  << ",  " << local_edge_count
				  << ",  " << local_cross_edge_count
				  << ",  comm_factor = " << Matching::communication_factor << std::endl;
#endif

	if(is_root)
	{
		std::cout << "Elapsed time: " << max_time_in_sec << std::endl;

		std::cout << "Rounds: " << global_depth << std::endl;

		std::cout << "Size of Matching: " << global_matching_size << std::endl;

		std::cout << "Weight of Matching: " << global_weight << std::endl;

		std::cout << "Is maximal Matching: " << maximal_matching << std::endl;
	}

	if(!matching_output_file.empty())
	{
		// matching edges already contain global vertex IDs
		std::vector<unsigned int> local_pairs;
		for(typename std::vector<Edge>::iterator it=matching.begin();
				it!=matching.end(); it++)
		{
			local_pairs.push_back(it->n1);
			local_pairs.push_back(it->n2);
		}

		int local_count = local_pairs.size();
		std::vector<int> recv_counts(procs);
		MPI_Gather(&local_count, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

		std::vector<int> displs(procs, 0);
		int total_count = 0;
		if(is_root)
		{
			for(int i=0; i<procs; i++)
			{
				displs[i] = total_count;
				total_count += recv_counts[i];
			}
		}

		std::vector<unsigned int> all_pairs(total_count);
		MPI_Gatherv(&local_pairs[0], local_count, MPI_UNSIGNED,
					&all_pairs[0], &recv_counts[0], &displs[0], MPI_UNSIGNED,
					0, MPI_COMM_WORLD);

		if(is_root)
		{
			std::ofstream outfile(matching_output_file.c_str());
			if(outfile.is_open())
			{
				for(int i=0; i<total_count; i+=2)
				{
					outfile << all_pairs[i] << " " << all_pairs[i+1] << std::endl;
				}
				outfile.close();
				std::cout << "Matching written to " << matching_output_file << std::endl;
			}
			else
			{
				std::cerr << "Error: could not open " << matching_output_file << " for writing." << std::endl;
			}
		}
	}

	return 0;
}

#endif

