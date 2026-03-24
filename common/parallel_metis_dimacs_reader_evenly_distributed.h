/*
 * This file can be used in parallel programs to read the necessary data
 * for each process. A single process doesn't have to store the complete
 * graph.
 *
 * This function isn't parallel (anyway I don't know yet how to read data
 * in "real" parallel. With a single disk it will always be sequential in
 * the end.)
 */

#ifndef DIPL_PARALLEL_METIS_DIMACS_READER_H_
#define DIPL_PARALLEL_METIS_DIMACS_READER_H_

#include <vector>
#include <utility>
#include <algorithm>

#ifdef _MPI
#include <mpi.h>
#endif

#include <boost/unordered_map.hpp>

#include "linear_buffered_reader.h"
#include "rating_functions.h"
#include "math_funcs.h"

namespace dipl
{

const int FILE_POS_TAG   = 13;
const int FIRST_NODE_TAG = 17;



template <typename NodeIDType, typename EdgeIDType>
void extract_header(std::string &line,
					NodeIDType &n, EdgeIDType &m,
					dipl::byte &fmt, unsigned int &ncon)
{
	std::istringstream line_stream(line);

	line_stream >> n;
	line_stream >> m;

	if(line_stream.good())
	{ // read fmt
		unsigned int tmp_fmt;
		line_stream >> tmp_fmt;

		switch(tmp_fmt)
		{// interpret the value of tmp_fmt as a binary number;
			case 0:   fmt = 0; break;
			case 1:   fmt = 1; break;
			case 10:  fmt = 2; break;
			case 11:  fmt = 3; break;
			case 100: fmt = 4; break;
			case 101: fmt = 5; break;
			case 110: fmt = 6; break;
			case 111: fmt = 7; break;
			default:  fmt = 0; break;
		}
	}
	else
	{
		fmt = 0;
	}

	if(line_stream.good())
	{ // read fmt
		line_stream >> ncon;
	}
	else
	{
		// if vertex weights are set, ncon must be at least 1 (that's the default value in this case)
		ncon = fmt & 1;
	}

}

// the src_vertex is always the first end point of an edge
// this function adds all edges, i.e. also edges to nodes with a smaller ID
// than src_vertex. That's a different behavior compared to the one of
// in metis_dimas_reader.h
template <typename WeightType, typename NodeIDType, typename Edge>
void extract_edges_of_vertex(const std::string &line, const NodeIDType src_vertex,
							 const NodeIDType ncon, const bool edge_weights,
							 const NodeIDType first_global_vertex,
							 const NodeIDType end_global_vertex,
							 const bool no_index_file,
							 WeightType (*wf)(),
							 std::vector<Edge> &edges,
							 boost::unordered_map<NodeIDType, WeightType> &node_weights,
							 NodeIDType &vertex_degree)
{
	std::istringstream line_stream(line);

	vertex_degree = 0;

	// advance to start of edges and store the weight if we only have a single node weight
	if(no_index_file)
	{ // we have to read the node weights from the graph-files
		if(ncon==1)
		{ // read the single node weight
			WeightType w;
			line_stream >> w;
			node_weights.insert(std::make_pair(src_vertex, w));
		}
		else
		{// no special treatment for multiple node weights
			node_weights.insert(std::make_pair(src_vertex, 1.));
			for(unsigned int i=0; i<ncon; i++)
			{
				unsigned int dummy;
				line_stream >> dummy;
			}
		}
	}
	else
	{
		for(unsigned int i=0; i<ncon; i++)
		{
			unsigned int dummy;
			line_stream >> dummy;
		}
	}

	while(line_stream.good())
	{
		NodeIDType tgt_vertex;
		WeightType weight = 0.;

		line_stream >> tgt_vertex;

		if(line_stream.fail())
		{ // couldn't read the last number, e.g. there wasn't one because the line
		  // ended with a whitespace, thus line_stream.good() didn't catch it
			break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
		}

		tgt_vertex--; // metis starts with 1 as the smallest vertex-ID, but we start with 0!!!

		if(edge_weights)
		{
			line_stream >> weight;

			if(line_stream.fail())
			{ // couldn't read the last number, e.g. there wasn't one because the line
			  // ended with a whitespace, thus line_stream.good() didn't catch it
				break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
			}
		}
		else
		{
			weight = wf();
		}

		// another incident edge, increase vertex degree by one
		vertex_degree++;

		// only add edges to vertices with an ID larger than first_global_vertex
		// we'll receive the other cross edges in another step
		if(tgt_vertex>=src_vertex/*first_global_vertex*/)
		{
			edges.push_back( Edge(src_vertex, tgt_vertex, weight) );
		}

	}
}


/*!
 * \brief Returns the number of incident edges of the provided vertex that are
 * added to the edge array.
 *
 * All cross edges will be count, but only local edges to higher vertices
 * will be count, otherwise we would count the local edges twice
 *
 */
template <typename NodeIDType, typename EdgeIDType>
EdgeIDType number_of_contributing_edges(const std::string &line, const NodeIDType src_vertex,
										const NodeIDType ncon, const bool edge_weights,
										const NodeIDType first_global_vertex,
										const NodeIDType end_global_vertex)
{
	EdgeIDType edge_count = 0;

	std::istringstream line_stream(line);

	// advance to start of edges
	for(unsigned int i=0; i<ncon; i++)
	{
		unsigned int dummy;
		line_stream >> dummy;
	}

	while(line_stream.good())
	{
		NodeIDType tgt_vertex;
		line_stream >> tgt_vertex;

		if(line_stream.fail())
		{ // couldn't read the last number, e.g. there wasn't one because the line
		  // ended with a whitespace, thus line_stream.good() didn't catch it
			break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
		}

		tgt_vertex--; // metis starts with 1 as the smallest vertex-ID, but we start with 0!!!

		if(edge_weights)
		{
			double weight; // dummy
			line_stream >> weight;

			if(line_stream.fail())
			{ // couldn't read the last number, e.g. there wasn't one because the line
			  // ended with a whitespace, thus line_stream.good() didn't catch it
				break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
			}
		}

		if(tgt_vertex<first_global_vertex || tgt_vertex>=end_global_vertex)
		{
			edge_count++;
		}
		else if(tgt_vertex>=src_vertex)
		{
			edge_count++;
		}

	}

	return edge_count;
}



template <typename WeightType, typename NodeIDType>
WeightType get_node_weight_from_file(const NodeIDType n,
									 std::ifstream &infile_index)
{
	// read  the start position
	std::pair<unsigned long, double> node_info;

	infile_index.seekg(sizeof(std::pair<unsigned long, double>)*n, std::ios::beg);
	infile_index.read((char*)&node_info, sizeof(std::pair<unsigned long, double>));

	return (WeightType)node_info.second;
}


template <typename NodeIDType>
unsigned long get_node_start_pos(const NodeIDType n,
								 std::ifstream &infile_index)
{
	// read  the start position
	std::pair<unsigned long, double> node_info;

	infile_index.seekg(sizeof(std::pair<unsigned long, double>)*n, std::ios::beg);
	infile_index.read((char*)&node_info, sizeof(std::pair<unsigned long, double>));

	return node_info.first;
}


template <typename WeightType, typename NodeIDType, typename EdgeIDType, class Edge, class ContainerType>
void post_process_edge_weights(ContainerType &edge_container,
							   const boost::unordered_map<NodeIDType, WeightType> &node_weights,
							   WeightType (*post_process_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight))
{

	for(EdgeIDType i=0; i<edge_container.size(); i++)
	{
		WeightType w1, w2;

		w1 = 1.;//node_weights.at(edge_container[i].n1);
		w2 = 1.;//node_weights.at(edge_container[i].n2);

		edge_container[i].weight = post_process_function(edge_container[i].weight, w1, w2);
	}
}


template <typename WeightType, typename NodeIDType, typename EdgeIDType, class Edge, class ContainerType>
void fill_node_weights_from_file(const std::string &index_filename,
								 ContainerType &edge_container,
								 boost::unordered_map<NodeIDType, WeightType> &node_weights)
{
	std::vector<NodeIDType> nodes;

	// add all nodes to node_weights
	for(EdgeIDType i=0; i<edge_container.size(); i++)
	{
		std::pair<typename boost::unordered_map<NodeIDType, WeightType>::iterator, bool> insert_result;

		insert_result = node_weights.insert(std::make_pair(edge_container[i].n1, 1.));

		if(insert_result.second)
		{// insertion took place
			nodes.push_back(edge_container[i].n1);
		}

		insert_result = node_weights.insert(std::make_pair(edge_container[i].n2, 1.));

		if(insert_result.second)
		{// insertion took place
			nodes.push_back(edge_container[i].n2);
		}
	}

	std::sort(nodes.begin(), nodes.end());

	NodeIDType lowest_node = nodes[0];
	NodeIDType highest_node = nodes[nodes.size()-1];

	NodeIDType block_size = 3000000; // read up to block_size many node-infos at once

	std::vector< std::pair<unsigned long, double> > in_buffer(block_size);

	std::ifstream infile_index;
	infile_index.open(index_filename.c_str(), std::ifstream::in | std::ifstream::binary);

	NodeIDType current_start = lowest_node;

	NodeIDType current_nodes_pos = 0;

	while(infile_index.good() && current_start<highest_node)
	{
		NodeIDType nodes_to_read;

		if( (current_start+block_size) < highest_node)
		{
			nodes_to_read = block_size;
		}
		else
		{
			nodes_to_read = highest_node - current_start;
		}


		infile_index.seekg(sizeof(std::pair<unsigned long, double>)*current_start, std::ios::beg);

		infile_index.read((char*)&in_buffer[0], sizeof(std::pair<unsigned long, double>)*nodes_to_read);


		NodeIDType current_node;

		if(current_nodes_pos<nodes.size())
		{
			current_node = nodes[current_nodes_pos];
		}
		else
		{
			current_node = current_start;
		}

		NodeIDType old_start = current_start;
		current_start += nodes_to_read;

		while(current_node < current_start)
		{
			NodeIDType pos = current_node - old_start;
			node_weights[current_node] = in_buffer[pos].second;

			current_nodes_pos++;
			if(current_nodes_pos<nodes.size())
			{
				current_node = nodes[current_nodes_pos];
			}
			else
			{
				current_node = current_start;
			}
		}
	}
}




#ifdef _MPI
template <typename NodeIDType, class Edge, class ContainerType>
void send_cross_edges_to_neighbors(const NodeIDType first_global_vertex,
								   const NodeIDType end_global_vertex,
								   std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
								   ContainerType &edge_container)
{
	MPI_Datatype MPI_EDGE_TYPE;
	MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &MPI_EDGE_TYPE);
	MPI_Type_commit(&MPI_EDGE_TYPE);

	int proc_id, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	int tag = 17;

	std::vector< std::vector<Edge> > msgs_for_proc(procs);
	std::vector<MPI_Request> reqs(procs, MPI_REQUEST_NULL);

	for(NodeIDType n=0; n<edge_container.size(); n++)
	{
		NodeIDType n1 = edge_container[n].n1;
		NodeIDType n2 = edge_container[n].n2;

		if(n1<first_global_vertex || n1>=end_global_vertex)
		{
			int partner = proc_id_of_vertex(n1, proc_ranges);
			msgs_for_proc[partner].push_back(edge_container[n]);
		}
		else if(n2<first_global_vertex || n2>=end_global_vertex)
		{
			int partner = proc_id_of_vertex(n2, proc_ranges);
			msgs_for_proc[partner].push_back(edge_container[n]);
		}
	}

	// send cross edge_container to partners
	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			if(!msgs_for_proc[p].empty())
			{// send cross edge_container to p
				MPI_Isend(&msgs_for_proc[p][0], msgs_for_proc[p].size(), MPI_EDGE_TYPE,
						  p, tag, MPI_COMM_WORLD, &reqs[p]);
			}
			else
			{// send empty message
				MPI_Isend(0, 0, MPI_EDGE_TYPE,
						  p, tag, MPI_COMM_WORLD, &reqs[p]);
			}
		}
	}


	//receive cross edge_container
	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			MPI_Status status;
			MPI_Probe(p, tag, MPI_COMM_WORLD, &status);

			int msg_count;
			MPI_Get_count(&status, MPI_EDGE_TYPE, &msg_count);

			std::vector<Edge> incoming_edges(msg_count);

			MPI_Recv(&incoming_edges[0], msg_count, MPI_EDGE_TYPE, p, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(int m=0; m<msg_count; m++)
			{
				edge_container.push_back(incoming_edges[m]);
			}
		}
	}

	MPI_Waitall(procs, &reqs[0], MPI_STATUS_IGNORE);
}
#endif


template <typename WeightType, typename NodeIDType, typename EdgeIDType, class Edge, class ContainerType>
void get_dimacs_graph(const std::string &filename,
					  NodeIDType &num_vertices,
					  EdgeIDType &num_edges,
					  ContainerType &edge_container,
					  NodeIDType &first_global_vertex,
					  NodeIDType &end_global_vertex,
					  std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
					  WeightType (*weight_func)(),
					  WeightType (*post_process_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight) = identity_edge_weight,
					  const bool verbose = false)
{
	int procs = 1;
	int proc_id = 0;

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

	if(proc_id==0)
	{
		std::cout << "Warning: parallel_metis_dimacs_reader_evenly_distributed.h doesn't support expansionstar2 rating!" << std::endl;
	}

	std::string line;

	LinearBufferedReader<1000000> infile(filename);

	unsigned int ncon;
	dipl::byte fmt;

	// find the header line
	infile.getline(line);
	while(line[0]=='%')
	{
		infile.getline(line);
	}

	extract_header(line, num_vertices, num_edges, fmt, ncon);

	EdgeIDType total_degree_per_proc = 2*num_edges/procs;

	bool is_multi_graph, contains_edge_weights, contains_vertex_weights;

	contains_vertex_weights = fmt & 1;
	contains_edge_weights   = fmt & 2;
	is_multi_graph    = fmt & 4;

	if(verbose)
	{
		std::cout << num_vertices << ", " << num_edges << ", " << (int)fmt << ", " << ncon << std::endl;
		std::cout << "contains_vertex_weights: " << contains_vertex_weights << ",  contains_edge_weights: " << contains_edge_weights << ",  is_multi_graph: " << is_multi_graph << std::endl;
	}

//	std::vector<WeightType> node_weights(num_vertices, 1);
	boost::unordered_map<NodeIDType, WeightType> node_weights;

	NodeIDType src_vertex = 0; // I numerate vertices starting with 0!!

//	EdgeIDType number_of_edges = 0;

	if(proc_id>0)
	{// first process must not receive any messages
		unsigned long pos = 0;
		NodeIDType first_node = 0;

		MPI_Recv(&pos, 1, MPI_UNSIGNED_LONG, proc_id-1, FILE_POS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&first_node, sizeof(NodeIDType), MPI_BYTE, proc_id-1, FIRST_NODE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		infile.seekg(pos);
		src_vertex = first_node;
		first_global_vertex = first_node;
	}
	else
	{
		first_global_vertex = 0;
		src_vertex = 0;
	}


	EdgeIDType degree_count = 0;
//
//	while(infile.good() && src_vertex<num_vertices &&  read_local_edges<local_edges_per_proc)
//	{
//		infile.getline(line);
//
//		if(line[0] != '%')
//		{ // not a comment-line
//			if(src_vertex>=first_global_vertex && src_vertex<end_global_vertex)
//			{
//				number_of_edges += number_of_contributing_edges<NodeIDType, EdgeIDType>(line, src_vertex,
//																ncon, contains_edge_weights,
//																first_global_vertex,
//																end_global_vertex);
//			}
//
//			src_vertex++;
//		}
//	}
//
//	edge_container.reserve(number_of_edges);
//
//
//	infile_index.open(index_filename.str().c_str(), std::ifstream::in | std::ifstream::binary);
//
//	// if there's an index-file: seek to the start position of node first_global_vertex
//	if(infile_index.good())
//	{
//		// seek to the start position of the first node of the process
//
//		// read  the start position
//		unsigned long start_pos = get_node_start_pos(first_global_vertex, infile_index);
//
//		infile.seekg(start_pos);
//		src_vertex = first_global_vertex;
//		infile_index.close();
//	}
//	else
//	{
//		src_vertex = 0;
//		infile.seekg(0);
//		infile.getline(line);
//		while(line[0]=='%')
//		{
//			infile.getline(line);
//		}
//	}

	std::vector<Edge> read_edges;

	// read edges
	while(infile.good() && src_vertex<num_vertices &&  (proc_id>=(procs-1) || degree_count<total_degree_per_proc) )
	{
		infile.getline(line);

		if(line[0] != '%')
		{ // not a comment-line
			bool no_index_file = true;
			NodeIDType vertex_degree;

			read_edges.clear();

			extract_edges_of_vertex<WeightType, NodeIDType, Edge>(line, src_vertex, ncon,
																  contains_edge_weights,
																  first_global_vertex,
																  end_global_vertex,
																  no_index_file,
																  weight_func,
																  read_edges,
																  node_weights,
																  vertex_degree);

			degree_count += vertex_degree;

			// insert edges to graph
			for(typename std::vector<Edge>::iterator e_it=read_edges.begin(); e_it!=read_edges.end(); e_it++)
			{
				if(e_it->n2>=src_vertex)
				{
					edge_container.push_back(*e_it);
				}
			}

			src_vertex++;
		}
	}

	std::cout << "termination: proc " << proc_id << ": src_vertex=" << src_vertex
			  << ", num_vertices=" << num_vertices << ", read_local_edges=" << degree_count
			  << ", local_edges_per_proc=" << total_degree_per_proc << ", procs=" << procs << std::endl;

	end_global_vertex = src_vertex;

	if(proc_id<procs-1)
	{// the last process must not send any messages
		unsigned long pos = infile.tellg();
		NodeIDType first_node = src_vertex;

		MPI_Send(&pos, 1, MPI_UNSIGNED_LONG, proc_id+1, FILE_POS_TAG, MPI_COMM_WORLD);
		MPI_Send(&first_node, sizeof(NodeIDType), MPI_BYTE, proc_id+1, FIRST_NODE_TAG, MPI_COMM_WORLD);
	}

	num_edges = edge_container.size();

	post_process_edge_weights<WeightType, NodeIDType, EdgeIDType, Edge>(edge_container, node_weights, post_process_function);


	// exchange process ranges
	proc_ranges.resize(procs);

	proc_ranges[proc_id] = std::make_pair(first_global_vertex, end_global_vertex);

#ifdef _MPI
	MPI_Allgather(&proc_ranges[proc_id], sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
				  &proc_ranges[0],       sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
				  MPI_COMM_WORLD);
#endif


#ifdef _MPI
	send_cross_edges_to_neighbors<NodeIDType, Edge, ContainerType>(first_global_vertex, end_global_vertex,
			proc_ranges, edge_container);
#endif




}


} // end of namespace dipl


#endif
