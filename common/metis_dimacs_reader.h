#ifndef _DIPL_METIS_DIMACS_READER_H_
#define _DIPL_METIS_DIMACS_READER_H_

#include <vector>

#include "linear_buffered_reader.h"
#include "rating_functions.h"

namespace dipl
{



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


template <typename WeightType, typename NodeIDType, typename Edge>
void extract_edges_of_vertex(const std::string &line, const NodeIDType src_vertex,
							 const NodeIDType ncon, const bool edge_weights,
							 WeightType (*wf)(),
							 std::vector<Edge> &edges,
							 std::vector<WeightType> &node_weights)
{
	std::istringstream line_stream(line);

	// advance to start of edges, we're not interested in vertex-weights
	if(ncon==1)
	{ // read the single node weight
		line_stream >> node_weights[src_vertex];
	}
	else
	{// no specially treatment for multiple node weights
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

		tgt_vertex--; // metis starts with 1 as the smallest vertex-ID, but I start with 0!!!

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

		// we already added edges to vertices with a smaller vertex-ID
		if(tgt_vertex >= src_vertex)
		{
			edges.push_back( Edge(src_vertex, tgt_vertex, weight) );
		}

	}
}


template <typename WeightType, typename EdgeIDType, class Edge>
void post_process_edge_weights(std::vector<Edge> &edges,
							   const std::vector<WeightType> &node_weights,
							   WeightType (*post_process_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight))
{
	for(EdgeIDType i=0; i<edges.size(); i++)
	{
		edges[i].weight = post_process_function(edges[i].weight, node_weights[edges[i].n1], node_weights[edges[i].n2]);
	}
}


template <typename WeightType, typename NodeIDType, typename EdgeIDType, class Edge>
void get_dimacs_graph(const std::string &filename,
					  NodeIDType &num_vertices,
					  EdgeIDType &num_edges,
					  std::vector<Edge> &edges,
					  WeightType (*weight_func)(),
					  WeightType (*post_process_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight) = identity_edge_weight,
					  const bool verbose = false)
{
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

	bool is_multi_graph, contains_edge_weights, contains_vertex_weights;

	contains_vertex_weights = fmt & 1;
	contains_edge_weights   = fmt & 2;
	is_multi_graph    = fmt & 4;

	if(verbose)
	{
		std::cout << num_vertices << ", " << num_edges << ", " << (int)fmt << ", " << ncon << std::endl;
		std::cout << "contains_vertex_weights: " << contains_vertex_weights << ",  contains_edge_weights: " << contains_edge_weights << ",  is_multi_graph: " << is_multi_graph << std::endl;
	}

	if(!is_multi_graph)
	{ // num_edges already contains the correct value, according to the standard
		edges.reserve(num_edges);
	}

	std::vector<WeightType> node_weights(num_vertices, 1);

	NodeIDType src_vertex = 0; // I numerate vertices starting with 0!!

	std::vector<Edge> read_edges;
	// read edges
	while(infile.good())
	{
		infile.getline(line);

		if(line[0] != '%')
		{ // not a comment-line
			read_edges.clear();
			extract_edges_of_vertex<WeightType, NodeIDType, Edge>(line, src_vertex, ncon, contains_edge_weights, weight_func, read_edges, node_weights);

			// insert edges to graph
			for(typename std::vector<Edge>::iterator e_it=read_edges.begin(); e_it!=read_edges.end(); e_it++)
			{
				edges.push_back(*e_it);
			}

			src_vertex++;
		}
	}

	num_edges = edges.size();

	post_process_edge_weights<WeightType, EdgeIDType, Edge>(edges, node_weights, post_process_function);
}


} // end of namespace dipl


#endif
