/*
 * Read graphs in the Matrix-Market-Format. Assumes Coordinate Format
 * and treats the matrices as bipartite graphs.
 */

#ifndef DIPL_READ_MATRIX_MARKET_H_
#define DIPL_READ_MATRIX_MARKET_H_

#include <vector>
#include <utility>
#include <algorithm>

#include <graphs/adjacency_list.h>

#include "linear_buffered_reader.h"

namespace dipl
{



template <typename NodeIDType, typename EdgeIDType>
void extract_header(std::string &line,
					NodeIDType &rows,
					NodeIDType &columns,
					EdgeIDType &lines)
{
	std::istringstream line_stream(line);

	line_stream >> rows;
	line_stream >> columns;
	line_stream >> lines;
}

template <typename WeightType, typename NodeIDType>
void add_edge_to_graph(const std::string &line,
					   const NodeIDType num_rows, const NodeIDType num_cols,
					   dipl::AdjacencyList<NodeIDType, WeightType> &graph,
					   WeightType &minimal_weight)
{
	std::istringstream line_stream(line);

	NodeIDType n1, n2;
	WeightType weight;

	line_stream >> n1; n1--; // format starts counting at 1
	line_stream >> n2; n2--; // format starts counting at 1
	line_stream >> weight;

	if(weight<minimal_weight)
	{
		minimal_weight = weight;
	}

	// interleave the row and column nodes
	if(num_rows<=num_cols)
	{
		n1 *= 2;

		if(n2>=num_rows)
		{
			n2 = (2*num_rows)+(n2-num_rows);
		}
		else
		{
			n2 = (n2*2)+1;
		}

	}
	else
	{// num_rows > num_cols
		n2 = (n2*2)+1;

		if(n1>=num_cols)
		{
			n1 = (num_cols*2) + (n1-num_cols);
		}
		else
		{
			n1 *= 2;
		}
	}

	graph.add_incident_edge(n1, n2, weight);
	graph.add_incident_edge(n2, n1, weight);
}



template <typename WeightType, typename NodeIDType, typename EdgeIDType>
void get_matrix_market_graph(const std::string &filename,
					  dipl::AdjacencyList<NodeIDType, WeightType> &graph,
					  EdgeIDType &num_edges, WeightType &minimal_negative_weight)
{
	// we're only interested in the minimal negative weight
	minimal_negative_weight = 0;

	std::string line;
	LinearBufferedReader<1000000> infile(filename);

	// find the header line
	infile.getline(line);
	while(line[0]=='%')
	{
		infile.getline(line);
	}

	NodeIDType rows, cols;
	EdgeIDType num_lines;

	extract_header(line, rows, cols, num_lines);

	// creating a bipartite graph => rows+cols nodes
	graph.resize(rows+cols);

	num_edges = 0;

	// read edges
	while(infile.good())
	{
		infile.getline(line);
//		std::cout << line << std::endl;

		if(line[0] != '%')
		{ // not a comment-line
			num_edges++;
			add_edge_to_graph(line, rows, cols, graph, minimal_negative_weight);
		}
	}
}


} // end of namespace dipl


#endif
