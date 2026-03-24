
#include <vector>
#include <utility>
#include <algorithm>

#include <common/linear_buffered_reader.h>
#include <graphs/adjacency_list.h>


template <typename NodeIDType>
void extract_header(std::string &line,
					NodeIDType &n, NodeIDType &m,
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
template <typename WeightType, typename NodeIDType>
void extract_edges_of_vertex(const std::string &line, const NodeIDType src_vertex,
							 const NodeIDType ncon, const bool edge_weights,
							 dipl::AdjacencyList<NodeIDType, WeightType> &graph)
{
	std::istringstream line_stream(line);

	WeightType src_vertex_weight = 1.;

	if(ncon==1)
	{ // read the single node weight
		line_stream >> src_vertex_weight;
	}
	else
	{// no special treatment for multiple node weights
		for(unsigned int i=0; i<ncon; i++)
		{
			unsigned int dummy;
			line_stream >> dummy;
		}
	}

	graph.set_vertex_weight(src_vertex, src_vertex_weight);

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
			weight = 1.;
		}

		graph.add_incident_edge(src_vertex,
								tgt_vertex, weight);

	}
}


template <typename WeightType, typename NodeIDType>
void get_dimacs_graph(const std::string &filename,
					  NodeIDType &num_edges,
					  dipl::AdjacencyList<NodeIDType, WeightType> &graph,
					  const bool verbose = false)
{
	std::string line;

	dipl::LinearBufferedReader<1000000> infile(filename);

	unsigned int ncon;
	dipl::byte fmt;

	// find the header line
	infile.getline(line);
	while(line[0]=='%')
	{
		infile.getline(line);
	}

	NodeIDType num_vertices;

	extract_header(line, num_vertices, num_edges, fmt, ncon);


	graph.resize(num_vertices);

	bool is_multi_graph, contains_edge_weights, contains_vertex_weights;

	contains_vertex_weights = fmt & 1;
	contains_edge_weights   = fmt & 2;
	is_multi_graph    = fmt & 4;

	if(verbose)
	{
		std::cout << num_vertices << ", " << num_edges << ", " << (int)fmt << ", " << ncon << std::endl;
		std::cout << "contains_vertex_weights: " << contains_vertex_weights << ",  contains_edge_weights: " << contains_edge_weights << ",  is_multi_graph: " << is_multi_graph << std::endl;
	}


	NodeIDType src_vertex = 0; // I numerate vertices starting with 0!!

	// read edges
	while(infile.good())
	{
		infile.getline(line);

		if(line[0] != '%')
		{ // not a comment-line

			extract_edges_of_vertex<WeightType, NodeIDType>(line, src_vertex, ncon,
														    contains_edge_weights,
															graph);
			src_vertex++;
		}
	}

}



unsigned long bit_interleave(const unsigned int n1, const unsigned int n2)
{
	unsigned long result = 0;

	for(unsigned int i=0; i<(sizeof(unsigned int)*8); i++)
	{
		unsigned long tmp1, tmp2;
		tmp1 = n1 & (1 << i);
		tmp2 = n2 & (1 << i);
		result |= (tmp1 << i) | (tmp2 << (i+1));
	}

	return result;
}


void set_curve_values(const std::string &filename,
					  std::vector<std::pair<unsigned int, unsigned long> > &curve_values)
{
	std::stringstream extractor;
	std::string line;
	dipl::LinearBufferedReader<1000000> infile(filename);

	for(unsigned int i=0; i<curve_values.size(); i++)
	{
		infile.getline(line);
		extractor << line;

		unsigned int x,y;
		extractor >> x;
		extractor >> y;

		curve_values[i] = std::make_pair(i, bit_interleave(x, y));
	}
}


void set_new_node_numbering(const std::vector<std::pair<unsigned int, unsigned long> > &curve_values,
							std::vector<unsigned int> &new_node_number)
{
	for(unsigned int i=0; i<curve_values.size(); i++)
	{
		new_node_number[curve_values[i].first] = i;
	}
}


template <typename NodeIDType, typename WeightType>
void write_graph_using_curve_order(const std::string &filename, unsigned int num_edges,
								   const dipl::AdjacencyList<NodeIDType, WeightType> &graph,
								   const std::vector<std::pair<unsigned int, unsigned long> > &curve_values,
								   const std::vector<unsigned int> &new_node_number)
{
	unsigned int one_per_cent = curve_values.size()/100;
	unsigned int percent = 0;

	std::ofstream outfile(filename.c_str());
	if(outfile.good())
	{
		// write header
		outfile << graph.num_vertices() << " " << num_edges << " 11" << std::endl;

		// write nodes
		for(unsigned int i=0; i<curve_values.size(); i++)
		{
			if(i%one_per_cent == 0)
			{
				std::cout << percent++ << "%  ";
				std::cout.flush();
			}
			NodeIDType current_node = curve_values[i].first;

			// write node weight
			outfile << graph.get_node_weight(current_node);

			for(typename dipl::AdjacencyList<NodeIDType, WeightType>::incident_edge_iterator e_it=graph.begin_incident_edges(current_node);
																					e_it!=graph.end_incident_edges(current_node);
																					e_it++)
			{
				outfile << " " << (new_node_number[e_it->first]+1) << " " << e_it->second;
			}
			outfile << std::endl;
		}

		std::cout << std::endl;
	}
}


bool pair_compare(const std::pair<unsigned int, unsigned long> &l,
				  const std::pair<unsigned int, unsigned long> &r)
{
	return l.second <= r.second;
}


int main(int argc, char **argv)
{
	std::string infile_graph, infile_coords, outfile_graph;

	if(argc!=4) return 1;

	infile_graph = argv[1];
	infile_coords = argv[2];
	outfile_graph = argv[3];


	dipl::AdjacencyList<unsigned int, double> g;

	unsigned int num_edges;

	std::cout << "Reading graph " << infile_graph << std::endl;

	get_dimacs_graph(infile_graph, num_edges, g);

	std::cout << "Read graph " << infile_graph << std::endl;

	std::vector<std::pair<unsigned int, unsigned long> > curve_values(g.num_vertices());

	set_curve_values(infile_coords, curve_values);

	// get new order for nodes
	std::sort(curve_values.begin(), curve_values.end(), pair_compare);


	std::vector<unsigned int> new_node_number(g.num_vertices());

	set_new_node_numbering(curve_values, new_node_number);


	write_graph_using_curve_order(outfile_graph, num_edges,
								  g,
								  curve_values,
								  new_node_number);

	return 0;
}


