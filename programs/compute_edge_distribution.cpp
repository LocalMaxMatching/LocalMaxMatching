#include <vector>
#include <utility>
#include <algorithm>

#include <common/linear_buffered_reader.h>



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


/*!
 * \brief Returns the number of incident edges of the provided vertex that are
 * added to the edge array.
 *
 * All cross edges will be count, but only local edges to higher vertices
 * will be count, otherwise we would count the local edges twice
 *
 */
template <typename NodeIDType, typename EdgeIDType>
void number_of_contributing_edges(const std::string &line, const NodeIDType src_vertex,
										const NodeIDType ncon, const bool edge_weights,
										const NodeIDType first_global_vertex,
										const NodeIDType end_global_vertex,
										EdgeIDType &local_edges,
										EdgeIDType &cross_edges)
{
	cross_edges = 0;
	local_edges = 0;

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
			cross_edges++;
		}
		else if(tgt_vertex>=src_vertex)
		{
			local_edges++;
		}

	}
}





template <typename NodeIDType, typename EdgeIDType>
void get_edge_distribution(const std::string &filename,
					  const int procs,
					  std::vector<std::pair<EdgeIDType, EdgeIDType> > &edges_of_proc,
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
	EdgeIDType num_edges;

	extract_header(line, num_vertices, num_edges, fmt, ncon);

	bool is_multi_graph, contains_edge_weights, contains_vertex_weights;

	contains_vertex_weights = fmt & 1;
	contains_edge_weights   = fmt & 2;
	is_multi_graph    = fmt & 4;

	// there are some problems if NodeIDType is a floating point type
	NodeIDType base_size = num_vertices/procs;

	if(verbose)
	{
		std::cout << num_vertices << ", " << num_edges << ", " << (int)fmt << ", " << ncon << std::endl;
		std::cout << "contains_vertex_weights: " << contains_vertex_weights << ",  contains_edge_weights: " << contains_edge_weights << ",  is_multi_graph: " << is_multi_graph << std::endl;
	}

	NodeIDType src_vertex = 0;

	edges_of_proc.resize(procs, std::make_pair<EdgeIDType, EdgeIDType>(0,0));

	int current_proc = 0;
	NodeIDType first_global_vertex = 0;
	NodeIDType end_global_vertex = first_global_vertex + base_size + (current_proc < (int)(num_vertices % procs));

	while(infile.good() && src_vertex<num_vertices)
	{
		infile.getline(line);

		if(line[0] != '%')
		{ // not a comment-line
			EdgeIDType local_edges = 0;
			EdgeIDType cross_edges = 0;
			number_of_contributing_edges<NodeIDType, EdgeIDType>(line, src_vertex,
															ncon, contains_edge_weights,
															first_global_vertex,
															end_global_vertex,
															local_edges,
															cross_edges);

			edges_of_proc[current_proc].first  += local_edges;
			edges_of_proc[current_proc].second += cross_edges;

			src_vertex++;
		}

		if(src_vertex>=end_global_vertex)
		{
			current_proc++;
			first_global_vertex = end_global_vertex;
			end_global_vertex = first_global_vertex + base_size + (current_proc < (int)(num_vertices % procs));
		}
	}


}


int main(int argc, char **argv)
{
	typedef unsigned int EdgeIDType;
	typedef unsigned int NodeIDType;

	std::string input_filename;
	int procs;

	if(argc==(1+2))
	{
		input_filename = argv[1];
		std::stringstream in;
		in << argv[2];
		in >> procs;
	}
	else
	{
		std::cout << "error: Wrong number of arguments. First argument must be input filename and second argument number of processes!" << std::endl;
		exit(1);
	}


	std::vector<std::pair<EdgeIDType, EdgeIDType> > edges_of_proc;


	get_edge_distribution<NodeIDType, EdgeIDType>(input_filename,
												  procs,
												  edges_of_proc);

	for(int p=0; p<procs;p++)
	{
		std::cout << p << "\t" << edges_of_proc[p].first << "\t" << edges_of_proc[p].second << "\t" << (edges_of_proc[p].first+edges_of_proc[p].second) << std::endl;
	}

	return 0;
}



