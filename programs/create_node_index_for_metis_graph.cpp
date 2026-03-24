
#include <vector>
#include <utility>

#include <common/linear_buffered_reader.h>
#include <common/math_funcs.h>



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



double get_node_weight(const std::string &line,
					   const int ncon)
{
	std::istringstream line_stream(line);

	// advance to start of edges and store the weight if we only have a single node weight

	if(ncon==1)
	{ // read the single node weight
		double w;
		line_stream >> w;

		return w;
	}

	return 1.;
}



void create_node_index_file(const std::string &input_filename,
					  	  	const std::string &output_filename,
					  	  	const bool verbose = false)
{
	std::string line;

	dipl::LinearBufferedReader<1000000> infile(input_filename);

	std::stringstream output_buf;

	std::ofstream outfile;
	outfile.open(output_filename.c_str(), std::ios_base::out | std::ios_base::binary);

	unsigned int ncon;
	dipl::byte fmt;

	// find/read the header line
	infile.getline(line);
	while(line[0]=='%')
	{
		infile.getline(line);
	}

	unsigned long num_vertices, num_edges;

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


	// read edges
	while(infile.good())
	{
		// store for each node line the start position
		if(line[0] != '%')
		{ // not a comment-line
			unsigned long current_byte_pos = infile.tellg();
			infile.getline(line);
			double weight = get_node_weight(line, ncon);

			std::pair<unsigned long, double> node_info = std::make_pair(current_byte_pos, weight);

			outfile.write((char*)&node_info, sizeof(std::pair<unsigned long, double>));
		}
	}

	outfile.close();
}



int main(int argc, char **argv)
{
	std::string input_filename, output_filename;

	if(argc==(1+1))
	{
		input_filename = argv[1];
	}
	else
	{
		std::cout << "error: Wrong number of arguments. First argument must be input filename and second argument output filename!" << std::endl;
		exit(1);
	}

	
	std::stringstream filename;
	filename << input_filename << ".index";
	output_filename = filename.str();

	create_node_index_file(input_filename, output_filename, true);

	return 0;
}




