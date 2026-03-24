
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

// the src_vertex is always the first end point of an edge
// this function adds all edges, i.e. also edges to nodes with a smaller ID
// than src_vertex. That's a different behavior compared to the one of
// in metis_dimas_reader.h
void store_rand_edge_weights_of_vertex(const std::string &line,
									   const unsigned int ncon,
									   const bool edge_weights,
									   std::ostream &output)
{
	std::istringstream line_stream(line);

	// store node weights, copying them to outfile
	for(unsigned int i=0; i<ncon; i++)
	{
		double dummy;
		line_stream >> dummy;
		output << dummy << " ";
	}


	while(line_stream.good())
	{
		long tgt_vertex;
		double weight = 0.;

		line_stream >> tgt_vertex;

		if(line_stream.fail())
		{ // couldn't read the last number, e.g. there wasn't one because the line
		  // ended with a whitespace, thus line_stream.good() didn't catch it
			break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
		}

		if(edge_weights)
		{
			line_stream >> weight;

			if(line_stream.fail())
			{ // couldn't read the last number, e.g. there wasn't one because the line
			  // ended with a whitespace, thus line_stream.good() didn't catch it
				break; // use break instead of continue, either it's the end of a line or the graph-file is corrupted
			}
		}

		weight = dipl::random_weight();

		output << tgt_vertex << " " << weight << " ";

	}

	output << std::endl;
}



void store_rand_weight_graph(const std::string &input_filename,
					  	  	 const std::string &output_filename,
					  	  	 const bool verbose = false)
{
	std::string line;

	dipl::LinearBufferedReader<1000000> infile(input_filename);

	unsigned long output_threshold = 100000;
	unsigned long signs_read = 0;

	std::stringstream output_buf;

	std::ofstream outfile;
	outfile.open(output_filename.c_str());

	unsigned int ncon;
	dipl::byte fmt;

	// find the header line
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


	unsigned int out_fmt;

	if(is_multi_graph && contains_vertex_weights)
	{
		out_fmt = 111;
	}
	else if(is_multi_graph)
	{
		out_fmt = 110;
	}
	if(contains_vertex_weights)
	{
		out_fmt = 11;
	}
	else
	{
		out_fmt = 10;
	}


	outfile << num_vertices << " " << num_edges
			<< " " << out_fmt << " " << ncon << std::endl;


	// read edges
	while(infile.good())
	{
		infile.getline(line);

		if(line[0] != '%')
		{ // not a comment-line
			store_rand_edge_weights_of_vertex(line, ncon, contains_edge_weights, output_buf);

			signs_read += line.length();

			if(signs_read>output_threshold)
			{
				outfile << output_buf.str();
				output_buf.str("");
				signs_read = 0;
			}
		}
	}

	outfile << output_buf.str();

}



int main(int argc, char **argv)
{
	std::string input_filename, output_filename;

	if(argc==(2+1))
	{
		input_filename = argv[1];
		output_filename = argv[2];
	}
	else
	{
		std::cout << "error: Wrong number of arguments. First argument must be input filename and second argument output filename!" << std::endl;
		exit(1);
	}

	store_rand_weight_graph(input_filename, output_filename, true);

	return 0;
}




