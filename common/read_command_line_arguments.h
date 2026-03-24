/*
 * You have to link against libboost_program_options_
 * 	-lboost_program_options
 */

#ifndef READ_COMMAND_LINE_ARGS_H_
#define READ_COMMAND_LINE_ARGS_H_

#include <iostream>

#include <boost/program_options.hpp>


namespace dipl
{

void read_command_line_args(int argc, char **argv,
							bool &read_graph, std::string &path,
							bool &create_kn, unsigned int &num_vertices,
							bool &create_grid, unsigned int &dim, unsigned int &dim_length,
							std::string &rating, int &seed, unsigned int &repetitions,
							bool is_root=true,
							std::string *matching_output_file=NULL)
{
	// set default values
	read_graph  = false;
	create_kn   = false;
	create_grid = false;

	// read command line args
	boost::program_options::options_description desc("Allowed Options");

	desc.add_options()
			("help,h", "Prints this help-message.")
			("grid", "Generates a grid to be used as the input graph. You also have to specify --dim and --dim_length .")
			("dim", boost::program_options::value<unsigned int>(&dim), "Number of the dimensions of the grid.")
			("dim_length", boost::program_options::value<unsigned int>(&dim_length), "The length of each direction of the grid.")
			("kn", "Generates a complete graph (K_n-Graph) to be used as the input graph. You also have to specify the number of vertices of this graph (--vertices).")
			("vertices", boost::program_options::value<unsigned int>(&num_vertices), "Number of vertices of the graph.")
			("read_graph", boost::program_options::value<std::string>(&path), "Reads a graph from the provided file.")
			("edge_rating", boost::program_options::value<std::string>(&rating)->default_value("const"), "Specifies the edge rating. Possible arguments: const, weight, expansionstar2, rand .")
			("seed", boost::program_options::value<int>(&seed), "Provides the seed used for computing random values.")
			("repetitions", boost::program_options::value<unsigned int>(&repetitions)->default_value(1), "How often to repeat the computation of the matching.")
			("matching_output_file", boost::program_options::value<std::string>(), "Write matched edge pairs to the specified file.")
			;

	boost::program_options::variables_map vm;

	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc,argv, desc), vm);
		boost::program_options::notify(vm);
	}
	catch (std::exception &e)
	{
		if(is_root)
		{
			std::cout << e.what() << std::endl << std::endl;
			std::cout << desc;
		}
		exit(1);
	}

	if(!vm.count("seed") && (rating.compare("rand")==0))
	{// no seed provided but using random weights
		seed = time(NULL);
	}

	if(matching_output_file != NULL && vm.count("matching_output_file"))
	{
		*matching_output_file = vm["matching_output_file"].as<std::string>();
	}

	if(vm.count("help"))
	{
		if(is_root)
		{
			std::cout << desc;
		}
		exit(1);
	}
	else if(vm.count("grid"))
	{
		if(vm.count("dim") && vm.count("dim_length"))
		{
			create_grid = true;
			if(is_root)
			{
				std::cout << "Creating grid with " << dim << " dimensions, and a dimension-length of " << dim_length << std::endl;
			}
			return ;
		}
		else
		{
			if(is_root)
			{
				std::cout << desc;
			}
			exit(1);
		}
	}
	else if(vm.count("kn"))
	{
		if(vm.count("vertices"))
		{
			create_kn = true;
			if(is_root)
			{
				std::cout << "Creating complete graph with " << num_vertices << " vertices." << std::endl;
			}
			return ;
		}
		else
		{
			if(is_root)
			{
				std::cout << desc;
			}
			exit(1);
		}
	}
	else if(vm.count("read_graph"))
	{
		read_graph = true;
		if(is_root)
		{
			std::cout << "Reading graph from file " << path << std::endl;
		}
		return ;
	}
	else
	{
		if(is_root)
		{
			std::cout << desc;
		}
		exit(1);
	}
}


}

#endif


