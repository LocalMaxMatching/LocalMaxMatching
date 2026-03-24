#ifndef GET_GRAPH_H_
#define GET_GRAPH_H_

#include <common/parallel_metis_dimacs_reader.h>
//#include <common/parallel_metis_dimacs_reader_evenly_distributed.h>
#include <common/k_n_graph.h>
#include <common/grid_graph.h>
#include <common/rating_functions.h>
#include <common/math_funcs.h>

namespace dipl
{


// non parallel case
template <typename WeightType, typename NodeIDType, typename EdgeIDType, typename Edge, class ContainerType>
void get_edge_graph(const bool read_graph, const std::string &path,
					const bool create_kn,
					const bool create_grid, const unsigned int dim, const unsigned int dim_length,
					const std::string &rating,
					const int seed,
					const bool is_root,
					NodeIDType &num_vertices,
					EdgeIDType &num_edges,
					ContainerType &edge_container)
{
	// only dummy varibales
	NodeIDType first_global_vertex, end_global_vertex;
	std::vector<std::pair<NodeIDType, NodeIDType> > proc_ranges;

	if(read_graph)
	{ // read the graph from the provided input file
		WeightType (*rating_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight);

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Using edge-weight for edge-rating." << std::endl;
			rating_function = dipl::identity_edge_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Using expansionstar2 for edge-rating." << std::endl;
			rating_function = dipl::expansionstar2;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::rand_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_1_weight;
		}

		dipl::get_dimacs_graph<WeightType, NodeIDType, EdgeIDType, Edge>(path,
																		 num_vertices, num_edges,
																		 edge_container,
																		 first_global_vertex, end_global_vertex,
																		 proc_ranges,
																		 dipl::const_weight, rating_function);
	}
	else if(create_kn)
	{
		WeightType (*rating_function)();

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Complete graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Complete graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::random_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_weight;
		}

		dipl::KNGraph<NodeIDType, WeightType, ContainerType>::get_graph(num_vertices,
														 edge_container,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 first_global_vertex, end_global_vertex,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 proc_ranges,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 rating_function);
	}
	else if(create_grid)
	{
		WeightType (*rating_function)();

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Grid graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Grid graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::random_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_weight;
		}

		dipl::get_grid<WeightType, NodeIDType, Edge, ContainerType>(dim, dim_length, edge_container, num_vertices, first_global_vertex, end_global_vertex, proc_ranges, rating_function);
	}
	else
	{
		if(is_root) std::cout << "Wrong arguments specified? Use --help for more information.";
		exit(1);
	}
}






template <typename WeightType, typename NodeIDType, typename EdgeIDType, typename Edge, class ContainerType>
void get_edge_graph(const bool read_graph, const std::string &path,
					const bool create_kn,
					const bool create_grid, const unsigned int dim, const unsigned int dim_length,
					const std::string &rating,
					const int seed,
					const bool is_root,
					NodeIDType &num_vertices,
					EdgeIDType &num_edges,
					ContainerType &edge_container,
					NodeIDType &first_global_vertex, NodeIDType &end_global_vertex,
					std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges)
{
	if(read_graph)
	{ // read the graph from the provided input file
		WeightType (*rating_function)(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight);

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Using edge-weight for edge-rating." << std::endl;
			rating_function = dipl::identity_edge_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Using expansionstar2 for edge-rating." << std::endl;
			rating_function = dipl::expansionstar2;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::rand_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_1_weight;
		}

		dipl::get_dimacs_graph<WeightType, NodeIDType, EdgeIDType, Edge>(path,
																		 num_vertices, num_edges,
																		 edge_container,
																		 first_global_vertex, end_global_vertex,
																		 proc_ranges,
																		 dipl::const_weight, rating_function);
	}
	else if(create_kn)
	{
		WeightType (*rating_function)();

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Complete graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Complete graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::random_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_weight;
		}

		dipl::KNGraph<NodeIDType, WeightType, ContainerType>::get_graph(num_vertices,
														 edge_container,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 first_global_vertex, end_global_vertex,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 proc_ranges,
	 	  	  	 	 	 	 	 	 	 	 	 	 	 rating_function);
	}
	else if(create_grid)
	{
		WeightType (*rating_function)();

		if(rating.compare("weight")==0)
		{
			if(is_root) std::cout << "Grid graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("expansionstar2")==0)
		{
			if(is_root) std::cout << "Grid graphs only support constant weight and random weight at the moment. Using constant weight (1) instead." << std::endl;
			rating_function = dipl::const_weight;
		}
		else if(rating.compare("rand")==0)
		{
			srand(seed); // set seed for random number generation

			if(is_root) std::cout << "Using random weight for edge-rating with seed=" << seed << " ." << std::endl;
			rating_function = dipl::random_weight;
		}
		else
		{
			if(is_root) std::cout << "Using constant weight (1) for edge-rating - solving the cardinality problem." << std::endl;
			rating_function = dipl::const_weight;
		}

		dipl::get_grid<WeightType, NodeIDType, Edge, ContainerType>(dim, dim_length, edge_container, num_vertices, first_global_vertex, end_global_vertex, proc_ranges, rating_function);
	}
	else
	{
		if(is_root) std::cout << "Wrong arguments specified? Use --help for more information.";
		exit(1);
	}
}

} // end of ns dipl

#endif
