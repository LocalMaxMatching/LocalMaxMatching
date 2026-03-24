#ifndef GRID_GRAPH_H_
#define GRID_GRAPH_H_

#include "grid_graph_general.h"
#include "grid_graph_2d.h"

namespace dipl
{



template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
void get_grid(const unsigned int dim, const unsigned int dim_length,
			  ContainerType &edge_container, NodeIDType &num_vertices,
			  NodeIDType &first_global_vertex, NodeIDType &end_global_vertex,
			  std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
			  WeightType (*weight_func)())
{
	if(dim==2)
	{
		Grid2D<WeightType, NodeIDType, Edge, ContainerType>::get_grid(dim_length,
													   edge_container, num_vertices,
													   first_global_vertex, end_global_vertex,
													   proc_ranges, weight_func);
	}
	else
	{
		get_general_grid<WeightType, NodeIDType, Edge, ContainerType>(dim, dim_length,
						 edge_container, num_vertices,
						 first_global_vertex, end_global_vertex,
						 proc_ranges,
						 weight_func);
	}
}



}// end ns dipl

#endif

