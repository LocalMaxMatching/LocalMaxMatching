#ifndef GRID_GRAPH_GENERAL_H_
#define GRID_GRAPH_GENERAL_H_

#include "math_funcs.h"

#ifdef _MPI
#include <mpi.h>
#endif


namespace dipl
{

const int grid_tag = 0;




template <typename NodeIDType>
NodeIDType my_pow(unsigned int dim_length, unsigned int exp)
{
	NodeIDType result = 1;

	for(unsigned int i=0; i<exp; i++)
	{
		result *= dim_length;
	}

	return result;
}

/*!
 *
 * Doesn't return an optimal process distribution!!!
 *
 * @param dl	dim length
 * @param edge_container
 * @param num_vertices
 * @param first_global_vertex
 * @param end_global_vertex
 * @param proc_ranges
 * @param weight_func
 */
template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
void get_general_grid(const unsigned int dim, const unsigned int dl,
					  ContainerType &edge_container, NodeIDType &num_vertices,
					  NodeIDType &first_global_vertex, NodeIDType &end_global_vertex,
					  std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
					  WeightType (*weight_func)())
{
	unsigned int dim_length = dl;

	int procs = 1;
	int proc_id = 0;

#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	std::vector< std::vector<Edge> > msgs_for_proc(procs);
	std::vector<MPI_Request> reqs(procs, MPI_REQUEST_NULL);

	MPI_Datatype MPI_EDGE_TYPE;
	MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &MPI_EDGE_TYPE);
	MPI_Type_commit(&MPI_EDGE_TYPE);
#endif

	num_vertices = my_pow<NodeIDType>(dim_length, dim);

	NodeIDType base_number_of_proc_vertices = num_vertices / procs;
	NodeIDType remainder = num_vertices % procs;


	first_global_vertex = base_number_of_proc_vertices*proc_id + (proc_id<(int)remainder ? proc_id : remainder);

	NodeIDType number_of_local_vertices = base_number_of_proc_vertices + (proc_id<(int)remainder);
	end_global_vertex = first_global_vertex + number_of_local_vertices;




//		std::vector<Edge> messages_to_right_process, messages_to_top_process;


	// edges won't be larger than the following size
	edge_container.reserve(number_of_local_vertices*3);

	for(NodeIDType nid=first_global_vertex; nid<end_global_vertex; nid++)
	{
		NodeIDType dist_to_neighbor = 1;
		// iterate over all neighbors with a larger id
		for(unsigned int i=0; i<dim; i++)
		{
			NodeIDType neighbor = nid + dist_to_neighbor;
			dist_to_neighbor *= dim_length;

			// check for valid neighbor
			if((neighbor%dist_to_neighbor) > (nid%dist_to_neighbor))
			{
				Edge new_edge = Edge(nid,
									 neighbor,
									 weight_func());

				edge_container.push_back(new_edge);

#ifdef _MPI
				if(neighbor<first_global_vertex || neighbor>=end_global_vertex)
				{// neighbor is on different process - send information about this edge to neighbor process
					int neighbor_proc_id = proc_id_of_vertex(neighbor, num_vertices, procs);
					msgs_for_proc[neighbor_proc_id].push_back(new_edge);
				}
#endif
			}
		}

	}

#ifdef _MPI
	// send cross edge_container to partners
	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			if(!msgs_for_proc[p].empty())
			{// send cross edge_container to p
				MPI_Isend(&msgs_for_proc[p][0], msgs_for_proc[p].size(), MPI_EDGE_TYPE,
						  p, grid_tag, MPI_COMM_WORLD, &reqs[p]);
			}
			else
			{// send empty message
				MPI_Isend(0, 0, MPI_EDGE_TYPE,
						  p, grid_tag, MPI_COMM_WORLD, &reqs[p]);
			}
		}
	}


	//receive cross edge_container
	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			MPI_Status status;
			MPI_Probe(p, grid_tag, MPI_COMM_WORLD, &status);

			int msg_count;
			MPI_Get_count(&status, MPI_EDGE_TYPE, &msg_count);

			std::vector<Edge> incoming_edges(msg_count);

			MPI_Recv(&incoming_edges[0], msg_count, MPI_EDGE_TYPE, p, grid_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(int m=0; m<msg_count; m++)
			{
				edge_container.push_back(incoming_edges[m]);
			}
		}
	}

	MPI_Waitall(procs, &reqs[0], MPI_STATUS_IGNORE);

	// exchange process ranges
	proc_ranges.resize(procs);

	proc_ranges[proc_id] = std::make_pair(first_global_vertex, end_global_vertex);

	MPI_Datatype PairType;

	MPI_Type_contiguous(sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE, &PairType);


	MPI_Allgather(&proc_ranges[proc_id], sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
				  &proc_ranges[0],       sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
				  MPI_COMM_WORLD);
#endif
}

}// end of ns dipl

#endif
