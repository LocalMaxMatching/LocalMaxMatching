#ifndef K_N_GRAPH_H_
#define K_N_GRAPH_H_

#include <vector>
#include <cmath>
#include <limits>
#include <utility>

#include <common/edge.h>
#include <common/math_funcs.h>


namespace dipl
{

template <typename NodeIDType, typename WeightType, class ContainerType>
class KNGraph
{
public:
	typedef dipl::Edge<WeightType, NodeIDType> Edge;

	static const int kn_graph_tag = 0;

public:


	static int proc_id_of_vertex(NodeIDType vid, NodeIDType num_vertices,
								 int procs)
	{
		int proc_id;

		NodeIDType base_number_of_proc_vertices = num_vertices / procs;
		NodeIDType remainder = num_vertices % procs;

		// Processes with proc-ranges smaller than threshold_id, have
		// local vertex count of base_number_of_proc_vertices+1
		NodeIDType threshold_id = remainder * (base_number_of_proc_vertices+1);

		if(vid<threshold_id)
		{
			proc_id = vid/(base_number_of_proc_vertices+1);
		}
		else
		{
			proc_id = remainder + ((vid-threshold_id)/base_number_of_proc_vertices);
		}

		return proc_id;
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
	static void get_graph(const NodeIDType num_vertices,
						 ContainerType &edge_container,
			  	  	  	 NodeIDType &first_global_vertex, NodeIDType &end_global_vertex,
			  	  	  	 std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
			  	  	  	 WeightType (*weight_func)())
	{
		int procs = 1;
		int proc_id = 0;

#ifdef _MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Datatype MPI_EDGE_TYPE;
		MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &MPI_EDGE_TYPE);
		MPI_Type_commit(&MPI_EDGE_TYPE);

		std::vector< std::vector<Edge> > msgs_for_proc(procs);
		std::vector<MPI_Request> reqs(procs, MPI_REQUEST_NULL);
#endif

		NodeIDType base_number_of_proc_vertices = num_vertices / procs;
		NodeIDType remainder = num_vertices % procs;


		first_global_vertex = base_number_of_proc_vertices*proc_id + (proc_id<(int)remainder ? proc_id : remainder);

		NodeIDType number_of_local_vertices = base_number_of_proc_vertices + (proc_id<(int)remainder);
		end_global_vertex = first_global_vertex + number_of_local_vertices;

		unsigned long number_of_edges = (number_of_local_vertices*(number_of_local_vertices-1)/2)
										+ number_of_local_vertices*(num_vertices-number_of_local_vertices);

		edge_container.reserve(number_of_edges);

		// add local edge_container
		for(NodeIDType nid=first_global_vertex; nid<end_global_vertex; nid++)
		{
			for(NodeIDType neighbor_id=nid+1; neighbor_id<end_global_vertex; neighbor_id++)
			{
				Edge new_edge = Edge(nid,
									 neighbor_id,
									 weight_func());

				edge_container.push_back(new_edge);
			}
		}

#ifdef _MPI
		std::vector<int> number_of_msgs(procs, 0);

		// add cross edges, only to neighbors with a larger ID - we'll receive the other cross edges in the next step
		// add local edge_container
		for(NodeIDType nid=first_global_vertex; nid<end_global_vertex; nid++)
		{
			for(NodeIDType neighbor_id=end_global_vertex; neighbor_id<num_vertices; neighbor_id++)
			{
				Edge new_edge = Edge(nid,
									 neighbor_id,
									 weight_func());

				edge_container.push_back(new_edge);

				int neighbor_proc_id = proc_id_of_vertex(neighbor_id, num_vertices, procs);
				number_of_msgs[neighbor_proc_id]++;
			}
		}


		for(int p=0; p<procs;p++)
		{
			msgs_for_proc[p].reserve(number_of_msgs[p]);
		}

		for(unsigned long i=0; i<edge_container.size(); i++)
		{
			const Edge &e = edge_container[i];

			if(e.n2>=end_global_vertex)
			{
				int neighbor_proc_id = proc_id_of_vertex(e.n2, num_vertices, procs);

				msgs_for_proc[neighbor_proc_id].push_back(e);
			}

		}


		// send cross edge_container to partners
		for(int p=0; p<procs; p++)
		{
			if(p!=proc_id)
			{
				if(!msgs_for_proc[p].empty())
				{// send cross edge_container to p
					MPI_Isend(&msgs_for_proc[p][0], msgs_for_proc[p].size(), MPI_EDGE_TYPE,
							  p, kn_graph_tag, MPI_COMM_WORLD, &reqs[p]);
				}
				else
				{// send empty message
					MPI_Isend(0, 0, MPI_EDGE_TYPE,
							  p, kn_graph_tag, MPI_COMM_WORLD, &reqs[p]);
				}
			}
		}


		//receive cross edge_container
		for(int p=0; p<procs; p++)
		{
			if(p!=proc_id)
			{
				MPI_Status status;
				MPI_Probe(p, kn_graph_tag, MPI_COMM_WORLD, &status);

				int msg_count;
				MPI_Get_count(&status, MPI_EDGE_TYPE, &msg_count);

				std::vector<Edge> incoming_edges(msg_count);

				MPI_Recv(&incoming_edges[0], msg_count, MPI_EDGE_TYPE, p, kn_graph_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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



};



}// end of namespace dipl


#endif

