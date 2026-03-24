#ifndef GRID_GRAPH_2D_H_
#define GRID_GRAPH_2D_H_

#include <iostream>

#ifdef _MPI
#include <mpi.h>
#endif

namespace dipl
{


// // only an estimation, return value in Byte
//template <class Graph>
//void mem_usg_for_hypercube_matching(const int procs, const int dim_length)
//{
//	typedef typename Graph::Edge Edge;
//	typedef typename Graph::NodeIDType NodeIDType;
//	typedef typename Graph::EdgeIDType EdgeIDType;
//	typedef typename Graph::WeightType WeightType;
//
//	int xprocs, yprocs;
//
//	get_2d_process_distribution(procs, xprocs, yprocs);
//
//	unsigned long max_x_length = dim_length/xprocs + 1;
//	unsigned long max_y_length = dim_length/yprocs + 1;
//
//
//	unsigned long size_of_creation_edge_vector = max_x_length*max_y_length*2*sizeof(Edge) + sizeof(Edge)*(max_x_length+max_y_length);
//
//
//	unsigned long size_of_graph_local_edge_vector = max_x_length*max_y_length*2*sizeof(Edge) - sizeof(Edge)*(max_x_length+max_y_length);
//	unsigned long size_of_graph_cross_edge_vector = 2 * sizeof(Edge)*(max_x_length+max_y_length);
//	unsigned long size_of_graph_proc_of_ghost_vertex = sizeof(int) * 2 * (max_x_length+max_y_length);
//	unsigned long size_of_graph_active_cross_edges_of_proc = sizeof(EdgeIDType) * procs;
//	unsigned long size_of_graph_ghost_global_to_local = sizeof(std::pair<NodeIDType, NodeIDType>) * 2 * (max_x_length+max_y_length);
//	unsigned long size_of_local_to_global = sizeof(NodeIDType) * 2 * (max_x_length+max_y_length);
//
//	unsigned long size_of_graph_tmp_ghost_vertices = sizeof(NodeIDType) * 2 * (max_x_length+max_y_length);
//
//	unsigned long graph_size =  size_of_graph_local_edge_vector + size_of_graph_cross_edge_vector + size_of_graph_proc_of_ghost_vertex +
//			size_of_graph_active_cross_edges_of_proc + size_of_graph_ghost_global_to_local +  size_of_local_to_global + size_of_graph_tmp_ghost_vertices;
//
//
//	// sizeof(Message) = sizeof(NodeIDType)*2
//	unsigned long size_of_matching_messages_of_proc = 2 * (max_x_length+max_y_length) * 2 * sizeof(NodeIDType) + procs*sizeof(MPI_Request);
//	unsigned long size_of_matching_matched = (max_x_length*max_y_length+2*(max_x_length+max_y_length)) / 8;
//	// sizeof(Candidate) = "sizeof(WeightType)+sizeof(NodeIDType)" = sizeof(std::pair<WeighType, NodeIDType>), because of allignment
//	unsigned long size_of_matching_candidate_of_node = sizeof(std::pair<WeightType, NodeIDType>)*(max_x_length*max_y_length+2*(max_x_length+max_y_length));
//	unsigned long size_of_matching_recv_requests = procs*sizeof(MPI_Request);
//	unsigned long size_of_matching_recv_msgs = 2*(max_x_length+max_y_length) * sizeof(NodeIDType)*2;
//
//	unsigned long matching_size = size_of_matching_messages_of_proc + size_of_matching_matched + size_of_matching_candidate_of_node + size_of_matching_candidate_of_node +
//			size_of_matching_recv_requests + size_of_matching_recv_msgs;
//
//	// approx, assuming matching has size one quarter of the total amount of edges
//	unsigned long size_of_matching = size_of_creation_edge_vector/4;
//
//
//
//	unsigned long init_plus_graph = size_of_creation_edge_vector + graph_size;
//	unsigned long graph_plus_matching = graph_size + matching_size + size_of_matching;
//	unsigned long total_size = size_of_creation_edge_vector + graph_size + matching_size + size_of_matching;
//
//
//	std::cout << "init + graph: " << init_plus_graph << ",  graph + matching: " << graph_plus_matching << ",  total: " << total_size << std::endl;
//}


void get_2d_process_distribution(const int procs, int &xprocs, int &yprocs)
{
	xprocs = sqrt(procs);

	while(procs%xprocs!=0)
	{
		xprocs--;
	}

	yprocs = procs/xprocs;
}


template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
class Grid2D
{
public:

	static const unsigned int dim = 2;

	static unsigned int dim_length;

	static int procs, xprocs, yprocs, proc_id;

	static unsigned int base_length_x, base_length_y, base_length_xp1, base_length_yp1;
	static unsigned int x_remainder, y_remainder;

	static NodeIDType dim_length_X_base_length_yp1, dim_length_X_base_length_y; // dim_length*(base_length_y+1) and dim_length*base_length_y, respectively
	static NodeIDType base_length_yp1_X_y_remainder, dim_length_X_base_length_yp1_X_y_remainder;
	static NodeIDType base_length_xp1_X_x_remainder;

public:

	/*!
	 *
	 * @param dl	dim length
	 * @param edges
	 * @param num_vertices
	 * @param first_global_vertex
	 * @param end_global_vertex
	 * @param proc_ranges
	 * @param weight_func
	 */
	static void get_grid(const unsigned int dl,
						 ContainerType &edge_container, NodeIDType &num_vertices,
			  	  	  	 NodeIDType &first_global_vertex, NodeIDType &end_global_vertex,
			  	  	  	 std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
			  	  	  	 WeightType (*weight_func)())
	{
		dim_length = dl;

		procs = 1;
		proc_id = 0;

#ifdef _MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Datatype MPI_EDGE_TYPE;
		MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &MPI_EDGE_TYPE);
		MPI_Type_commit(&MPI_EDGE_TYPE);

		std::vector<Edge> messages_to_right_process, messages_to_top_process;
#endif

		// set the process distribution
		get_2d_process_distribution(procs, xprocs, yprocs);

		int xp = proc_id % xprocs;
		int yp = proc_id / xprocs;

#ifdef _MPI
		int right_process  = yp*xprocs+xp+1;
		int top_process    = (yp+1)*xprocs+xp;

		int left_process   = yp*xprocs+xp-1;
		int bottom_process = (yp-1)*xprocs+xp;
#endif

		/// \todo use init-function to set those values, and other static members
		/*unsigned int*/ base_length_x = dim_length / (unsigned int)xprocs;
		/*unsigned int*/ base_length_y = dim_length / (unsigned int)yprocs;
		base_length_xp1 = base_length_x+1;
		base_length_yp1 = base_length_y+1;

		dim_length_X_base_length_y   = dim_length * base_length_y;
		dim_length_X_base_length_yp1 = dim_length * base_length_yp1;

		base_length_yp1_X_y_remainder = base_length_yp1 * y_remainder;
		dim_length_X_base_length_yp1_X_y_remainder = dim_length_X_base_length_yp1 * y_remainder;

		x_remainder = dim_length%xprocs;
		y_remainder = dim_length%yprocs;

		const NodeIDType proc_x_start = xp*base_length_x + ( xp<(int)(dim_length%xprocs) ? xp : (dim_length%xprocs) );
		const NodeIDType proc_y_start = yp*base_length_y + ( yp<(int)(dim_length%yprocs) ? yp : (dim_length%yprocs) );

		const NodeIDType proc_x_dim_length = base_length_x + (xp<(int)(dim_length%xprocs));
		const NodeIDType proc_y_dim_length = base_length_y + (yp<(int)(dim_length%yprocs));

		first_global_vertex = std::numeric_limits<NodeIDType>::max();
		end_global_vertex   = std::numeric_limits<NodeIDType>::min();

		unsigned long number_of_edges = proc_x_dim_length*proc_y_dim_length*2+proc_x_dim_length+proc_y_dim_length;

		// edges won't be larger than the following size
		edge_container.reserve(number_of_edges);

		for(NodeIDType y=proc_y_start; y<proc_y_start+proc_y_dim_length; y++)
		{
			for(NodeIDType x=proc_x_start; x<proc_x_start+proc_x_dim_length; x++)
			{
				NodeIDType nvid = y*dim_length + x;
				NodeIDType gvid = nvid_to_gvid(nvid);

				if(x<dim_length-1)
				{// right neighbor exists
					NodeIDType rp = nvid+1; // right neighbor
					NodeIDType neighbor_gvid = nvid_to_gvid(rp);

					Edge new_edge = Edge(neighbor_gvid,
										 gvid,
										 weight_func());

					edge_container.push_back(new_edge);

#ifdef _MPI
					// we only add edges to neighbors with a smaller ID or if the
					// neighbor is located on a different process. Local edes
					// should only be added once
					if(x==proc_x_start+proc_x_dim_length-1 && xp<xprocs-1)
					{
						messages_to_right_process.push_back(new_edge);
					}
#endif
				}


				if(y<dim_length-1)
				{// top neighbor exists
					NodeIDType tp = nvid+dim_length; // top neighbor
					NodeIDType neighbor_gvid = nvid_to_gvid(tp);

					Edge new_edge = Edge(neighbor_gvid,
										 gvid,
										 weight_func());

					edge_container.push_back(new_edge);

#ifdef _MPI
					// we only add edges to neighbors with a smaller ID or if the
					// neighbor is located on a different process. Local edes
					// should only be added once
					if(y==proc_y_start+proc_y_dim_length-1 && yp<yprocs-1)
					{
						messages_to_top_process.push_back(new_edge);
					}
#endif
				}


				// adjust first_global_vertex and end_global_vertex if necessary
				if(first_global_vertex>gvid) first_global_vertex = gvid;
				if(end_global_vertex  <gvid) end_global_vertex   = gvid;

			}
		}

		end_global_vertex++; // we have to increase this value by one

#ifdef _MPI
		MPI_Request req;

		// send cross edge_container to partners
		if(!messages_to_right_process.empty())
		{
			MPI_Isend(&messages_to_right_process[0], messages_to_right_process.size(),
					  MPI_EDGE_TYPE, right_process, 0, MPI_COMM_WORLD, &req);
		}

		if(!messages_to_top_process.empty())
		{
			MPI_Isend(&messages_to_top_process[0], messages_to_top_process.size(),
					  MPI_EDGE_TYPE, top_process, 0, MPI_COMM_WORLD, &req);
		}

		//receive cross edge_container
		if(xp>0)
		{// receive from left
			std::vector<Edge> messages(proc_y_dim_length);
			MPI_Recv(&messages[0], proc_y_dim_length, MPI_EDGE_TYPE, left_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(unsigned int i=0; i<messages.size(); i++)
			{
				edge_container.push_back(messages[i]);
			}
		}

		if(yp>0)
		{// receive from bottom
			std::vector<Edge> messages(proc_x_dim_length);
			MPI_Recv(&messages[0], proc_x_dim_length*sizeof(Edge), MPI_EDGE_TYPE, bottom_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(unsigned int i=0; i<messages.size(); i++)
			{
				edge_container.push_back(messages[i]);
			}
		}


		// exchange process ranges
		proc_ranges.resize(procs);

		proc_ranges[proc_id] = std::make_pair(first_global_vertex, end_global_vertex);

		MPI_Allgather(&proc_ranges[proc_id], sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
					  &proc_ranges[0],       sizeof(std::pair<NodeIDType, NodeIDType>), MPI_BYTE,
					  MPI_COMM_WORLD);
#endif

		num_vertices = dim_length*dim_length; // that's the global number of vertices!!!
	}


	// satic dim_length, xprocs and yprocs must be set before calling this function, e.g. call get_cube at first
	static inline NodeIDType nvid_to_gvid(const NodeIDType nvid)
	{
		int p = proc_of_nvid(nvid);

		int xp = p % xprocs;
		int yp = p / xprocs;

		unsigned int base_length_x = dim_length / xprocs;
		unsigned int base_length_y = dim_length / yprocs;

		NodeIDType proc_x_dim_length = base_length_x + (xp < (int)(dim_length%xprocs));
		NodeIDType proc_x_start      = (xp * base_length_x) + ((xp < (int)(dim_length%xprocs)) ? xp : dim_length%xprocs);

		NodeIDType proc_y_dim_length = base_length_y + (yp < (int)(dim_length%yprocs));
		NodeIDType proc_y_start      = (yp * base_length_y) + ((yp < (int)(dim_length%yprocs)) ? yp : dim_length%yprocs);

		NodeIDType proc_fst_vertex = dim_length*proc_y_start + proc_x_start*proc_y_dim_length;

		NodeIDType gvid = proc_fst_vertex +
						  proc_x_dim_length * (nvid/dim_length - proc_y_start) +
						  nvid%dim_length - proc_x_start;

		return gvid;
	}


	// satic dim_length, xprocs and yprocs must be set before calling this function, e.g. call get_cube at first
	static inline int proc_of_nvid(const NodeIDType nvid)
	{
		int p;

		unsigned int base_length_x = dim_length / (unsigned int)xprocs;
		unsigned int base_length_y = dim_length / (unsigned int)yprocs;

		// add full rows to p
		if(nvid > (dim_length*(base_length_y+1)*(dim_length%yprocs)))
		{
			p = xprocs*(dim_length%yprocs); // base_columns with length: (base_length+1)

			// base columns with length: base_length
			p += xprocs * ( (nvid-(dim_length*(base_length_y+1)*(dim_length%yprocs))) / (dim_length*base_length_y) );
		}
		else
		{
			p = xprocs * (nvid / (dim_length*(base_length_y+1)));
		}

		// add missing columns to p

		// base rows with length: base_length+1
		if((nvid%dim_length) > ((base_length_x+1)*(dim_length%xprocs)) )
		{
			p += dim_length%xprocs;
			p += ((nvid%dim_length) - ((base_length_x+1)*(dim_length%xprocs))) / base_length_x;
		}
		else
		{
			p += (nvid%dim_length) / (base_length_x+1);
		}

		return p;
	}


	// satic dim_length, xprocs and yprocs must be set before calling this function, e.g. call get_cube at first
	static inline int proc_of_gvid(NodeIDType gvid)
	{
		int p;

		unsigned int base_length_x = dim_length / (unsigned int)xprocs;
		unsigned int base_length_y = dim_length / (unsigned int)yprocs;
		unsigned int y_length;

		// add full rows to p
		if(gvid > (dim_length*(base_length_y+1)*(dim_length%yprocs)))
		{
			p = xprocs*(dim_length%yprocs); // base_columns with length: (base_length+1)

			gvid -= (dim_length%yprocs) * (base_length_y+1) * dim_length; // adjust gvid (subtract the rows of processes that we already added to p)

			// base columns with length: base_length
			p += xprocs * ( gvid / (dim_length*base_length_y) );

			gvid -= dim_length * base_length_y * (gvid/(dim_length*base_length_y)); // adjust gvid (subtract the rows of processes that we already added to p)
			y_length = base_length_y;
		}
		else
		{
			p = xprocs * (gvid / (dim_length*(base_length_y+1)));
			gvid -= dim_length * (base_length_y+1) * (gvid / (dim_length*(base_length_y+1)));
			y_length = base_length_y+1;
		}


		// add missing columns to p

		// base rows with length: base_length+1
		if(gvid > ((base_length_x+1)*y_length*(dim_length%xprocs)) )
		{
			p += dim_length%xprocs;

			gvid -= ((base_length_x+1)*y_length*(dim_length%xprocs));

			p += gvid / (y_length*base_length_x);
		}
		else
		{
			p += (gvid) / (y_length*(base_length_x+1));
		}

		return p;
	}


	static inline int proc_of_neighbor_gvid(NodeIDType gvid)
	{
		int p = proc_id;

		const int xp = proc_id % xprocs;
		const int yp = proc_id / xprocs;

		const bool x_extra = xp<x_remainder;
		const bool y_extra = yp<y_remainder;

		const NodeIDType x_length = base_length_x + x_extra;
		const NodeIDType x_start  = (xp * base_length_x) + (x_extra ? xp : x_remainder);

		const NodeIDType y_length = base_length_y + y_extra;
		const NodeIDType y_start  = (yp * base_length_y) + (y_extra ? yp : y_remainder);


		const NodeIDType nodes_below = y_start*dim_length;
		const NodeIDType nodes_left  = nodes_below + x_start*y_length;
		const NodeIDType nodes_right = nodes_left + x_length*y_length;
		const NodeIDType nodes_above = (y_start+y_length)*dim_length;

		p -= (gvid<nodes_below)*xprocs; // if bottom neighbor
		p -= (gvid>=nodes_below && gvid<nodes_left); // if left neighbor
		p += (gvid>=nodes_right && gvid<nodes_above); // if right neighbor
		p += (gvid>=nodes_above)*xprocs; // if top neighbor

		return p;
	}




	// satic dim_length, xprocs and yprocs must be set before calling this function, e.g. call get_cube at first
	static inline NodeIDType gvid_to_nvid(NodeIDType gvid)
	{
		unsigned int y_length, y_start, x_start;

		// add full rows to p
		if(gvid > (dim_length_X_base_length_yp1_X_y_remainder))
		{
			y_start = base_length_yp1_X_y_remainder;

			gvid -= y_start * dim_length; // adjust gvid (subtract the rows of processes that we already added to p)

			y_start += ( gvid / (dim_length_X_base_length_y) ) * base_length_y;

			gvid -= dim_length_X_base_length_y * (gvid/(dim_length_X_base_length_y)); // adjust gvid (subtract the rows of processes that we already added to p)

			y_length = base_length_y;
		}
		else
		{
			y_start = (gvid / (dim_length_X_base_length_yp1)) * (base_length_yp1);
			gvid -= dim_length * y_start;
			y_length = base_length_yp1;
		}

		// base rows with length: base_length+1
		if(gvid > (base_length_xp1_X_x_remainder*y_length) )
		{
			x_start = base_length_xp1_X_x_remainder;

			gvid -= x_start*y_length;

			x_start += (gvid/(base_length_x*y_length))*base_length_x;
			gvid -= (gvid/(base_length_x*y_length))*base_length_x*y_length;

			x_start += gvid%base_length_x;
			y_start += gvid/base_length_x;
		}
		else
		{
			x_start = (gvid / ((base_length_xp1)*y_length))*(base_length_xp1);
			gvid -= x_start*y_length;

			x_start += gvid%(base_length_xp1);
			y_start += gvid/(base_length_xp1);
		}

		return y_start*dim_length + x_start;
	}


};




template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::xprocs;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::yprocs;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::procs;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::dim_length;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::proc_id;


template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_x;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_y;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_xp1;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_yp1;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::x_remainder;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
unsigned int Grid2D<WeightType, NodeIDType, Edge, ContainerType>::y_remainder;


template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
NodeIDType Grid2D<WeightType, NodeIDType, Edge, ContainerType>::dim_length_X_base_length_y;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
NodeIDType Grid2D<WeightType, NodeIDType, Edge, ContainerType>::dim_length_X_base_length_yp1;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
NodeIDType Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_yp1_X_y_remainder;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
NodeIDType Grid2D<WeightType, NodeIDType, Edge, ContainerType>::dim_length_X_base_length_yp1_X_y_remainder;

template <typename WeightType, typename NodeIDType, class Edge, class ContainerType>
NodeIDType Grid2D<WeightType, NodeIDType, Edge, ContainerType>::base_length_xp1_X_x_remainder;

} // end ns dipl

#endif
