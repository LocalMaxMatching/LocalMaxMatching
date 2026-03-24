#ifndef _PARALLEL_EDGE_GRAPH_
#define _PARALLEL_EDGE_GRAPH_


#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
#include <utility>

#include <boost/unordered_map.hpp>

#include <mpi.h>

#include <common/edge.h>
#include <common/union_find.h>
#include <common/math_funcs.h>


namespace dipl
{

template <typename T>
void swap(const unsigned long &pos1, const unsigned long &pos2,
		  std::vector<T> &vec)
{
	T tmp = vec[pos1];
	vec[pos1] = vec[pos2];
	vec[pos2] = tmp;
}

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1 nodes
/// _EdgeIDType should be at least unsigned long if there are more than 2^32 - 1 edges
template <typename _WeightType=double,
		  typename _NodeIDType=unsigned int,
		  typename _EdgeIDType=unsigned int>
class ParallelEdgeGraph
{

public:

	typedef _WeightType WeightType;
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

	typedef dipl::Edge<WeightType, NodeIDType> Edge;

	typedef typename std::vector<Edge>::iterator		edge_iterator;
	typedef typename std::vector<Edge>::const_iterator	const_edge_iterator;

//	class edge_iterator;

	typedef std::list<NodeIDType> 			ConnectedComponent;
	typedef std::list<ConnectedComponent> 	ConnectedComponents;


public:

	/*!
	 * \brief Constructor that initializes the graph from the given edge list.
	 *
	 * The constructor also ensures that the first node (n1) of each cross-edge
	 * is local and the second node (n2) is a ghost node!!
	 *
	 * Runtime: O(m+n)
	 *
	 * @param node_count	Number of vertices of the graph, defined by the edge list.
	 * @param edge_list		List of edges that define the graph. The vertices of the
	 * 						graph have to be numbered from 0 to node_count-1
	 */
	ParallelEdgeGraph(const NodeIDType node_count, const std::vector<Edge> &local_and_cross_edges,
					  const NodeIDType first_global_vertex, const NodeIDType end_global_vertex,
					  const std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
					  const bool delete_provided_edges=false);

	ParallelEdgeGraph()
	: initialized(false), dummy_vertex(std::numeric_limits<NodeIDType>::max())
	{
		dummy_edge.n1 = std::numeric_limits<NodeIDType>::max();
		dummy_edge.n2 = std::numeric_limits<NodeIDType>::max();
		dummy_edge.weight = min_val<WeightType>();
	}

	NodeIDType num_all_local_vertices() const
	{
		return _num_all_local_vertices;
	}


	NodeIDType num_local_vertices() const
	{
		return _num_local_vertices;
	}


	NodeIDType num_local_ghost_vertices() const
	{
		return _num_local_ghost_vertices;
	}

	EdgeIDType num_all_local_edges() const
	{
		return _num_all_local_edges;
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	const Edge& get_edge(const edge_iterator &edge_it) const
	{
		return *edge_it;
	}

	const Edge& get_edge(const_edge_iterator &edge_it) const
	{
		return *edge_it;
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	Edge& get_edge(const edge_iterator &edge_it)
	{
		return *edge_it;
	}



	/*!
	 * \brief Returns the edge corresponding to the given edge_id.
	 *
	 * Make sure that the edge_id is still valid. It's valid as long
	 * as you don't deactivate any edges. I.e. deactivating edges
	 * makes currently use edge_ids invalid (it least some of them).
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	const Edge& get_edge(const EdgeIDType &edge_id) const
	{
		return edges[edge_id];
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id.
	 *
	 * Make sure that the edge_id is still valid. It's valid as long
	 * as you don't deactivate any edges. I.e. deactivating edges
	 * makes currently use edge_ids invalid (it least some of them).
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	Edge& get_edge(const EdgeIDType &edge_id)
	{
		return edges[edge_id];
	}


	EdgeIDType get_edge_id(const edge_iterator &edge_it) const
	{
		return edge_it-edges.begin();
	}

	EdgeIDType get_edge_id(const const_edge_iterator &edge_it) const
	{
		return edge_it-edges.begin();
	}

	EdgeIDType get_edge_id(const edge_iterator &edge_it)
	{
		return edge_it-edges.begin();
	}


	/*!
	 * \brief Returns the weight of the edge referenced by the given iterator.
	 *
	 * Don't use edge references!!!
	 *
	 * @param edge_it
	 * @return			Returns the weight of the edge referenced by the given iterator.
	 */
	WeightType get_edge_weight(const edge_iterator &edge_it) const
	{
		return edge_it->weight;
	}

	WeightType get_edge_weight(const edge_iterator &edge_it)
	{
		return edge_it->weight;
	}


	edge_iterator begin_active_local_edges()
	{
		return edges.begin();//get_edge_iterator(0);//edge_iterator(0, edges);
	}

	const_edge_iterator begin_active_local_edges() const
	{
		return edges.begin();//get_edge_iterator(0);//edge_iterator(0, edges);
	}



	edge_iterator end_active_local_edges()
	{
		return get_edge_iterator(_end_active_local_edges);//edge_iterator(_end_active_local_edges, edges);
	}

	const_edge_iterator end_active_local_edges() const
	{
		return get_edge_iterator(_end_active_local_edges);//edge_iterator(_end_active_local_edges, edges);
	}


	edge_iterator begin_active_local_cross_edges()
	{
		return get_edge_iterator(_start_cross_edges);//edge_iterator(_start_cross_edges, edges);
	}

	const_edge_iterator begin_active_local_cross_edges() const
	{
		return get_edge_iterator(_start_cross_edges);//edge_iterator(_start_cross_edges, edges);
	}



	edge_iterator end_active_local_cross_edges()
	{
		return get_edge_iterator(_end_active_local_cross_edges);//edge_iterator(_end_active_local_cross_edges, edges);
	}

	const_edge_iterator end_active_local_cross_edges() const
	{
		return get_edge_iterator(_end_active_local_cross_edges);//edge_iterator(_end_active_local_cross_edges, edges);
	}



	edge_iterator begin_local_edges()
	{
		return edges.begin();//edge_iterator(0, edges);
	}

	const_edge_iterator begin_local_edges() const
	{
		return edges.begin();//edge_iterator(0, edges);
	}



	edge_iterator end_local_edges()
	{
		return get_edge_iterator(_num_local_edges);//edge_iterator(_num_local_edges, edges);
	}

	const_edge_iterator end_local_edges() const
	{
		return get_edge_iterator(_num_local_edges);//edge_iterator(_num_local_edges, edges);
	}



	edge_iterator begin_local_cross_edges()
	{
		return get_edge_iterator(_start_cross_edges);//edge_iterator(_start_cross_edges, edges);
	}

	const_edge_iterator begin_local_cross_edges() const
	{
		return get_edge_iterator(_start_cross_edges);//edge_iterator(_start_cross_edges, edges);
	}



	edge_iterator end_local_cross_edges()
	{
		return get_edge_iterator(_end_cross_edges);//edge_iterator(_end_cross_edges, edges);
	}

	const_edge_iterator end_local_cross_edges() const
	{
		return get_edge_iterator(_end_cross_edges);//edge_iterator(_end_cross_edges, edges);
	}


	/*!
	 * \brief Deactivates the local edge, that is referenced by the given iterator.
	 *
	 * The given edge iterator will reference a new active edge after this operation,
	 * unless all edges are inactive.
	 *
	 * The new referenced edge will be an edge that hasn't been iterated yet!!!
	 * (unless all edges are inactive).
	 *
	 * @param edge_it	An iterator to the local edge to be deactivated.
	 *
	 * @return			Returns an iterator-reference to the deactivated local edge.
	 */
	edge_iterator deactivate_local_edge(const edge_iterator edge_it);


	/*!
	 * \brief Deactivates the cross edge, that is referenced by the given iterator.
	 *
	 * The given edge iterator will reference a new active edge after this operation,
	 * unless all edges are inactive.
	 *
	 * The new referenced edge will be an edge that hasn't been iterated yet!!!
	 * (unless all edges are inactive).
	 *
	 * @param edge_it	An iterator to the cross edge to be deactivated.
	 *
	 * @return			Returns an iterator-reference to the deactivated cross edge.
	 */
	edge_iterator deactivate_local_cross_edge(const edge_iterator edge_it);


	bool correct_state() const;


	/*!
	 * \brief Returns true if there are active edges, otherwise false.
	 *
	 * @return Returns true if there are active edges, otherwise false.
	 */
	bool has_active_edges() const
	{
//		std::cout << "proc " << proc_id << ": " << _end_active_local_edges << ", " << _end_active_local_cross_edges << std::endl;
		return _end_active_local_edges | (_end_active_local_cross_edges-_start_cross_edges);
	}

	/*!
	 * \brief Returns the number of active edges.
	 *
	 * @return Returns the number of active edges.
	 */
	EdgeIDType num_active_edges() const
	{
		return _end_active_local_edges + (_end_active_local_cross_edges-_start_cross_edges);
	}

	WeightType get_weight() const;


	// if proc_id == self_id -> return false
	bool is_active_partner(const int proc_id) const
	{
		return active_cross_edges_of_proc[proc_id];
	}

	EdgeIDType get_active_cross_edges_count_of_proc(const int proc_id) const
	{
		return active_cross_edges_of_proc[proc_id];
	}


	bool global_id_is_local(const NodeIDType n_id) const
	{
		return n_id>=first_global_vertex && n_id<end_global_vertex;
	}

	/*!
	 *
	 * @param n_id local id
	 * @return
	 */
	bool is_local(const NodeIDType n_id) const
	{
		return n_id<_num_local_vertices;
	}

	/*!
	 *
	 * @param n_id local id
	 * @return
	 */
	bool is_ghost(const NodeIDType n_id) const
	{
		return n_id>=_num_local_vertices && n_id<_num_all_local_vertices;
	}


	NodeIDType global_vertex_id_to_local_id(const NodeIDType n_id) const
	{
		if(global_id_is_local(n_id))
		{
			return global_vertex_id_to_local_id_of_local_vertex(n_id);
		}
		else
		{
			return global_vertex_id_to_local_id_of_ghost_vertex(n_id);
		}
	}

	NodeIDType global_vertex_id_to_local_id_of_local_vertex(const NodeIDType n_id) const
	{
		return n_id-first_global_vertex;
	}


	/*!
	 * \brief Returns the local ID of a (local) ghost vertex, identified by its
	 * global ID.
	 *
	 * Runtime: O(log(number of local ghost vertices))
	 *
	 *	\todo Implement a second version of this function, but it get's the ID
	 *	of an adjacent vertex as a hint => thus should be able to implement the
	 *	function in O(vertex degree) or even O(log(vertex degree))
	 *
	 * @param n_id	The global ID of the ghost vertex
	 * @return		The local ID of the ghost vertex
	 */
	NodeIDType global_vertex_id_to_local_id_of_ghost_vertex(const NodeIDType n_id) const
	{
//		// do binary search
//		NodeIDType first = 0;
//		NodeIDType last  = _num_local_ghost_vertices;
//
//		while(last-first>1) // !!! >1
//		{
//			if(n_id < ghost_global_to_local[(first+last)/2].first)
//			{
//				last = (first+last)/2;
//			}
//			else
//			{
//				first = (first+last)/2;
//			}
//		}
//
//		// what if last==first?? -> id doesn't exist
//
//		return ghost_global_to_local[first].second;

		return ghost_global_to_local_hash.at(n_id);
	}

	// n_id must be the local id of a ghost vertex
	int get_proc_of_ghost_vertex(const NodeIDType n_id) const
	{
		return proc_of_ghost_vertex[n_id-_num_local_vertices];
	}

	/*!
	 *
	 * if n_id corresponds to the dummy node the value 0 is returned.
	 *
	 * @param n_id
	 * @return
	 */
	NodeIDType local_vertex_id_to_global_id(const NodeIDType n_id) const
	{
		return local_to_global[n_id];
	}

	int get_active_partner_count() const
	{
		return active_partners;
	}



	void get_local_conntected_components(ConnectedComponents &connected_components);


	void print() const
	{
		for(EdgeIDType pos=0; pos<edges.size(); pos++)
		{
			std::cout << edges[pos] << ",  ";
		}
		std::cout << std::endl;
	}

	bool is_maximal_matching(std::vector<Edge> &matching);


	void reserve(unsigned long n)
	{
		edges.reserve(n+1);
	}

	void push_back(const Edge &e)
	{
		edges.push_back(e);
	}

	NodeIDType size() const
	{
		// we only have a dummy edge if edges.size>0 and if the graph
		// has been initialized
		return edges.size()-((edges.size()>0) && initialized);
	}

	Edge& operator[](const NodeIDType &n)
	{
		return edges[n];
	}

	Edge operator[](const NodeIDType &n) const
	{
		return edges[n];
	}

	void initialize(const NodeIDType first_global_vertex,
			const NodeIDType end_global_vertex,
			const std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges);


	void activate_edges()
	{
		_end_active_local_edges = _num_local_edges;
		_end_active_local_cross_edges = _num_all_local_edges;

		// count the active cross edges for each proc
		active_cross_edges_of_proc.assign(procs, 0);

		for(EdgeIDType e_id=_start_cross_edges; e_id<_end_cross_edges; e_id++)
		{
			// .n2 is ghost node
			active_cross_edges_of_proc[get_proc_of_ghost_vertex(edges[e_id].n2)]++;
		}

		active_partners = 0;
		for(int p=0; p<procs; p++)
		{
			if(active_cross_edges_of_proc[p]>0)
			{
				active_partners++;
			}
		}
	}

private:

	static bool pair_comp(const std::pair<NodeIDType, NodeIDType> fst, const std::pair<NodeIDType, NodeIDType> snd)
	{
		return fst.first < snd.first;
	}

	/*!
	 * \brief returns an edge iterator pointing to the edge starting at
	 * position edge_pos.
	 *
	 * @param edge_pos	The position of the edge
	 * @return			Edge iterator pointing to the specified edge.
	 */
	const_edge_iterator get_edge_iterator(const EdgeIDType &edge_pos) const
	{
		return edges.begin()+edge_pos;
	}

	edge_iterator get_edge_iterator(const EdgeIDType &edge_pos)
	{
		return edges.begin()+edge_pos;
	}

	/// \todo make private
//private:
public:

	NodeIDType first_global_vertex;
	NodeIDType end_global_vertex;

	NodeIDType _num_local_vertices;
	NodeIDType _num_local_ghost_vertices;
	NodeIDType _num_all_local_vertices;

	EdgeIDType _num_local_edges;
	EdgeIDType _num_local_cross_edges;
	EdgeIDType _num_all_local_edges;

	EdgeIDType _end_active_local_edges; // should point to the first position of inactive edges
	EdgeIDType _end_active_local_cross_edges;

	EdgeIDType _start_cross_edges, _end_cross_edges;

//	std::vector<Edge> local_edges;
//	std::vector<Edge> local_cross_edges;

	// first part local edges, second part cross edges
	std::vector<Edge> edges;

	std::vector<int> proc_of_ghost_vertex;

	std::vector<EdgeIDType> active_cross_edges_of_proc;

	// first element of the pair represents the global id of the node and the second element the local id
//	std::vector< std::pair<NodeIDType, NodeIDType> > ghost_global_to_local;

	// stores the global id for each local vertex
	// also add one dummy node
	std::vector<NodeIDType> local_to_global;

	boost::unordered_map<NodeIDType, NodeIDType> ghost_global_to_local_hash;

	int active_partners;

	int proc_id, procs;

	bool initialized;

	Edge dummy_edge;
	NodeIDType dummy_vertex;
};


//template <typename W, typename N,typename E>
//class ParallelEdgeGraph<W, N, E>::edge_iterator
//{
//	typedef typename ParallelEdgeGraph<W, N, E>::edge_iterator self_type;
//
//private:
//	// index/refrence to an edge within the edges-vector
//	EdgeIDType edge_ref;
//
//	const std::vector<Edge> *edge_vector;
//
//public:
//
//	edge_iterator(const EdgeIDType &ref, const std::vector<Edge> &edge_vector)
//	: edge_ref(ref), edge_vector(&edge_vector)
//	{}
//
//	self_type operator++(int )
//	{
//		self_type old = *this;
//		edge_ref++;
//		return old;
//	}
//
//	self_type& operator++()
//	{
//		++edge_ref;
//		return *this;
//	}
//
//	bool operator==(const self_type &r) const
//	{
//		return edge_ref==r.edge_ref;
//	}
//
//	bool operator!=(const self_type &r) const
//	{
//		return edge_ref!=r.edge_ref;
//	}
//
//	bool operator<(const self_type &r) const
//	{
//		return edge_ref<r.edge_ref;
//	}
//
//	bool operator<=(const self_type &r) const
//	{
//		return edge_ref<=r.edge_ref;
//	}
//
//	bool operator>(const self_type &r) const
//	{
//		return edge_ref>r.edge_ref;
//	}
//
//	bool operator>=(const self_type &r) const
//	{
//		return edge_ref>=r.edge_ref;
//	}
//
//	const Edge& operator*() const
//	{
//		return (*edge_vector)[edge_ref];
//	}
//
//	const Edge* operator->() const
//	{
//		return &(*edge_vector)[edge_ref];
//	}
//
//	friend class ParallelEdgeGraph<W, N, E>;
//};



template <typename W, typename N, typename E>
ParallelEdgeGraph<W, N, E>::ParallelEdgeGraph(const NodeIDType global_node_count, const std::vector<Edge> &local_and_cross_edges,
											  const NodeIDType first_global_vertex, const NodeIDType end_global_vertex,
											  const std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges,
											  const bool delete_provided_edges)
: first_global_vertex(first_global_vertex), end_global_vertex(end_global_vertex),
  initialized(false), dummy_vertex(std::numeric_limits<NodeIDType>::max())
{
	// add edges to graph
	edges.reserve(local_and_cross_edges.size()+1); // don't forget the dummy edge

	for(typename std::vector<Edge>::const_iterator e_it=local_and_cross_edges.begin(); e_it!=local_and_cross_edges.end(); e_it++)
	{
		edges.push_back(*e_it);
	}

	if(delete_provided_edges)
	{// delete the given edge vector
		delete &local_and_cross_edges;
	}

	initialize(first_global_vertex, end_global_vertex, proc_ranges);

	dummy_edge.n1 = std::numeric_limits<NodeIDType>::max();
	dummy_edge.n2 = std::numeric_limits<NodeIDType>::max();
	dummy_edge.weight = min_val<WeightType>();

	return;
}




template <typename W, typename N, typename E>
typename ParallelEdgeGraph<W, N, E>::edge_iterator ParallelEdgeGraph<W, N, E>::deactivate_local_edge(const edge_iterator edge_it)
{
	EdgeIDType edge_ref = get_edge_id(edge_it);

	if(edge_ref>=_end_active_local_edges)
	{
		return edge_it;
	}

	_end_active_local_edges--;

	Edge swap_edge = edges[_end_active_local_edges];
	edges[_end_active_local_edges] = edges[edge_ref];
	edges[edge_ref] = swap_edge;


	return get_edge_iterator(_end_active_local_edges);//  (_end_active_local_edges, edges);
}


template <typename W, typename N, typename E>
typename ParallelEdgeGraph<W, N, E>::edge_iterator ParallelEdgeGraph<W, N, E>::deactivate_local_cross_edge(const edge_iterator edge_it)
{
	EdgeIDType edge_ref = get_edge_id(edge_it);

	if(edge_ref>=_end_active_local_cross_edges)
	{
		return edge_it;
	}

	// .n2 is ghost - we ensure this in the constructor
	active_cross_edges_of_proc[get_proc_of_ghost_vertex(edges[edge_ref].n2)]--;

	active_partners -= (active_cross_edges_of_proc[get_proc_of_ghost_vertex(edges[edge_ref].n2)]==0);

	_end_active_local_cross_edges--;

	Edge swap_edge = edges[_end_active_local_cross_edges];
	edges[_end_active_local_cross_edges] = edges[edge_ref];
	edges[edge_ref] = swap_edge;


	return get_edge_iterator(_end_active_local_cross_edges);// edge_iterator(_end_active_local_cross_edges, edges);
}




// stores local vertex IDs in connected_components
template <typename W, typename N, typename E>
void ParallelEdgeGraph<W, N, E>::get_local_conntected_components(ConnectedComponents &connected_components)
{
	NodeIDType num_vertices = num_all_local_vertices();

	UnionFind<NodeIDType> uf(num_vertices);


	for(EdgeIDType e_id=0; e_id<_num_all_local_edges; e_id++)
	{
		uf.combine(edges[e_id].n1, edges[e_id].n2);
	}

	std::vector<ConnectedComponent> sets(num_vertices);

	for(NodeIDType n=0; n<num_vertices; n++)
	{
		sets[uf.find(n)].push_back(n);
	}


	for(NodeIDType n=0; n<num_vertices; n++)
	{
		if(! sets[n].empty())
		{
			connected_components.push_back(sets[n]);
		}
	}

}


template <typename W, typename N, typename E>
bool ParallelEdgeGraph<W, N, E>::is_maximal_matching(std::vector<Edge> &matching)
{
	int procs, proc_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	int result = true;

	std::vector<bool> matched(num_all_local_vertices(), false);
	std::vector< std::vector<Edge> > msgs_for_proc(procs);
	std::vector<MPI_Request> req_for_proc(procs, MPI_REQUEST_NULL);

	// set matched nodes
//	for(EdgeIDType i=0; i<(EdgeIDType)matching.size(); i++)
	for(typename std::vector<Edge>::iterator e_it=matching.begin(); e_it!=matching.end(); e_it++)
	{
		const Edge &e = *e_it;//matching[i];
		NodeIDType local_n1 = global_vertex_id_to_local_id(e.n1);
		NodeIDType local_n2 = global_vertex_id_to_local_id(e.n2);

		if(matched[local_n1])
		{
			std::cout << "proc " << proc_id << ": Vertex " << e.n1 << " incident to more than one matched edge!" << std::endl;
			result = false;
		}

		if(matched[local_n2])
		{
			std::cout << "proc " << proc_id << ": Vertex " << e.n2 << " incident to more than one matched edge!" << std::endl;
			result = false;
		}

		matched[local_n1] = true;
		matched[local_n2] = true;

		// check for cross-edge -> send cross edge to partner
		// only n2 might be a ghost node
		if(is_ghost(local_n2))
		{
			msgs_for_proc[get_proc_of_ghost_vertex(local_n2)].push_back(e);
		}
	}

	int test_for_matching_tag = 13;

	MPI_Datatype EDGE_TYPE;
	MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &EDGE_TYPE);
	MPI_Type_commit(&EDGE_TYPE);

	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			if(!msgs_for_proc[p].empty())
			{
				MPI_Isend(&msgs_for_proc[p][0], msgs_for_proc[p].size(), EDGE_TYPE, p, test_for_matching_tag, MPI_COMM_WORLD, &req_for_proc[p]);
			}
			else
			{
				MPI_Isend(0, 0, EDGE_TYPE, p, test_for_matching_tag, MPI_COMM_WORLD, &req_for_proc[p]);
			}
		}
	}

	std::vector<Edge> matched_edges_of_partners;

	for(int p=0; p<procs; p++)
	{
		if(p!=proc_id)
		{
			MPI_Status status;
			MPI_Probe(p, test_for_matching_tag, MPI_COMM_WORLD, &status);

			int count;
			MPI_Get_count(&status, EDGE_TYPE, &count);

			if(count>0)
			{
				std::vector<Edge> recv_edges(count);
				MPI_Recv(&recv_edges[0], count, EDGE_TYPE, p, test_for_matching_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				for(EdgeIDType i=0; i<(EdgeIDType)count; i++)
				{
					matched_edges_of_partners.push_back(recv_edges[i]);
				}
			}
		}
	}



	for(EdgeIDType i=0; i<(EdgeIDType)matched_edges_of_partners.size(); i++)
	{
		const Edge &e = matched_edges_of_partners[i];

		NodeIDType local_n1 = global_vertex_id_to_local_id(e.n1);
		NodeIDType local_n2 = global_vertex_id_to_local_id(e.n2);

		if(matched[local_n1])
		{
			std::cout << "proc " << proc_id << ": Vertex " << e.n1 << " incident to more than one matched edge!" << std::endl;
			result = false;
		}

		if(matched[local_n2])
		{
			std::cout << "proc " << proc_id << ": Vertex " << e.n2 << " incident to more than one matched edge!" << std::endl;
			result = false;
		}

		matched[local_n1] = true;
		matched[local_n2] = true;
	}


	// check for edges not incident to matched nodes
	for(edge_iterator e_it=begin_active_local_edges();
			e_it!=end_active_local_edges(); e_it++)
	{
		const Edge &e = *e_it;

		if(!matched[e.n1] && !matched[e.n2])
		{
			std::cout << "proc " << proc_id << ": Edge " << e << " not matched, although non of the incident edges is matched!" << std::endl;
			result = false;
		}
	}

	for(edge_iterator e_it=begin_active_local_cross_edges();
			e_it!=end_active_local_cross_edges(); e_it++)
	{
		const Edge &e = *e_it;

		if(!matched[e.n1] && !matched[e.n2])
		{
			std::cout << "proc " << proc_id << ": Edge " << e << " not matched, although non of the incident edges is matched!" << std::endl;
			result = false;
		}
	}

	MPI_Waitall(procs, &req_for_proc[0], MPI_STATUSES_IGNORE);

	int global_result;

	MPI_Allreduce(&result, &global_result, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

	return global_result;
}


template <typename W, typename N, typename E>
void ParallelEdgeGraph<W, N, E>::initialize(const NodeIDType _first_global_vertex, const NodeIDType _end_global_vertex,
											const std::vector<std::pair<NodeIDType, NodeIDType> > &proc_ranges)
{
	if(initialized)
	{// we only allow to initialize the graph once
		return;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	first_global_vertex = _first_global_vertex;
	end_global_vertex = _end_global_vertex;

	//	// there are some problems if NodeIDType is a floating point type
	_num_local_vertices = end_global_vertex-first_global_vertex;

	std::vector<NodeIDType> ghost_vertices; // might contain multiple entries of the same ghost vertex


	EdgeIDType num_local_edges = 0;
	EdgeIDType num_cross_edges = 0; // also corresponds to the size of ghost_vertices


	// first run to get number of local edges and cross edges -> for better memory usage
	for(typename std::vector<Edge>::const_iterator e_it=edges.begin(); e_it!=edges.end(); e_it++)
	{

		if(global_id_is_local(e_it->n1) && global_id_is_local(e_it->n2))
		{
			num_local_edges++;
		}
		else if(!global_id_is_local(e_it->n1) && global_id_is_local(e_it->n2))
		{
			num_cross_edges++;
		}
		else if(global_id_is_local(e_it->n1) && !global_id_is_local(e_it->n2))
		{
			num_cross_edges++;
		}
	}

	ghost_vertices.reserve(num_cross_edges);

	_num_local_edges       = num_local_edges;
	_num_local_cross_edges = num_cross_edges;
	_num_all_local_edges   = _num_local_edges + _num_local_cross_edges;

	_start_cross_edges = _num_local_edges;
	_end_cross_edges = _num_all_local_edges;

	_end_active_local_edges       = _num_local_edges;
	_end_active_local_cross_edges = _num_all_local_edges;

	// we add a dummy edge
	edges.resize(_num_all_local_edges+1);

	// set dummy edge
	edges[_num_all_local_edges].n1 = 0;
	edges[_num_all_local_edges].n2 = 0;
	edges[_num_all_local_edges].weight = min_val<WeightType>();

	EdgeIDType cross_edge_pos = _num_local_edges;


	// rearrange edges, at first local edges, then cross edges
	for(NodeIDType n=0; n<_num_local_edges; )
	{
		const Edge &e = edges[n];

		if(global_id_is_local(e.n1) && global_id_is_local(e.n2))
		{// local edge is in correct part
			n++;
		}
		else
		{// move cross edge to second part
			swap(n, cross_edge_pos, edges);
			cross_edge_pos++;
		}
	}

	//make sure the cross edges have the correct structure
	for(NodeIDType n=_num_local_edges; n<edges.size()-1; n++) // we don't care about the dummy edge
	{
		if(!global_id_is_local(edges[n].n1))
		{// first vertex is not supposed to be the ghost vertex
			// the first node of cross edges is supposed to be local!!!
			NodeIDType ghost_vertex = edges[n].n1;
			edges[n].n1 = edges[n].n2;
			edges[n].n2 = ghost_vertex;

			ghost_vertices.push_back(ghost_vertex);
		}
		else
		{// don't forget to add the ghost vertex
			ghost_vertices.push_back(edges[n].n2);
		}
	}

	// sort ghost vertices, hopefully it's inplace
	std::sort(ghost_vertices.begin(), ghost_vertices.end());

	NodeIDType new_local_ghost_id = _num_local_vertices; // the first local id of the local ghost vertices starts with the number _num_local_vertices
	NodeIDType last_global_id = -1;

//	NodeIDType num_ghost_vertices = 0;
//
//	// get number of ghost vertices
//	for(NodeIDType i=0; i<ghost_vertices.size(); i++)
//	{
//		if(last_global_id!=ghost_vertices[i])
//		{
//			last_global_id = ghost_vertices[i];
//			num_ghost_vertices++;
//		}
//	}

	//	ghost_global_to_local.reserve(num_ghost_vertices);

	// we have to reset last_global_id or we might miss a few nodes!!!
	last_global_id = -1;

	_num_local_ghost_vertices = 0;


	// fill ghost_global_to_local vector
	for(NodeIDType i=0; i<ghost_vertices.size(); i++)
	{
		if(last_global_id!=ghost_vertices[i])
		{
			last_global_id = ghost_vertices[i];

			//			ghost_global_to_local.push_back(std::make_pair(last_global_id, new_local_ghost_id));

			ghost_global_to_local_hash.insert(std::make_pair(last_global_id, new_local_ghost_id));

			_num_local_ghost_vertices++;

			new_local_ghost_id++;
		}
	}


	//	std::sort(ghost_global_to_local.begin(), ghost_global_to_local.end(), pair_comp);


	//	_num_local_ghost_vertices = ghost_global_to_local.size();
	_num_all_local_vertices = _num_local_vertices + _num_local_ghost_vertices;

	// add one dummy node
	local_to_global.resize(_num_all_local_vertices+1);
	local_to_global[_num_all_local_vertices] = 0; // set the dummy node some value that's not totally weird

	// set the global IDs of the local vertices
	for(NodeIDType i=0; i< _num_local_vertices; i++)
	{
		local_to_global[i] = i+first_global_vertex;
	}

	// set the global IDs of the ghost vertices
	for(typename boost::unordered_map<NodeIDType, NodeIDType>::iterator it=ghost_global_to_local_hash.begin();
			it!=ghost_global_to_local_hash.end(); it++)
	{
		// first value of the unordered_map::iterator is the key (global id)
		// and the second is the value (local id)
		local_to_global[it->second] = it->first;
	}


	// rename the IDs of the local edges
	for(EdgeIDType i=0; i<_num_local_edges; i++)
	{
		edges[i].n1 = global_vertex_id_to_local_id_of_local_vertex(edges[i].n1);
		edges[i].n2 = global_vertex_id_to_local_id_of_local_vertex(edges[i].n2);
	}

	// rename the IDs of the cross edges, they start right after the local edges
	for(EdgeIDType i=_start_cross_edges; i<_end_cross_edges; i++)
	{
		edges[i].n1 = global_vertex_id_to_local_id(edges[i].n1);
		edges[i].n2 = global_vertex_id_to_local_id(edges[i].n2);
	}


	proc_of_ghost_vertex.resize(_num_local_ghost_vertices);


	// set the proc for each ghost vertex
	for(NodeIDType g_id=_num_local_vertices; g_id<_num_all_local_vertices; g_id++)
	{
		for(unsigned int p=0; p<proc_ranges.size(); p++)
		{
			if(proc_ranges[p].first<=local_to_global[g_id] && local_to_global[g_id]<proc_ranges[p].second)
			{
				proc_of_ghost_vertex[g_id-_num_local_vertices] = p;
			}
		}
	}


	// count the active cross edges for each proc
	active_cross_edges_of_proc.resize(procs, 0);

	for(EdgeIDType e_id=_start_cross_edges; e_id<_end_cross_edges; e_id++)
	{
		// .n2 is ghost node
		active_cross_edges_of_proc[get_proc_of_ghost_vertex(edges[e_id].n2)]++;
	}

	active_partners = 0;
	for(int p=0; p<procs; p++)
	{
		if(active_cross_edges_of_proc[p]>0)
		{
			active_partners++;
		}
	}

	initialized = true;
}


} // end namespace dipl


#endif













