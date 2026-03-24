#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>

#include <common/linear_buffered_reader.h>
#include <common/edge.h>

namespace dipl
{

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1 nodes
/// _EdgeIDType should be at least unsigned long if there are more than 2^32 - 1 edges
template <typename _WeightType=double,
		  typename _NodeIDType=unsigned int,
		  typename _EdgeIDType=unsigned int>
class EdgeGraph
{

public:

	typedef _WeightType WeightType;
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

	typedef dipl::Edge<WeightType, NodeIDType> 					 EdgeWithoutID;
	typedef dipl::EdgeWithID<WeightType, NodeIDType, EdgeIDType> Edge;

//	class edge_iterator;
	typedef typename std::vector<Edge>::iterator 		edge_iterator;
	typedef typename std::vector<Edge>::const_iterator 	const_edge_iterator;

public:

	/*!
	 * \brief Constructor that initializes the graph from the given edge list.
	 *
	 * Runtime: O(m+n)
	 *
	 * @param node_count	Number of vertices of the graph, defined by the edge list.
	 * @param edge_list		List of edges that define the graph. The vertices of the
	 * 						graph have to be numbered from 0 to node_count-1
	 */
	EdgeGraph(const NodeIDType node_count, const std::list<EdgeWithoutID > &edge_list);

	/*!
	 * \brief Constructor that initializes the graph from the given edge list.
	 *
	 * Runtime: O(m+n)
	 *
	 * @param node_count	Number of vertices of the graph, defined by the edge list.
	 * @param edge_list		List of edges that define the graph. The vertices of the
	 * 						graph have to be numbered from 0 to node_count-1
	 */
	EdgeGraph(const NodeIDType node_count, const std::vector<EdgeWithoutID > &edge_vector);



	/*!
	 * \brief Returns the number of vertices of the graph.
	 *
	 * Runtime: O(1)
	 *
	 * @return Number of vertices of the graph.
	 */
	NodeIDType num_vertices() const
	{
		return _num_vertices;
	}


	/*!
	 * \brief Returns the number of edges of the graph.
	 *
	 * Runtime: O(1)
	 *
	 * @return Number of edges of the Graph
	 */
	EdgeIDType num_edges() const
	{
		return _num_edges;
	}

	/*!
	 * \brief Returns the edge with the given edge_pos.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_ref
	 * @return			Returns a const reference to the edge at the given pos.
	 */
	const Edge& get_edge(const EdgeIDType &edge_pos) const
	{
		return edges[edge_pos];
	}


	/*!
	 * \brief Returns the edge with the given edge_pos.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_ref
	 * @return			Returns a reference to the edge at the given pos.
	 */
	Edge& get_edge(const EdgeIDType &edge_pos)
	{
		return edges[edge_pos];
	}

	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	const Edge& get_edge(const_edge_iterator &edge_it) const
	{
		return *edge_it;//get_edge(edge_management[edge_it.edge_ref]);
	}

	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	Edge& get_edge(edge_iterator &edge_it)
	{
		return *edge_it;//edges[edge_it.edge_ref];
	}


	/*!
	 * \brief Returns the edge id corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	EdgeIDType get_edge_id(const_edge_iterator &edge_it) const
	{
		return edge_it->edge_id;
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	EdgeIDType get_edge_id(edge_iterator &edge_it) const
	{
		return edge_it->edge_id;
	}

	/*!
	 * \brief Returns the current position of the edge referenced
	 * by the given iterator.
	 *
	 * @param edge_it
	 * @return
	 */
	EdgeIDType get_edge_pos(const_edge_iterator &edge_it) const
	{
		return edge_it-edges.begin();//edge_it.edge_ref;
	}

	EdgeIDType get_edge_pos(edge_iterator &edge_it) const
	{
		return edge_it-edges.begin();//edge_it.edge_ref;
	}



	/*!
	 * \brief Returns the weight of the edge referenced by the given iterator.
	 *
	 * Don't use edge references!!!
	 *
	 * @param edge_it
	 * @return			Returns the weight of the edge referenced by the given iterator.
	 */
	WeightType get_edge_weight(const_edge_iterator &edge_it) const
	{
		return edge_it->weight;
	}

	WeightType get_edge_weight(edge_iterator &edge_it) const
	{
		return edge_it->weight;
	}


	/*!
	 * \brief Returns an edge_id_iterator that points to the first edge of the graph.
	 *
	 * Runtime: O(1)
	 *
	 * @return	edge_id_iterator of the first edge of the graph.
	 */
	edge_iterator begin_edges()
	{
		return edges.begin();//edge_iterator(0, this);
	}

	const_edge_iterator begin_edges() const
	{
		return edges.begin();//edge_iterator(0, this);
	}


	/*!
	 * \brief Returns an edge_id_iterator that indicates the end.
	 *
	 * Runtime: O(1)
	 *
	 * @return	edge_id_iterator of the end.
	 */
	edge_iterator end_edges()
	{
		return get_edge_iterator(_num_edges);//edge_iterator(_num_edges, this);
	}

	const_edge_iterator end_edges() const
	{
		return get_edge_iterator(_num_edges);//edge_iterator(_num_edges, this);
	}



	edge_iterator begin_active_edges()
	{
		return edges.begin();//edge_iterator(0, this);
	}

	const_edge_iterator begin_active_edges() const
	{
		return edges.begin();//edge_iterator(0, this);
	}



	edge_iterator end_active_edges()
	{
		return get_edge_iterator(_end_active_edges);//edge_iterator(_end_active_edges, this);
	}

	const_edge_iterator end_active_edges() const
	{
		return get_edge_iterator(_end_active_edges);//edge_iterator(_end_active_edges, this);
	}



	/*!
	 * \brief Returns whether the given edge is inactive.
	 *
	 * Don't use edge references!!!!
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id	The edge id.
	 * @return			Returns true if the given edge is inactive, otherwise false.
	 */
	inline bool is_inactive(EdgeIDType edge_id) const;


	/*!
	 * \brief Returns whether the given edge is active.
	 *
	 * Don't use edge references!!!!
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id	The edge id.
	 * @return			Returns true if the given edge is active, otherwise false.
	 */
	inline bool is_active(EdgeIDType edge_id) const;

	/*!
	 * \brief Deactivates the edge, that is referenced by the given iterator.
	 *
	 * The given edge iterator will reference a new active edge after this operation,
	 * unless all edges are inactive.
	 *
	 * The new referenced edge will be an edge that hasn't been iterated yet!!!
	 * (unless all edges are inactive).
	 *
	 * @param edge_it	An iterator to the edge to be deactivated.
	 *
	 * @return			Returns an iterator-reference to the deactivated edge.
	 */
	edge_iterator deactivate_edge(edge_iterator &edge_it)
	{
		return get_edge_iterator(deactivate_edge(get_edge_pos(edge_it)));
//				edge_iterator(deactivate_edge(edge_it.edge_ref), this);
	}

	const_edge_iterator deactivate_edge(const_edge_iterator &edge_it)
	{
		return get_edge_iterator(deactivate_edge(get_edge_pos(edge_it)));
//				edge_iterator(deactivate_edge(edge_it.edge_ref), this);
	}

	bool correct_state() const;


	/*!
	 * \brief Returns true if there are active edges, otherwise false.
	 *
	 * @return Returns true if there are active edges, otherwise false.
	 */
	bool has_active_edges() const
	{
		return _end_active_edges; // as long _end_active_edges != 0 there are active edges, so it's just a cast to bool
	}

	/*!
	 * \brief Returns the number of active edges.
	 *
	 * @return Returns the number of active edges.
	 */
	EdgeIDType num_active_edges() const
	{
		return _end_active_edges;
	}

	WeightType get_weight() const;

	void activate_edges()
	{
		_end_active_edges = _num_edges;

		for(EdgeIDType i=0; i<_num_edges; i++)
		{
			active_edges[i] = true;
		}
	}


	void activate_and_shuffle_edges()
	{
		activate_edges();
		_end_active_edges = _num_edges;

		std::random_shuffle(edges.begin(),edges.begin()+=_num_edges);
	}


	static bool edge_comp(const Edge &fst, const Edge &snd)
	{
		return (fst.n1 < snd.n1) || ((fst.n1 == snd.n1) && (fst.n2 < snd.n2));
	}

	void activate_and_sort_edges()
	{
		activate_edges();
		_end_active_edges = _num_edges;

		std::sort(edges.begin(),edges.begin()+=_num_edges, edge_comp);
	}

private:

	/*!
	 * \brief Deactivates the edge, that is referenced by edge_ref.
	 *
	 * The referenced position stores afterwards an active edge, unless
	 * all edges are inactive.
	 *
	 * The new referenced edge comes from a position >= edge_ref.
	 *
	 * @param edge_ref	A reference to the edge to be deactivated.
	 *
	 * @return			Returns a reference to the deactivated edge.
	 */
	EdgeIDType deactivate_edge(const EdgeIDType edge_ref);

	/*!
	 * \brief returns an edge iterator pointing to the edge starting at
	 * position edge_pos.
	 *
	 * @param edge_pos	The position of the edge
	 * @return			Edge iterator pointing to the specified edge.
	 */
	edge_iterator get_edge_iterator(const EdgeIDType &edge_pos)
	{
		return edges.begin()+edge_pos;
	}

private:
	NodeIDType _num_vertices;
	EdgeIDType _num_edges;

	EdgeIDType _end_active_edges; // should point to the first position of inactive edges

	std::vector<Edge> edges;	// stores edges, the first part contains of active
								// edges and the second part (starting at _end_active_edges
								// contains the inactive edges

	std::vector<bool> active_edges; // needs size: _num_edges+1, simulates a dummy edge that's always inactive

};


//template <typename W, typename N,typename E>
//class EdgeGraph<W, N, E>::edge_iterator
//{
//	typedef typename EdgeGraph<W, N, E>::edge_iterator self_type;
//
//private:
//	// index/refrence to an edge within the edges-vector
//	EdgeIDType edge_ref;
//
//	const EdgeGraph<W, N, E> *edge_graph;
//
//public:
//
//	edge_iterator(const EdgeIDType &ref, const EdgeGraph<W, N, E> *edge_graph)
//	: edge_ref(ref), edge_graph(edge_graph)
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
//		return !(*this==r);
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
//		return edge_graph->edges[edge_ref];
//	}
//
//	const Edge* operator->() const
//	{
//		return &(edge_graph->edges[edge_ref]);
//	}
//
//	friend class EdgeGraph<W, N, E>;
//};



template <typename W, typename N, typename E>
EdgeGraph<W, N, E>::EdgeGraph(const NodeIDType node_count, const std::list<EdgeWithoutID> &edge_list)
: _num_vertices(node_count), _num_edges(edge_list.size()), _end_active_edges(edge_list.size()),
  active_edges(_num_edges+1, true)
{
	edges.reserve(_num_edges+1);

	// set the edge weights
	EdgeIDType current_id = 0;
	for(typename std::list<dipl::Edge<WeightType, NodeIDType> >::const_iterator it=edge_list.begin(); it!=edge_list.end(); it++, current_id++)
	{
		edges.push_back(Edge(it->n1, it->n2, edges.size(), it->weight));
	}

	// add dummy_edge
	edges.push_back(Edge(-1, -1, _num_edges, min_val<WeightType>()));

	active_edges[_num_edges] = false; //the last dummy edges is always inactive!!!

}


template <typename W, typename N, typename E>
EdgeGraph<W, N, E>::EdgeGraph(const NodeIDType node_count, const std::vector<EdgeWithoutID> &edge_vector)
: _num_vertices(node_count), _num_edges(edge_vector.size()), _end_active_edges(edge_vector.size()),
  active_edges(_num_edges+1, true)
{
	edges.reserve(_num_edges+1);

	// set the edge weights
	EdgeIDType current_id = 0;
	for(typename std::vector<dipl::Edge<WeightType, NodeIDType> >::const_iterator it=edge_vector.begin(); it!=edge_vector.end(); it++, current_id++)
	{
		edges.push_back(Edge(it->n1, it->n2, edges.size(), it->weight));
	}

	// add dummy_edge
	edges.push_back(Edge(-1, -1, _num_edges, min_val<WeightType>()));

	active_edges[_num_edges] = false; //the last dummy edges is always inactive!!!

}



template <typename W, typename N, typename E>
inline bool EdgeGraph<W, N, E>::is_inactive(const EdgeIDType edge_id) const
{
	return !active_edges[edge_id];
}


template <typename W, typename N, typename E>
inline bool EdgeGraph<W, N, E>::is_active(const EdgeIDType edge_id) const
{
	return active_edges[edge_id];
}


template <typename W, typename N, typename E>
typename EdgeGraph<W, N, E>::EdgeIDType EdgeGraph<W, N, E>::deactivate_edge(const EdgeIDType edge_ref)
{
	if(is_inactive(edges[edge_ref].edge_id))
	{
		return edge_ref;
	}

	_end_active_edges--;

	Edge swap_edge = edges[_end_active_edges];
	edges[_end_active_edges] = edges[edge_ref];
	edges[edge_ref] = swap_edge;

	active_edges[edges[_end_active_edges].edge_id] = false;

	return _end_active_edges;
}


template <typename W, typename N, typename E>
bool EdgeGraph<W, N, E>::correct_state() const
{
	for(EdgeIDType e_ref=0; e_ref<_num_edges; e_ref++)
	{
		if( ( is_inactive(edges[e_ref].edge_id) && e_ref<_end_active_edges ) ||
			( is_active(edges[e_ref].edge_id)   && e_ref>=_end_active_edges ) )
		{
			return false;
		}
	}

	return true;
}


template <typename W, typename N, typename E>
W EdgeGraph<W, N, E>::get_weight() const
{
	W result = 0;

	for(edge_iterator e_it=begin_edges(); e_it!=end_edges(); e_it++)
	{
		result += e_it->weight;
	}

	return result;
}


} // end namespace dipl
















