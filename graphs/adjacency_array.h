#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <common/edge.h>

namespace dipl
{

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1 nodes
/// _EdgeIDType should be at least unsigned long if there are more than 2^31 - 1 edges
template <typename _WeightType=double,
		  typename _NodeIDType=unsigned int,
		  typename _EdgeIDType=unsigned int>
class AdjacencyArray
{
public:

	typedef _WeightType WeightType;
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

private:
	struct NodeStart
	{
		EdgeIDType start;
		EdgeIDType inactive_edges_start;

		NodeStart(const EdgeIDType start, EdgeIDType inactive_edges_start)
		: start(start), inactive_edges_start(inactive_edges_start)
		{}
	};

	typedef unsigned char byte;

public:
	struct Edge
	{
		NodeIDType n1;
		NodeIDType n2;

		EdgeIDType n1_edge_ref_pos;
		EdgeIDType n2_edge_ref_pos;

		WeightType weight;

		Edge()
		: n1(0), n2(0), n1_edge_ref_pos(0), n2_edge_ref_pos(0), weight(0)
		{}

		Edge(const NodeIDType n1, const NodeIDType n2, const WeightType weight)
		: n1(n1), n2(n2), n1_edge_ref_pos(0), n2_edge_ref_pos(0), weight(weight)
		{}

		friend std::ostream& operator<<(std::ostream &os, const Edge &e)
		{
			os << "(" << e.n1 << " <-> " << e.n2 << ", " << e.weight << ")";
			return os;
		}

	};

	typedef dipl::Edge<WeightType, NodeIDType> EdgeWithoutID;

	class edge_id_iterator;
	class incident_edges_iterator;

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
	AdjacencyArray(const NodeIDType node_count, const std::vector<Edge> &edge_vector);


	AdjacencyArray(const NodeIDType node_count, const std::vector<EdgeWithoutID> &edge_vector);


	virtual ~AdjacencyArray(){};

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
	 * \brief Returns the edge corresponding to the given edge_id.
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
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a refrence to the edge with the given ID.
	 */
	Edge& get_edge(const EdgeIDType &edge_id)
	{
		return edges[edge_id];
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	const Edge&  get_edge(const edge_id_iterator &edge_id) const
	{
		return get_edge(*edge_id);
	}


	/*!
	 * \brief Returns the edge corresponding to the given edge_id-iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	Edge&  get_edge(const edge_id_iterator &edge_id)
	{
		return get_edge(*edge_id);
	}


	/*!
	 * \brief Returns the edge corresponding to the given incident_edges_iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a const reference to the edge with the given ID.
	 */
	const Edge&  get_edge(const incident_edges_iterator &edge_ref) const
	{
		return get_edge(*edge_ref);
	}


	/*!
	 * \brief Returns the edge corresponding to the given incident_edges_iterator.
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id
	 * @return			Returns a reference to the edge with the given ID.
	 */
	Edge&  get_edge(const incident_edges_iterator &edge_ref)
	{
		return get_edge(*edge_ref);
	}


	/*!
	 * \brief Returns an edge_id_iterator that points to the first edge of the graph.
	 *
	 * Runtime: O(1)
	 *
	 * @return	edge_id_iterator of the first edge of the graph.
	 */
	edge_id_iterator begin_edges() const
	{
		return edge_id_iterator(0);
	}


	/*!
	 * \brief Returns an edge_id_iterator that indicates the end.
	 *
	 * Runtime: O(1)
	 *
	 * @return	edge_id_iterator of the end.
	 */
	edge_id_iterator end_edges() const
	{
		return edge_id_iterator(_num_edges);
	}


	/*!
	 * \brief Returns an incident_edges_iterator that points to the first edge of the
	 * 		incident edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator of the start of the incident edges.
	 */
	incident_edges_iterator begin_incident_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id].start, edge_refrences);
	}


	/*!
	 * \brief Returns an incident_edges_iterator that indicates the end of the incident
	 * 		edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator that indicates the end of the incident edges.
	 */
	incident_edges_iterator end_incident_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id+1].start, edge_refrences);
	}


	/*!
	 * \brief Returns an incident_edges_iterator that points to the first edge of the
	 * 		active incident edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator of the start of the active incident edges.
	 */
	incident_edges_iterator begin_incident_active_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id].start, edge_refrences);
	}


	/*!
	 * \brief Returns an incident_edges_iterator that indicates the end of the active incident
	 * 		edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator that indicates the end of the active incident edges.
	 */
	incident_edges_iterator end_incident_active_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id].inactive_edges_start, edge_refrences);
	}


	/*!
	 * \brief Returns an incident_edges_iterator that points to the first edge of the
	 * 		inactive incident edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator of the start of the inactive incident edges.
	 */
	incident_edges_iterator begin_incident_inactive_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id].inactive_edges_start, edge_refrences);
	}

	/*!
	 * \brief Returns an incident_edges_iterator that indicates the end of the inactive incident
	 * 		edges of the given node.
	 *
	 * Runtime: O(1)
	 *
	 * @param node_id	The node ID.
	 * @return			incident_edges_iterator that indicates the end of the inactive incident edges.
	 */
	incident_edges_iterator end_incident_inactive_edges(const NodeIDType &node_id) const
	{
		return incident_edges_iterator(nodes[node_id+1].start, edge_refrences);
	}


	/*!
	 * \brief Returns whether the given edge is inactive.
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
	 * Runtime: O(1)
	 *
	 * @param edge_id	The edge id.
	 * @return			Returns true if the given edge is active, otherwise false.
	 */
	inline bool is_active(EdgeIDType edge_id) const;


	/*!
	 * \brief Deactives the given edge ("adds it to the list of inactive edges");
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id	The ID of the edge to be deactivated.
	 */
	virtual void deactivate_edge(EdgeIDType edge_id);


	EdgeIDType degree(const NodeIDType node_id) const
	{
		return nodes[node_id+1].start - nodes[node_id].start;
	}

	EdgeIDType active_degree(const NodeIDType node_id) const
	{
		return nodes[node_id].inactive_edges_start - nodes[node_id].start;
	}


	void print_node_info(const NodeIDType node_id) const
	{
		std::cout << node_id << ":  start=" << nodes[node_id].start << ", inactive_edges_start=" << nodes[node_id].inactive_edges_start << "   (" << (node_id+1) << ":  start=" << nodes[node_id+1].start << ", inactive_edges_start=" << nodes[node_id+1].inactive_edges_start << ")" << std::endl;
	}


	bool correct_state() const
	{
		for(EdgeIDType e=0; e<_num_edges; e++)
		{
			if( ((edges[e].n1_edge_ref_pos >= nodes[edges[e].n1].inactive_edges_start) && (edges[e].n2_edge_ref_pos < nodes[edges[e].n2].inactive_edges_start)) ||
				((edges[e].n1_edge_ref_pos < nodes[edges[e].n1].inactive_edges_start) && (edges[e].n2_edge_ref_pos >= nodes[edges[e].n2].inactive_edges_start)) )
			{
				std::cout << "incorrect state of edge " << e << ":" << std::endl;
				std::cout << edges[e].n1_edge_ref_pos << ">=" << nodes[edges[e].n1].inactive_edges_start << " && " << edges[e].n2_edge_ref_pos << "<" << nodes[edges[e].n2].inactive_edges_start << std::endl;
				std::cout << edges[e].n1_edge_ref_pos << "<" << nodes[edges[e].n1].inactive_edges_start << " && " << edges[e].n2_edge_ref_pos << ">=" << nodes[edges[e].n2].inactive_edges_start << std::endl;

				return false;
			}
		}

		return true;
	}


	WeightType get_weight() const;


private:

	static void extract_header(std::string &line, NodeIDType &n, EdgeIDType &m, byte &fmt, unsigned int &ncon);

	static void extract_edges_of_vertex(const std::string &line, const NodeIDType src_vertex,
									    const unsigned int ncon, const bool edge_weights,
									    std::vector<Edge> &edges, EdgeIDType &vertex_degree);
private:
	NodeIDType _num_vertices;
	EdgeIDType _num_edges;


	std::vector<NodeStart> nodes;
	std::vector<EdgeIDType> edge_refrences;
	std::vector<Edge> edges;
};


template <typename W, typename N, typename E>
class AdjacencyArray<W, N, E>::edge_id_iterator
{
	typedef typename AdjacencyArray<WeightType, NodeIDType, EdgeIDType>::edge_id_iterator self_type;

	// index of an edge within the edges-vector
	EdgeIDType edge_ref;

public:

	edge_id_iterator()
	: edge_ref(0)
	{}

	edge_id_iterator(const EdgeIDType &id)
	: edge_ref(id)
	{}

	self_type operator++(int )
	{
		self_type old = *this;
		edge_ref++;
		return old;
	}

	self_type& operator++()
	{
		++edge_ref;
		return *this;
	}

	bool operator==(const self_type &r) const
	{
		return edge_ref==r.edge_ref;
	}

	bool operator!=(const self_type &r) const
	{
		return !(*this==r);
	}

	bool operator<(const self_type &r) const
	{
		return edge_ref<r.edge_ref;
	}

	bool operator<=(const self_type &r) const
	{
		return edge_ref<=r.edge_ref;
	}

	bool operator>(const self_type &r) const
	{
		return edge_ref>r.edge_ref;
	}

	bool operator>=(const self_type &r) const
	{
		return edge_ref>=r.edge_ref;
	}

	EdgeIDType operator*() const
	{
		return edge_ref;
	}

};


template <typename W, typename N, typename E>
class AdjacencyArray<W, N, E>::incident_edges_iterator
{
	typedef typename AdjacencyArray<WeightType, NodeIDType, EdgeIDType>::incident_edges_iterator self_type;

	// index of an edge within the edges-vector
	EdgeIDType edge_reference;
	const std::vector<EdgeIDType> *edge_refrences;

public:

	incident_edges_iterator()
	: edge_reference(0), edge_refrences(NULL)
	{}

	incident_edges_iterator(const EdgeIDType &id, const std::vector<EdgeIDType> &edge_refrences)
	: edge_reference(id), edge_refrences(&edge_refrences)
	{}

	self_type operator++(int )
	{
		self_type old = *this;
		edge_reference++;
		return old;
	}

	self_type& operator++()
	{
		++edge_reference;
		return *this;
	}

	bool operator==(const self_type &r) const
	{
		return edge_reference==r.edge_reference;
	}

	bool operator!=(const self_type &r) const
	{
		return !(*this==r);
	}

	bool operator<(const self_type &r) const
	{
		return edge_reference<r.edge_reference;
	}

	bool operator<=(const self_type &r) const
	{
		return edge_reference<=r.edge_reference;
	}

	bool operator>(const self_type &r) const
	{
		return edge_reference>r.edge_reference;
	}

	bool operator>=(const self_type &r) const
	{
		return edge_reference>=r.edge_reference;
	}

	// return edge id
	EdgeIDType operator*() const
	{
		return (*edge_refrences)[edge_reference];
	}


	// return edge reference pos
	EdgeIDType current_pos() const
	{
		return edge_reference;
	}

};


template <typename W, typename N, typename E>
AdjacencyArray<W, N, E>::AdjacencyArray(const NodeIDType node_count, const std::vector<Edge> &edge_vector)
: _num_vertices(node_count), _num_edges(edge_vector.size())
{
	for(typename std::vector<EdgeWithoutID>::const_iterator it=edge_vector.begin(); it!=edge_vector.end(); it++)
	{
		edges.push_back(Edge(it->n1, it->n2, it->weight));
	}

	std::vector<EdgeIDType> node_degrees(node_count, 0);

	// compute the degrees of the nodes
	for(typename std::vector<Edge>::const_iterator edge_it=edges.begin(); edge_it!=edges.end(); edge_it++)
	{
		node_degrees[edge_it->n1]++;
		node_degrees[edge_it->n2]++;
	}

	// fill the nodes vector (initially the values of a node n are too large, deg(n))
	nodes.reserve(node_count+1);

	NodeIDType current_node_count = 0;

	for(NodeIDType i=0; i<node_count; i++)
	{
		current_node_count += node_degrees[i];
		nodes.push_back(NodeStart(current_node_count, current_node_count));
	}

	nodes.push_back(NodeStart(current_node_count, current_node_count));

	edge_refrences.assign(current_node_count, 0);

	// fill the edge_refrence vector and adjust the nodes array
	for(EdgeIDType i=0; i<_num_edges; i++)
	{
		edge_refrences[--nodes[edges[i].n1].start] = i;
		edges[i].n1_edge_ref_pos = nodes[edges[i].n1].start;

		edge_refrences[--nodes[edges[i].n2].start] = i;
		edges[i].n2_edge_ref_pos = nodes[edges[i].n2].start;
	}

}


template <typename W, typename N, typename E>
AdjacencyArray<W, N, E>::AdjacencyArray(const NodeIDType node_count, const std::vector<EdgeWithoutID> &edge_vector)
: _num_vertices(node_count), _num_edges(edge_vector.size())
{
	for(typename std::vector<EdgeWithoutID>::const_iterator it=edge_vector.begin(); it!=edge_vector.end(); it++)
	{
		edges.push_back(Edge(it->n1, it->n2, it->weight));
	}

	std::vector<EdgeIDType> node_degrees(node_count, 0);

	// compute the degrees of the nodes
	for(typename std::vector<Edge>::const_iterator edge_it=edges.begin(); edge_it!=edges.end(); edge_it++)
	{
		node_degrees[edge_it->n1]++;
		node_degrees[edge_it->n2]++;
	}

	// fill the nodes vector (initially the values of a node n are too large, deg(n))
	nodes.reserve(node_count+1);

	NodeIDType current_node_count = 0;

	for(NodeIDType i=0; i<node_count; i++)
	{
		current_node_count += node_degrees[i];
		nodes.push_back(NodeStart(current_node_count, current_node_count));
	}

	nodes.push_back(NodeStart(current_node_count, current_node_count));

	edge_refrences.assign(current_node_count, 0);

	// fill the edge_refrence vector and adjust the nodes array
	for(EdgeIDType i=0; i<_num_edges; i++)
	{
		edge_refrences[--nodes[edges[i].n1].start] = i;
		edges[i].n1_edge_ref_pos = nodes[edges[i].n1].start;

		edge_refrences[--nodes[edges[i].n2].start] = i;
		edges[i].n2_edge_ref_pos = nodes[edges[i].n2].start;
	}

}





template <typename W, typename N, typename E>
inline bool AdjacencyArray<W, N, E>::is_inactive(EdgeIDType edge_id) const
{
	return (edges[edge_id].n1_edge_ref_pos >= nodes[edges[edge_id].n1].inactive_edges_start);
}


template <typename W, typename N, typename E>
inline bool AdjacencyArray<W, N, E>::is_active(EdgeIDType edge_id) const
{
	return (edges[edge_id].n1_edge_ref_pos < nodes[edges[edge_id].n1].inactive_edges_start);
}


template <typename W, typename N, typename E>
void AdjacencyArray<W, N, E>::deactivate_edge(EdgeIDType edge_id)
{
	if(is_inactive(edge_id))
	{
		return ;
	}

	Edge current_edge = edges[edge_id];
	EdgeIDType new_edge_ref_pos;
	EdgeIDType swap_edge_ref_pos;
	EdgeIDType swap_edge_id;

	// set n1-edge inactive
	new_edge_ref_pos = --nodes[current_edge.n1].inactive_edges_start;

	swap_edge_ref_pos = current_edge.n1_edge_ref_pos;
	swap_edge_id = edge_refrences[new_edge_ref_pos];

	edge_refrences[swap_edge_ref_pos] = swap_edge_id;

	// we have to compare the edge_ref_pos, otherwise there might be some problems when swap_edge.n1==swqp_edge.n2
	if(edges[swap_edge_id].n1_edge_ref_pos==new_edge_ref_pos)
	{
		edges[swap_edge_id].n1_edge_ref_pos = swap_edge_ref_pos;
	}
	else
	{
		edges[swap_edge_id].n2_edge_ref_pos = swap_edge_ref_pos;
	}

	edge_refrences[new_edge_ref_pos] = edge_id;
	edges[edge_id].n1_edge_ref_pos = new_edge_ref_pos;


	current_edge = edges[edge_id]; // might have been changed in the last step if (n1==n2 and n2_ref was the new_edge_ref_pos)

	// set n2-edge inactive
	new_edge_ref_pos = --nodes[current_edge.n2].inactive_edges_start;

	swap_edge_ref_pos = current_edge.n2_edge_ref_pos;
	swap_edge_id = edge_refrences[new_edge_ref_pos];

	edge_refrences[swap_edge_ref_pos] = swap_edge_id;

	// we have to compare the edge_ref_pos, otherwise there might be some problems when swap_edge.n1==swqp_edge.n2
	/// \todo is that correct or could the changes one step earlier break the state of the swap-edge (n1==n2)?
	if(edges[swap_edge_id].n2_edge_ref_pos==new_edge_ref_pos)
	{
		edges[swap_edge_id].n2_edge_ref_pos = swap_edge_ref_pos;
	}
	else
	{
		edges[swap_edge_id].n1_edge_ref_pos = swap_edge_ref_pos;
	}

	edge_refrences[new_edge_ref_pos] = edge_id;
	edges[edge_id].n2_edge_ref_pos = new_edge_ref_pos;

}


template <typename W, typename N, typename E>
W AdjacencyArray<W, N, E>::get_weight() const
{
	W result = 0;


	for(EdgeIDType i=0; i<_num_edges; i++)
	{
		result += edges[i].weight;
	}

	return result;
}


} // end namespace dipl
















