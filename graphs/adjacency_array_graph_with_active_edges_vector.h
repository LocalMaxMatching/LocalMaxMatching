
#include <stdlib.h>
#include <vector>

#include "adjacency_array.h"

namespace dipl
{


template <typename _WeightType=double,
		  typename _NodeIDType=unsigned int,
		  typename _EdgeIDType=unsigned int>
class AdjacencyArrayGraphWithActiveEdgesVector: public AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>
{
public:

	typedef _WeightType WeightType;
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

	typedef typename AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>::Edge Edge;
	typedef typename AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>::incident_edges_iterator incident_edges_iterator;

	typedef typename AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>::EdgeWithoutID EdgeWithoutID;

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
	AdjacencyArrayGraphWithActiveEdgesVector(const NodeIDType node_count, const std::vector<Edge> &edge_vector)
	: AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>(node_count, edge_vector)
	{
		edge_management.reserve(this->num_edges());
		active_edge_pos.reserve(this->num_edges());

		for(EdgeIDType i=0; i<this->num_edges(); i++)
		{
			edge_management.push_back(i);
			active_edge_pos.push_back(i);
		}

		_end_active_edges = this->num_edges();
	}


	/*!
	 * \brief Constructor that initializes the graph from the given edge list.
	 *
	 * Runtime: O(m+n)
	 *
	 * @param node_count	Number of vertices of the graph, defined by the edge list.
	 * @param edge_list		List of edges that define the graph. The vertices of the
	 * 						graph have to be numbered from 0 to node_count-1
	 */
	AdjacencyArrayGraphWithActiveEdgesVector(const NodeIDType node_count, const std::vector<EdgeWithoutID> &edge_vector)
	: AdjacencyArray<_WeightType, _NodeIDType, _EdgeIDType>(node_count, edge_vector)
	{
		edge_management.reserve(this->num_edges());
		active_edge_pos.reserve(this->num_edges());

		for(EdgeIDType i=0; i<this->num_edges(); i++)
		{
			edge_management.push_back(i);
			active_edge_pos.push_back(i);
		}

		_end_active_edges = this->num_edges();
	}


	virtual ~AdjacencyArrayGraphWithActiveEdgesVector() {}


//	/*!
//	 * \brief Returns an edge_id_iterator that points to the first edge of the graph.
//	 *
//	 * Runtime: O(1)
//	 *
//	 * @return	edge_id_iterator of the first edge of the graph.
//	 */
//	edge_iterator begin_edges() const
//	{
//		return edge_id_iterator(0, this);
//	}
//
//
//	/*!
//	 * \brief Returns an edge_id_iterator that indicates the end.
//	 *
//	 * Runtime: O(1)
//	 *
//	 * @return	edge_id_iterator of the end.
//	 */
//	edge_iterator end_edges() const
//	{
//		return edge_id_iterator(_num_edges, this);
//	}
//
//
//
//	edge_iterator begin_active_edges() const
//	{
//		return edge_id_iterator(0, this);
//	}
//
//
//	edge_iterator end_active_edges() const
//	{
//		return edge_id_iterator(_end_active_edges, this);
//	}


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
	inline bool is_active(EdgeIDType edge_id) const
	{
		return active_edge_pos[edge_id]<_end_active_edges;
	}


	/*!
	 * \brief Deactives the given edge ("adds it to the list of inactive edges");
	 *
	 * Runtime: O(1)
	 *
	 * @param edge_id	The ID of the edge to be deactivated.
	 */
	virtual void deactivate_edge(EdgeIDType edge_id)
	{
		AdjacencyArray<WeightType, NodeIDType, EdgeIDType>::deactivate_edge(edge_id);

		_end_active_edges--;

		EdgeIDType swap_edge_id = edge_management[_end_active_edges];

		edge_management[active_edge_pos[edge_id]] = swap_edge_id;
		active_edge_pos[swap_edge_id] = active_edge_pos[edge_id];

		edge_management[_end_active_edges] = edge_id;
		active_edge_pos[edge_id] = _end_active_edges;
	}


	EdgeIDType random_active_edge() const
	{
		return edge_management[rand()%_end_active_edges];
	}

	bool has_active_edges() const
	{
		return _end_active_edges;
	}

	EdgeIDType num_active_edges() const
	{
		return _end_active_edges;
	}


//	bool correct_state() const
//	{
//		for(EdgeIDType e=0; e<_num_edges; e++)
//		{
//			if( ((edges[e].n1_edge_ref_pos >= nodes[edges[e].n1].inactive_edges_start) && (edges[e].n2_edge_ref_pos < nodes[edges[e].n2].inactive_edges_start)) ||
//				((edges[e].n1_edge_ref_pos < nodes[edges[e].n1].inactive_edges_start) && (edges[e].n2_edge_ref_pos >= nodes[edges[e].n2].inactive_edges_start)) )
//			{
//				std::cout << "incorrect state of edge " << e << ":" << std::endl;
//				std::cout << edges[e].n1_edge_ref_pos << ">=" << nodes[edges[e].n1].inactive_edges_start << " && " << edges[e].n2_edge_ref_pos << "<" << nodes[edges[e].n2].inactive_edges_start << std::endl;
//				std::cout << edges[e].n1_edge_ref_pos << "<" << nodes[edges[e].n1].inactive_edges_start << " && " << edges[e].n2_edge_ref_pos << ">=" << nodes[edges[e].n2].inactive_edges_start << std::endl;
//
//				return false;
//			}
//		}
//
//		return true;
//	}

private:

	std::vector<EdgeIDType> edge_management;
	std::vector<EdgeIDType> active_edge_pos;


	EdgeIDType _end_active_edges;
};




}
