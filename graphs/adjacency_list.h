#ifndef ADJACENCY_LIST_H_
#define ADJACENCY_LIST_H_


#include <list>
#include <vector>
#include <utility>



namespace dipl
{

template <typename _NodeIDType, typename _WeightType>
class AdjacencyList
{
public:
	typedef _NodeIDType NodeIDType;
	typedef _WeightType WeightType;

	typedef typename std::vector< std::pair<NodeIDType, WeightType> >::const_iterator incident_edge_iterator;

private:
	std::vector< std::pair<WeightType, std::vector< std::pair<NodeIDType, WeightType> > > > nodes;

public:

	AdjacencyList()
	{}

	AdjacencyList(const NodeIDType num_vertices)
	: nodes(num_vertices, std::pair<WeightType, std::vector< std::pair<NodeIDType, WeightType> > >())
	{}

	~AdjacencyList()
	{
	}

	void resize(const NodeIDType num_vertices)
	{
		nodes.resize(num_vertices, std::pair<WeightType, std::vector< std::pair<NodeIDType, WeightType> > >());
	}

	void set_vertex_weight(const NodeIDType src_node,
						   const WeightType weight)
	{
		nodes[src_node].first = weight;
	}

	void add_incident_edge(const NodeIDType src_node,
						   const NodeIDType tgt_node, const WeightType edge_weight)
	{
		nodes[src_node].second.push_back(std::make_pair(tgt_node, edge_weight));
	}

	NodeIDType num_vertices() const
	{
		return nodes.size();
	}


	incident_edge_iterator begin_incident_edges(const NodeIDType current_node) const
	{
		return nodes[current_node].second.begin();
	}

	incident_edge_iterator end_incident_edges(const NodeIDType current_node) const
	{
		return nodes[current_node].second.end();
	}

	WeightType get_node_weight(const NodeIDType n) const
	{
		return nodes[n].first;
	}

	NodeIDType out_degree(const NodeIDType node) const
	{
		return nodes[node].second.size();
	}

	void clear(const NodeIDType n_id)
	{
		nodes[n_id].second.clear();
	}

};



}




#endif





