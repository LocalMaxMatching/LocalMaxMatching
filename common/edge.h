#ifndef _DIPL_EDGE_H
#define _DIPL_EDGE_H

#include <iostream>

namespace dipl
{

template <typename WeightType, typename NodeIDType>
struct Edge
{
	NodeIDType n1;
	NodeIDType n2;

	WeightType weight;

	Edge()
	: n1(0), n2(0), weight(0)
	{}

	Edge(const NodeIDType n1, const NodeIDType n2, const WeightType weight)
	: n1(n1), n2(n2), weight(weight)
	{}

	friend std::ostream& operator<<(std::ostream &os, const Edge &e)
	{
		os << "(" << e.n1 << " <-> " << e.n2 << ", " << e.weight << ")";
		return os;
	}

};


template <typename WeightType, typename NodeIDType, typename EdgeIDType>
struct EdgeWithID
{
	NodeIDType n1;
	NodeIDType n2;

	EdgeIDType edge_id;

	WeightType weight;

	EdgeWithID()
	: n1(0), n2(0), edge_id(0), weight(0)
	{}

	EdgeWithID(const NodeIDType n1, const NodeIDType n2, const EdgeIDType edge_id, const WeightType weight)
	: n1(n1), n2(n2), edge_id(edge_id), weight(weight)
	{}

	friend std::ostream& operator<<(std::ostream &os, const EdgeWithID &e)
	{
		os << "(" << e.n1 << " <-> " << e.n2 << ", " << e.edge_id << ", " << e.weight << ")";
		return os;
	}

};

}


#endif

