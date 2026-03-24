#ifndef _UNION_FIND_H_
#define _UNION_FIND_H_


#include <vector>

#include "log_n_type.h"


namespace dipl
{


template <typename ElementType>
class UnionFind
{
	typedef typename LogNType<ElementType>::Type DepthType;

public:

	UnionFind(ElementType size)
	: size(size), parent(size), depth(size), visited(size, false)
	{
		for(unsigned int i=0; i<parent.size(); i++)
		{
			parent[i] = i;
			depth[i] = 0;
		}
	}

	/**
	 * return the representative of the set element belongs to
	 *
	 * @param
	 * @return
	 */
	ElementType find(const ElementType element)
	{
		ElementType representative = element;
		while(parent[representative]!=representative)
		{
			representative = parent[representative];
		}

		// path-compression
		ElementType current_element = element, next_element;

		while(parent[current_element]!=representative)
		{
			next_element = parent[current_element];
			parent[current_element] = representative;
			current_element = next_element;
		}

		return representative;
	}


	/*!
	 * Union of the 2 sets of element1 and element2
	 *
	 * @param element1
	 * @param element2
	 */
	void combine(const ElementType element1, const ElementType element2)
	{
		ElementType p1 = find(element1);
		ElementType p2 = find(element2);

		if(p1!=p2)
		{
			link(p1, p2);
		}
	}


	void reset(const ElementType e)
	{
		parent[e] = e;
		depth[e] = 0;
		visited[e] = false;
	}

	bool get_visited(const ElementType e) const
	{
		return visited[e];
	}

	void set_visited(const ElementType e, const bool val)
	{
		visited[e] = val;
	}

private:

	void link(const ElementType element1, const ElementType element2)
	{
		// element1 and element 2 must be representatives of 2 different sets

		if(depth[element1]>depth[element2])
		{
			parent[element2] = element1;
		}
		else
		{
			parent[element1] = element2;

			depth[element2] += (depth[element1]==depth[element2]);
		}
	}


private:
	ElementType size;

	std::vector<ElementType> parent;

	// stores the highest possible depth of each set
	// not unlikely that the actual depth is smaller
	// because of path compression
	std::vector<DepthType> depth;

	std::vector<bool> visited;
};


} // end of namespace dipl


#endif

