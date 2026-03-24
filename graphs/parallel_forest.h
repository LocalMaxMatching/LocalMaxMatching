#ifndef _PARALLEL_FOREST_H_
#define _PARALLEL_FOREST_H_

#include <list>
#include <vector>
#include <utility>

#include <mpi.h>


template <typename _NodeIDType, typename _EdgeIDType>
class AdjacencyList
{
public:
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

	typedef typename std::list< std::pair<EdgeIDType, NodeIDType> >::const_iterator incident_edge_iterator;

private:
	std::vector< std::list< std::pair<EdgeIDType, NodeIDType> > > nodes;

	/// \todo rename visited to is_root (see notes)
	std::vector<bool> visited_nodes; // initially not visited

	std::vector<NodeIDType> root_of_tree;

public:

	AdjacencyList(const NodeIDType num_vertices)
	: nodes(num_vertices, std::list< std::pair<EdgeIDType, NodeIDType> >()),
	  visited_nodes(num_vertices, false),
	  root_of_tree(num_vertices, -1)
	{}

	~AdjacencyList()
	{
	}

	void add_incident_edge(const NodeIDType src_node, const NodeIDType tgt_node, const EdgeIDType edge_id)
	{
		nodes[src_node].push_back(std::make_pair(edge_id, tgt_node));
	}

	NodeIDType num_vertices() const
	{
		return nodes.size();
	}


	incident_edge_iterator begin_incident_edges(const NodeIDType current_node) const
	{
		return nodes[current_node].begin();
	}

	incident_edge_iterator end_incident_edges(const NodeIDType current_node) const
	{
		return nodes[current_node].end();
	}

	NodeIDType out_degree(const NodeIDType node) const
	{
		return nodes[node].size();
	}

	void clear(const NodeIDType n_id)
	{
		nodes[n_id].clear();
	}

	bool visited(const NodeIDType &nid) const
	{
		return visited_nodes[nid];
	}

	void set_visited(const NodeIDType &nid, const bool val)
	{
		visited_nodes[nid] = val;
	}

	NodeIDType get_root_of_tree(const NodeIDType n) const
	{
		return root_of_tree[n];
	}

	template <class Graph>
	friend void make_tree(const Graph &g,
						  const typename Graph::NodeIDType current_node,
						  const typename Graph::NodeIDType start_parent_node,
						  AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest);

	template <class Graph>
	friend void print_local_tree(const Graph &g, typename Graph::NodeIDType root,
								 const AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest);


	template <class Graph>
		friend void print_local_tree(const Graph &g, typename Graph::NodeIDType root,
									 const AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest,
									 std::string indent);

};



/*!
 * This function traverses the given tree (not necessarily yet a tree) using a BFS or DFS to compute
 * a spanning directed tree.
 *
 * is start_node is the real root of the tree then start_parent_node must be equal to start node.
 * start_parent_node is the root of the tree
 *
 * @param g
 * @param start_node
 * @param start_parent_node
 * @param tree
 * @param visited
 */
template <class Graph>
void make_tree(const Graph &g,
			   const typename Graph::NodeIDType start_node,
			   const typename Graph::NodeIDType start_parent_node,
			   AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest)
{
	typedef typename Graph::EdgeIDType 	EdgeIDType;
	typedef typename Graph::NodeIDType 	NodeIDType;

	NodeIDType current_node, parent_node;


//	// list implementation -> BFS
//	forest.set_visited(start_node, true);

	// first ID is ID of active node and second ID is the ID of the parent node
	std::list< std::pair<NodeIDType, NodeIDType> > active_nodes_list;

	active_nodes_list.push_back(std::make_pair(start_node, start_parent_node));

	while(!active_nodes_list.empty())
	{
		current_node = active_nodes_list.front().first;
		parent_node = active_nodes_list.front().second;
		forest.root_of_tree[current_node] = start_parent_node;
		active_nodes_list.pop_front();

		for(typename std::list< std::pair<EdgeIDType, NodeIDType> >::iterator e_it=forest.nodes[current_node].begin(); e_it!=forest.nodes[current_node].end(); )
		{
			if(e_it->second != parent_node)
			{// we don't want to go back!!
				if(g.is_local(current_node))
				{
					// we don't follow edges starting at ghostnodes
//					forest.set_visited(e_it->second, true); // well not yet visited but it will be visited
					active_nodes_list.push_back(std::make_pair(e_it->second, current_node));
				}
				e_it++; // advance to next edge
			}
			else
			{// delete edge to parent node
				e_it = forest.nodes[current_node].erase(e_it); // delete edge, edge has been visited by another path
			}
		}
	}

}


template <class Graph>
void print_local_tree(const Graph &g,
						typename Graph::NodeIDType root,
							 const AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest)
{
	typedef typename AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType>::incident_edge_iterator incident_edge_iterator;

	std::cout << g.local_vertex_id_to_global_id(root) << std::endl;

	for(incident_edge_iterator e_it = forest.begin_incident_edges(root);
			e_it!=forest.end_incident_edges(root); e_it++)
	{
		print_local_tree(g, e_it->second, forest, std::string("\t"));
	}
}


template <class Graph>
void print_local_tree(const Graph &g,
		typename Graph::NodeIDType root,
							 const AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &forest,
							 std::string indent)
{
	typedef typename AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType>::incident_edge_iterator incident_edge_iterator;

	std::cout << indent << g.local_vertex_id_to_global_id(root) << std::endl;

	if(g.is_ghost(root))
	{
		return;
	}

	for(incident_edge_iterator e_it = forest.begin_incident_edges(root);
			e_it!=forest.end_incident_edges(root); e_it++)
	{
		print_local_tree(g, e_it->second, forest, indent+"\t");
	}
}


#endif


