
#include <list>
#include <vector>
#include <utility>


template <typename _NodeIDType, typename _EdgeIDType>
class AdjacencyList
{
public:
	typedef _NodeIDType NodeIDType;
	typedef _EdgeIDType EdgeIDType;

	typedef typename std::list< std::pair<EdgeIDType, NodeIDType> >::const_iterator incident_edge_iterator;

private:
	std::vector< std::list< std::pair<EdgeIDType, NodeIDType> > > nodes;
	std::vector<bool> visited_nodes; // initially not visited

public:

	AdjacencyList(const NodeIDType num_vertices)
	: nodes(num_vertices, std::list< std::pair<EdgeIDType, NodeIDType> >()),
	  visited_nodes(num_vertices, false)
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

	template <class Graph>
	friend void build_forest(const Graph &g,
							 const std::list<typename Graph::EdgeIDType> &edges,
							 AdjacencyList<typename Graph::NodeIDType,
							 typename Graph::EdgeIDType> &tree,
							 std::list<typename Graph::NodeIDType> &roots);


	template <class Graph>
	friend void make_tree(const Graph &g,
						  const typename Graph::NodeIDType current_node,
						  AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &tree);

};



/*!
 * This function traverses the given tree (not necessarily yet a tree) using a BFS or DFS to compute
 * a spanning directed tree.
 *
 * @param g
 * @param start_node
 * @param tree
 * @param visited
 */
template <class Graph>
void make_tree(const Graph &g, const typename Graph::NodeIDType start_node, AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &tree)
{
	typedef typename Graph::EdgeIDType 	EdgeIDType;
	typedef typename Graph::NodeIDType 	NodeIDType;

	NodeIDType current_node;


	// list implementation -> BFS
	tree.set_visited(start_node, true);

	std::list<NodeIDType> active_nodes_list;

	active_nodes_list.push_back(start_node);

	while(!active_nodes_list.empty())
	{
		current_node = active_nodes_list.front();
		active_nodes_list.pop_front();

		for(typename std::list< std::pair<EdgeIDType, NodeIDType> >::iterator e_it=tree.nodes[current_node].begin(); e_it!=tree.nodes[current_node].end(); )
		{
			if(!tree.visited(e_it->second))
			{
				tree.set_visited(e_it->second, true); // well not yet visited but it will be visited
				active_nodes_list.push_back(e_it->second);
				e_it++; // advance to next edge
			}
			else
			{
				e_it = tree.nodes[current_node].erase(e_it); // delete edge, edge has been visited by another path
			}
		}
	}



//	// stack implementation -> DFS
//	std::vector<NodeIDType> active_nodes_stack;
//	visited[start_node] = true;
//	active_nodes_stack.push_back(start_node);
//
//	while(!active_nodes_stack.empty())
//	{
//		current_node = active_nodes_stack.back();
//		active_nodes_stack.pop_back();
//
//		for(typename std::list< std::pair<EdgeIDType, NodeIDType> >::iterator e_it=tree.nodes[current_node].begin(); e_it!=tree.nodes[current_node].end(); )
//		{
//			if(!visited[e_it->second])
//			{
//				visited[e_it->second] = true; // well not yet visited but it will be visited
//				active_nodes_stack.push_back(e_it->second);
//				e_it++; // advance to next edge
//			}
//			else
//			{
//				e_it = tree.nodes[current_node].erase(e_it); // delete edge, edge has been visited by another path
//			}
//		}
//	}

//	// recursive stack implementation -> DFS
//	current_node = start_node;
//	visited[current_node] = true;
//
//	for(typename std::list< std::pair<EdgeIDType, NodeIDType> >::iterator e_it=tree.nodes[current_node].begin(); e_it!=tree.nodes[current_node].end(); )
//	{
//		if(!visited[e_it->second])
//		{
//			make_tree(g,e_it->second, tree, visited);
//			e_it++; // advance to next edge
//		}
//		else
//		{
//			e_it = tree.nodes[current_node].erase(e_it); // delete edge, edge has been visited by another path
//		}
//	}

}


// given tree must be set to the correct number of possible nodes (i.e. the hight possible node-id + 1)
template <class Graph>
void build_forest(const Graph &g,
				  const std::list<typename Graph::EdgeIDType> &edges,
				  AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &tree,
				  std::list<typename Graph::NodeIDType> &roots)
{
	typedef typename Graph::Edge 		Edge;
	typedef typename Graph::EdgeIDType 	EdgeIDType;
	typedef typename Graph::NodeIDType 	NodeIDType;

//	std::vector<bool> visited(tree.num_vertices(), false);

/// \todo build the forest using a union-find like structure?? - could be used to identify connected components (trees).
/// but I assume not directly to create the trees. because when combining two subtrees the nodes that are used to link
/// those subtrees are most likely not the representatives (i.e. the "directions" of the trees is wrong)

	for(typename std::list<typename Graph::EdgeIDType>::const_iterator e_it=edges.begin(); e_it!=edges.end(); e_it++)
	{
		const Edge &e = g.get_edge(*e_it);

		tree.add_incident_edge(e.n1, e.n2, *e_it);
		tree.add_incident_edge(e.n2, e.n1, *e_it);

		tree.set_visited(e.n1, false);
		tree.set_visited(e.n2, false);
	}

	for(typename std::list<typename Graph::EdgeIDType>::const_iterator e_it=edges.begin(); e_it!=edges.end(); e_it++)
	{
		const Edge &e = g.get_edge(*e_it);

		if(!tree.visited(e.n1))
		{
			tree.set_visited(e.n1, true);
			roots.push_back(e.n1);

			for(typename std::list< std::pair<EdgeIDType, NodeIDType> >::iterator subtree_e_it=tree.nodes[e.n1].begin(); subtree_e_it!=tree.nodes[e.n1].end(); )
			{
				if(!tree.visited(subtree_e_it->second))
				{
					make_tree(g,subtree_e_it->second, tree);
					subtree_e_it++; // advance to next edge
				}
				else
				{
					subtree_e_it = tree.nodes[e.n1].erase(subtree_e_it); // delete edge, edge has been visited by another path
				}
			}
		}
	}

}













