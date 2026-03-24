
#include <utility>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <graphs/tree.h>

#include <common/math_funcs.h>

/*!
 * create the (dynamic) table of costs for each edge (that's the incoming edge of the tgt_node)
 * if this edge is a matched or not
 *
 * @param g
 * @param tree
 * @param tgt_node
 * @param weight_to_node
 * @param local_results
 * @return
 */
template <class TreeGraph, class Graph>
//std::pair<typename Graph::WeightType, typename Graph::WeightType>
inline void
tree_max_weighted_matching(const Graph &g,
						   const TreeGraph &tree,
						   const typename TreeGraph::NodeIDType &tgt_node,
						   const typename Graph::WeightType &weight_to_node,
						   std::vector<std::pair<std::pair<typename Graph::WeightType, typename Graph::WeightType>, typename Graph::NodeIDType> > &local_results)
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::WeightType WeightType;

	std::list<NodeIDType> active_nodes;
	std::vector<std::pair<NodeIDType, WeightType> > bfs_traversal;
	active_nodes.push_back(tgt_node);
	bfs_traversal.push_back(std::make_pair(tgt_node, weight_to_node));

	// create BFS-traversal
	while(!active_nodes.empty())
	{
		NodeIDType current_node   = active_nodes.front();
		active_nodes.pop_front();

		const typename TreeGraph::incident_edge_iterator end = tree.end_incident_edges(current_node);
		for(typename TreeGraph::incident_edge_iterator
				e_it=tree.begin_incident_edges(current_node);
				e_it!=end; e_it++)
		{
			active_nodes.push_back(e_it->second);
			bfs_traversal.push_back(std::make_pair(e_it->second, g.get_edge(e_it->first).weight));
		}
	}

	// taverse the tree bottom up
	for(NodeIDType i=bfs_traversal.size(); i>0; )
	{
		i--; // decrementing it here -> because NodeIDType might be unsigned

		NodeIDType current_node   = bfs_traversal[i].first;
		WeightType current_weight = bfs_traversal[i].second;

		if(tree.out_degree(current_node)==0)
		{
			local_results[current_node].first = std::make_pair(current_weight, 0);
			local_results[current_node].second = -1; // there are no outgoing edges to be used!!
		}
		else
		{
			// case: edge to this subtree is matched
			WeightType matching_weight1 = current_weight;

			const typename TreeGraph::incident_edge_iterator end = tree.end_incident_edges(current_node);
			for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(current_node);
					e_it!=end; e_it++)
			{
				// already computed the local result from each following node -> traversing the tree bottom up
				matching_weight1 += local_results[e_it->second].first.second;
			}

			// case: edge to this subtree is not matched
			WeightType matching_weight2 = 0;

			WeightType old_diff = 0;
			local_results[current_node].second = -1; // set the outgoing partner initially to -1,
													 // in case we're not using any outgoing partner.

			for(typename TreeGraph::incident_edge_iterator
					e_it=tree.begin_incident_edges(current_node);
					e_it!=end; e_it++)
			{
				const std::pair<WeightType, WeightType> current_pair = local_results[e_it->second].first;
				const WeightType current_diff = current_pair.first-current_pair.second;

				if(current_diff>old_diff)
				{
					matching_weight2 += current_pair.first - old_diff;
					local_results[current_node].second = e_it->second; // set the target-node, who's edge has been used by the matching
					old_diff = current_diff;
				}
				else
				{
					matching_weight2 += current_pair.second;
				}
			}

			local_results[current_node].first = std::make_pair(matching_weight1, matching_weight2);
		}
	}

}


/*!
 * use the cost table to compute the max weighted matching of the tree.
 * recursively traverse the tree from the root.
 *
 * @param current_node
 * @param used_incoming_edge
 * @param tree
 * @param local_results
 * @param matching
 */
template <class TreeGraph, typename WeightType, typename NodeIDType, typename EdgeIDType>
inline void add_matching_edges(const NodeIDType &current_node,
							   const bool &used_incoming_edge,
							   const TreeGraph &tree,
							   std::vector<std::pair<std::pair<WeightType, WeightType>, NodeIDType> > &local_results,
							   std::vector<EdgeIDType> &matching)
{
	std::list<std::pair<NodeIDType, bool> > active_nodes;
	active_nodes.push_back(std::make_pair(current_node, used_incoming_edge));

	while(!active_nodes.empty())
	{
		NodeIDType current_node = active_nodes.front().first;
		bool used_incoming_edge = active_nodes.front().second;
		active_nodes.pop_front();

		if(used_incoming_edge)
		{
			for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(current_node);
					e_it!=tree.end_incident_edges(current_node); e_it++)
			{
				active_nodes.push_back(std::make_pair(e_it->second, false));
			}
		}
		else
		{
			for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(current_node); e_it!=tree.end_incident_edges(current_node); e_it++)
			{
				if(local_results[current_node].second == e_it->second)
				{
					matching.push_back(e_it->first);
					active_nodes.push_back(std::make_pair(e_it->second, true));
				}
				else
				{
					active_nodes.push_back(std::make_pair(e_it->second, false));
				}
			}
		}

	}



//	if(used_incoming_edge)
//	{
//		for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(current_node); e_it!=tree.end_incident_edges(current_node); e_it++)
//		{
//			add_matching_edges(e_it->second, false, tree, local_results, matching);
//		}
//	}
//	else
//	{
//		for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(current_node); e_it!=tree.end_incident_edges(current_node); e_it++)
//		{
// //			unsigned int e_id = e_it->first;
// //			unsigned int next_node = e_it->second;
// //			std::pair<WeightType, WeightType> tmp = local_results[current_node].first;
// //			unsigned int tmp_node = local_results[current_node].second;
//
//			if(local_results[current_node].second == e_it->second)
//			{
//				matching.push_back(e_it->first);
//				add_matching_edges(e_it->second, true, tree, local_results, matching);
//			}
//			else
//			{
//				add_matching_edges(e_it->second, false, tree, local_results, matching);
//			}
//		}
//	}
}



/*
 * In the initial call of this function (node==root), tree_size must be set to 0
 */
template <class TreeGraph>
void tree_size_and_depth(const TreeGraph &tree,
						 const typename TreeGraph::NodeIDType &node,
						 typename TreeGraph::NodeIDType &tree_size,
						 typename TreeGraph::NodeIDType &max_depth)
{
	typename TreeGraph::NodeIDType subtree_depth;

	tree_size++;
	max_depth=0;

	for(typename TreeGraph::incident_edge_iterator e_it=tree.begin_incident_edges(node); e_it!=tree.end_incident_edges(node); e_it++)
	{
		tree_size_and_depth(tree, e_it->second, tree_size, subtree_depth);

		if(subtree_depth+1>max_depth)
		{
			max_depth = subtree_depth+1; // we have to add 1 because of the edge going to the subtree
		}
	}

}

/// \todo this function needs another name
template <class Graph>
void max_weighted_matching_of_edge_list(const Graph &g,
										const std::list<typename Graph::EdgeIDType> &edges,
										std::vector<typename Graph::EdgeIDType> &matching,
										AdjacencyList<typename Graph::NodeIDType, typename Graph::EdgeIDType> &tree,
										/* first weight of a node specifies the accumulated weight of the subtree if the incoming edge is matched,
										 * the second weight weight specified the accumulated weight of the subtree if the incoming edge is not matched -> the node-entry specifies
										 * which outgoing edge is matched (if any) */
										std::vector<std::pair<std::pair<typename Graph::WeightType, typename Graph::WeightType>, typename Graph::NodeIDType> > &local_results
#ifdef MORE_INFORMATION
		,
		std::list<typename Graph::NodeIDType> &depths_of_trees,
		std::list<typename Graph::NodeIDType> &sizes_of_trees
#endif
		)
{
	std::list<typename Graph::NodeIDType> roots;

	build_forest(g, edges, tree, roots);

#ifdef MORE_INFORMATION
	// compute lists of the sizes and depths of the trees

	for(typename std::list<typename Graph::NodeIDType>::iterator n_it=roots.begin(); n_it!=roots.end(); n_it++)
	{
		typename Graph::NodeIDType size=0, depth=0;

		tree_size_and_depth(tree, *n_it, size, depth);

		sizes_of_trees.push_back(size);
		depths_of_trees.push_back(depth);
	}
#endif

	// compute the tree matchings
	for(typename std::list<typename Graph::NodeIDType>::iterator n_it=roots.begin(); n_it!=roots.end(); n_it++)
	{
		typename Graph::NodeIDType current_node = *n_it;

		tree_max_weighted_matching(g, tree, current_node, 0, local_results);

		add_matching_edges(current_node, local_results[current_node].first.first>=local_results[current_node].first.second, tree, local_results, matching);

	}

	for(typename std::list<typename Graph::EdgeIDType>::const_iterator it=edges.begin(); it!=edges.end(); it++)
	{
		const typename Graph::Edge &e = g.get_edge(*it);
		tree.clear(e.n1);
		tree.clear(e.n2);



//		local_results(g.num_vertices(), std::pair<std::pair<WeightType, WeightType>, NodeIDType>() );
//		typedef typename Graph::WeightType WeightType;
//		typedef typename Graph::NodeIDType NodeIDType;
//		local_results[e.n1] = std::pair<std::pair<WeightType, WeightType>, NodeIDType>();
//		local_results[e.n1] = std::pair<std::pair<WeightType, WeightType>, NodeIDType>();
	}

}




