#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>


#include <matching_algorithms/tree_matching.h>
#include <common/hash_functions.h>


namespace dipl
{

template <class Graph>
class LocalTreeMaximumMatching
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;

	typedef typename Graph::Edge	   Edge;

public:

	static bool better_than_old_edge(const Graph &g,
									 const WeightType old_weight, const EdgeIDType old_edge_id,
									 const WeightType new_weight, const EdgeIDType new_edge_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash2(old_edge_id);
		unsigned int new_hash = hash2(new_edge_id);

		return /*g.is_inactive(old_edge_id) ||*/
			   (new_weight>old_weight) ||
			   (new_weight==old_weight && new_hash<old_hash) ||
			   (new_weight==old_weight && new_hash==old_hash && new_edge_id<old_edge_id);

//		return g.is_inactive(old_edge_id) ||
//			   (new_weight>old_weight) ||
//			   (new_weight==old_weight && hash(new_edge_id)<hash(old_edge_id)) ||
//			   (new_weight==old_weight && hash(new_edge_id)==hash(old_edge_id) && new_edge_id<old_edge_id);
	}

	/*!
	 * \brief Computes a maximal weighted matching of the given graph.
	 *
	 * Runtime: ??
	 *
	 * @param g			The Graph
	 * @param matching	Return list of the resulting matching
	 */
	static void compute_weighted_matching(Graph &g,
									     std::list<Edge> &matching,
									     unsigned int &depth
#ifdef MORE_INFORMATION
									     ,
									     std::list<typename Graph::EdgeIDType> &active_edge_count,
									     std::list< std::list<typename Graph::NodeIDType> > &depths_of_trees,
									     std::list< std::list<typename Graph::NodeIDType> > &sizes_of_trees
#endif
									     )
	{
		typedef typename Graph::edge_iterator edge_iterator;


		std::vector<bool> is_active(g.num_vertices(), true);

		AdjacencyList<NodeIDType, EdgeIDType> tree(g.num_vertices());
		std::vector<std::pair<std::pair<WeightType, WeightType>, NodeIDType> > local_results(g.num_vertices(), std::pair<std::pair<WeightType, WeightType>, NodeIDType>() );

		// stores for each vertex the position of the local dominant edge
		std::vector<EdgeIDType> max_edge_of_node(g.num_vertices(), g.num_edges());

		depth = 0;

		while(g.has_active_edges())
		{
#ifdef MORE_INFORMATION
			active_edge_count.push_back(g.num_active_edges());
#endif

			depth++;

			// iterate over each active edge, to set local max edges
			for(edge_iterator e_it=g.begin_active_edges(); e_it<g.end_active_edges(); e_it++)
			{
				// check first endpoint of current edge
				if(better_than_old_edge(g,
										g.get_edge(max_edge_of_node[e_it->n1]).weight, g.get_edge(max_edge_of_node[e_it->n1]).edge_id,
										e_it->weight,                      e_it->edge_id))
				{
					max_edge_of_node[e_it->n1] = g.get_edge_pos(e_it);
				}

				// check second endpoint of current edge
				if(better_than_old_edge(g,
										g.get_edge(max_edge_of_node[e_it->n2]).weight, g.get_edge(max_edge_of_node[e_it->n2]).edge_id,
										e_it->weight,                      e_it->edge_id))
				{
					max_edge_of_node[e_it->n2] = g.get_edge_pos(e_it);
				}

			}

			std::list<EdgeIDType> tree_edges; // use vector instead of list - maybe use a counter somewhere for the correct size, if possible
			std::vector<EdgeIDType> tree_matching;

			// get tree edges, those are the incident edges that are maximal
			for(edge_iterator e_it=g.begin_active_edges(); e_it<g.end_active_edges(); e_it++)
			{
				if(g.get_edge_id(e_it)==g.get_edge(max_edge_of_node[e_it->n1]).edge_id)
				{
//					std::cout << "adding tree edge: " << g.get_edge(max_edge_of_node[e_it->n1]) << std::endl;
					tree_edges.push_back(max_edge_of_node[e_it->n1]);
				}
				else if(g.get_edge_id(e_it)==g.get_edge(max_edge_of_node[e_it->n2]).edge_id)
				{
//					std::cout << "adding tree edge: " << g.get_edge(max_edge_of_node[e_it->n2]) << std::endl;
					tree_edges.push_back(max_edge_of_node[e_it->n2]);
				}
			}


			// compute tree-matching
#ifdef MORE_INFORMATION
			std::list<typename Graph::NodeIDType> round_depths_of_trees;
			std::list<typename Graph::NodeIDType> round_sizes_of_trees;

			max_weighted_matching_of_edge_list(g, tree_edges, tree_matching, tree, local_results, round_depths_of_trees, round_sizes_of_trees);

			depths_of_trees.push_back(round_depths_of_trees);
			sizes_of_trees.push_back(round_sizes_of_trees);
#else
			max_weighted_matching_of_edge_list(g, tree_edges, tree_matching, tree, local_results);
#endif

			// add matched edges and set the incident nodes to inactive (matched)
			for(typename std::vector<EdgeIDType>::iterator e_it=tree_matching.begin(); e_it!=tree_matching.end(); e_it++)
			{
				const Edge &e = g.get_edge(*e_it);
				matching.push_back(e);
				is_active[e.n1] = false;
				is_active[e.n2] = false;
			}


			// deactivate edges that are incident to matched nodes
			edge_iterator e_it=g.begin_active_edges();
			while(e_it<g.end_active_edges())
			{
				const Edge &e = *e_it;

				// check first endpoint of current edge
				if(!is_active[e.n1] || !is_active[e.n2])
				{
					// deactivate edge
					g.deactivate_edge(e_it);
					continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
				}

				// reset the max edges of nodes
				// make sure that node still incident to active edges
				// don't store old information
				max_edge_of_node[e.n1] = g.num_edges();
				max_edge_of_node[e.n2] = g.num_edges();

				e_it++;
			}

		}

	}


	/*!
	 * \brief Checks if the given matching is a correct and maximal matching of the given graph.
	 *
	 * Runtime: O(m)
	 * Also true for incorrect matchings, because we abort as soon as we find an incorrect matching-edge.
	 *
	 * @param matchings
	 * @param g
	 * @return
	 */
	static bool is_maximal_matching(std::list<Edge> &matchings,
							 	    Graph &g,
							 	    NodeIDType &unmatched_nodes)
	{
		bool result = true;

		std::vector<bool> matched_edge(g.num_edges(), false);
		std::vector<Edge> incident_edge(g.num_vertices(), Edge(0,0,0, std::numeric_limits<EdgeIDType>::max()) );
		std::vector<bool> matched_node(g.num_vertices(), false);

		for(typename std::list<Edge>::const_iterator it=matchings.begin(); it!=matchings.end(); it++)
		{
			matched_node[it->n1] = true;
			matched_node[it->n2] = true;

			incident_edge[it->n1] = *it;
			incident_edge[it->n2] = *it;

			matched_edge[it->edge_id] = true;
		}

		// check if the given matching is a correct matching
		for(typename Graph::edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
		{
			if(matched_edge[e_it->edge_id] && (incident_edge[e_it->n1].edge_id != e_it->edge_id || incident_edge[e_it->n2].edge_id != e_it->edge_id) )
			{
				result = false;

				std::cerr << "Not a correct matching, " << *e_it << " is adjacent to another matched edge. " << incident_edge[e_it->n1] << " or " << incident_edge[e_it->n2] << std::endl;
			}
		}

		// check if the given matching is maximal
		for(typename Graph::edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
		{
			if(!matched_edge[e_it->edge_id])
			{
				if(!matched_node[e_it->n1] && !matched_node[e_it->n2])
				{
					result = false;

					std::cerr << "Matching isn't maximal, end points of " << g.get_edge(e_it) << " haven't been matched!" << std::endl;
				}
			}
		}

		unmatched_nodes = g.num_vertices();
		for(NodeIDType i=0; i<g.num_vertices(); i++)
		{
			unmatched_nodes -= matched_node[i];
		}

		return result;
	}


	static double get_weight(const Graph &g, const std::list<Edge> &matching)
	{
		double result = 0.;

		for(typename std::list<Edge>::const_iterator it=matching.begin(); it!=matching.end(); it++)
		{
			result += it->weight;
		}

		return result;
	}

};

}

