#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <boost/functional/hash.hpp>

#include <common/hash_functions.h>
#include <common/math_funcs.h>


namespace dipl
{

template <class Graph>
class LocalMaximumMatching
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;

	typedef typename Graph::Edge	   Edge;

public:




	static bool better_than_old_edge(const WeightType old_weight, const EdgeIDType old_edge_id,
									 const WeightType new_weight, const EdgeIDType new_edge_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = simple_random_hash(old_edge_id);//hash2(old_edge_id);
		unsigned int new_hash = simple_random_hash(new_edge_id);//hash2(new_edge_id);

		return (new_weight>old_weight) ||
			   (new_weight==old_weight && new_hash<old_hash) ||
			   (new_weight==old_weight && new_hash==old_hash && new_edge_id<old_edge_id);
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
									     std::list<typename Graph::Edge> &matching,
									     unsigned int &depth
#ifdef MORE_INFORMATION
									     ,
									     std::list<typename Graph::EdgeIDType> &active_edge_count
#endif
									    )
	{
		typedef typename Graph::edge_iterator edge_iterator;

		std::vector<bool> is_active(g.num_vertices(), true);

		// stores for each vertex the incident edge id, of the edge with the largest weight
		Edge dummy_edge = Edge(0,0, g.num_edges(), min_val<WeightType>());
		std::vector<Edge> max_edge_of_node(g.num_vertices(), Edge(0,0, g.num_edges(), min_val<WeightType>())); // initialize with dummy edge id, this edge is always inactive

		depth = 0;

		while(g.has_active_edges())
		{
			depth++;

#ifdef MORE_INFORMATION
			// set random edge weights just for testing if the initial random weights
			// correspond to random weight through the rounds

			for(edge_iterator e_it=g.begin_active_edges(); e_it<g.end_active_edges(); e_it++)
			{
				g.get_edge(e_it).weight = random_weight();
			}

			active_edge_count.push_back(g.num_active_edges());
#endif

			// iterate over each active edge, to set local max edges
			for(edge_iterator e_it=g.begin_active_edges(); e_it<g.end_active_edges(); e_it++)
			{
				// check first endpoint of current edge
				if(better_than_old_edge(max_edge_of_node[e_it->n1].weight, max_edge_of_node[e_it->n1].edge_id,
										e_it->weight,                      e_it->edge_id) )
				{
					max_edge_of_node[e_it->n1] = *e_it;
				}

				// check first endpoint of current edge
				if(better_than_old_edge(max_edge_of_node[e_it->n2].weight, max_edge_of_node[e_it->n2].edge_id,
										e_it->weight,                      e_it->edge_id))
				{
					max_edge_of_node[e_it->n2] = *e_it;
				}

			}


			// add matching-edges to matchings-list and set them inactive, also set the incident nodes inactive
			edge_iterator e_it=g.begin_active_edges();
			while(e_it<g.end_active_edges())
			{
				// check first endpoint of current edge
				if(g.get_edge_id(e_it)==max_edge_of_node[e_it->n1].edge_id && g.get_edge_id(e_it)==max_edge_of_node[e_it->n2].edge_id)
				{
//					std::cout << (*e_it) << std::endl;
					matching.push_back(*e_it);
					is_active[e_it->n1] = false;
					is_active[e_it->n2] = false;

					g.deactivate_edge(e_it);
					continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
				}

				e_it++;
			}


			// deactivate edges that are incident to matched nodes
			e_it=g.begin_active_edges();
			while(e_it<g.end_active_edges())
			{
				// check first endpoint of current edge
				if(!is_active[e_it->n1] || !is_active[e_it->n2])
				{
					g.deactivate_edge(e_it);
					continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
				}

				max_edge_of_node[e_it->n1] = dummy_edge;
				max_edge_of_node[e_it->n2] = dummy_edge;

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

				std::cerr << "Not a correct matching, " << *e_it << " is adjacent to another matched edge." << std::endl;
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
			result += it->weight;//g.get_edge(*it).weight;
		}

		return result;
	}

};

}

