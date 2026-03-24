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

public:




	static bool better_than_old_edge(const Graph &g,
									 const WeightType old_weight, const EdgeIDType old_edge_id,
									 const WeightType new_weight, const EdgeIDType new_edge_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash2(old_edge_id);
		unsigned int new_hash = hash2(new_edge_id);

//		std::cout << "old: " << old_edge_id << " -> " << old_hash << ",  new: " << new_edge_id << " -> " << new_hash << ",   (" << old_weight << ", " << new_weight << ")" << std::endl;

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
									     std::list<typename Graph::EdgeIDType> &matching,
									     unsigned int &depth
#ifdef MORE_INFORMATION
									     ,
									     std::list<typename Graph::EdgeIDType> &active_edge_count
#endif
									    )
	{
		typedef typename Graph::edge_iterator edge_iterator;

		std::vector<bool> is_active(g.num_vertices()+1, true);
		is_active[g.num_vertices()] = false;

		// stores for each vertex the incident edge id, of the edge with the largest weight
		std::vector<EdgeIDType > max_edge_of_node(g.num_vertices(), g.num_edges()); // initialize with dummy edge id, this edge is always inactive

		depth = 0;

		while(g.has_active_edges())
		{
#ifdef MORE_INFORMATION
			active_edge_count.push_back(g.num_active_edges());
#endif

			depth++;
//			std::cout << "round " << depth << std::endl;
//			std::cout << g.num_active_edges() << std::endl;

			// iterate over each active edge, to set local max edges
			for(edge_iterator e_it=g.begin_active_edges(); e_it<g.end_active_edges(); e_it++)
			{
//				std::cout << *e_it << std::endl;
				NodeIDType current_partner;
				WeightType current_weight;

				// i only use this comparison to get the same behavior as the parallel version
				if(max_edge_of_node[e_it->n1]<g.num_edges())
				{
					current_partner = g.get_edge(max_edge_of_node[e_it->n1]).n1==e_it->n1 ? g.get_edge(max_edge_of_node[e_it->n1]).n2 : g.get_edge(max_edge_of_node[e_it->n1]).n1;
					current_weight  = g.get_edge_weight(max_edge_of_node[e_it->n1]);
				}
				else
				{
					current_partner = 0;
					current_weight  = min_val<WeightType>();
				}


//				std::cout << "vertex " << e_it->n1 << ": ";
				// check first endpoint of current edge
				if(better_than_old_edge(g,
										current_weight, current_partner ^ e_it->n1,
										e_it->weight,                       e_it->n2 ^ e_it->n1) ||
				   !is_active[current_partner])
				{
					max_edge_of_node[e_it->n1] = g.get_edge_id(e_it);
				}


				if(max_edge_of_node[e_it->n2]<g.num_edges())
				{
					current_partner = g.get_edge(max_edge_of_node[e_it->n2]).n1==e_it->n2 ? g.get_edge(max_edge_of_node[e_it->n2]).n2 : g.get_edge(max_edge_of_node[e_it->n2]).n1;
					current_weight  = g.get_edge_weight(max_edge_of_node[e_it->n2]);
				}
				else
				{
					current_partner = 0;
					current_weight  = min_val<WeightType>();
				}


//				std::cout << "vertex " << e_it->n2 << ": ";
				// check first endpoint of current edge
				if(better_than_old_edge(g,
										current_weight, current_partner ^ e_it->n2,
										e_it->weight,                       e_it->n1 ^ e_it->n2) ||
					!is_active[current_partner])
				{
					max_edge_of_node[e_it->n2] = g.get_edge_id(e_it);
				}

			}


			// add matching-edges to matchings-list and set them inactive, also set the incident nodes inactive
			edge_iterator e_it=g.begin_active_edges();
			while(e_it<g.end_active_edges())
			{
				// check first endpoint of current edge
				if(g.get_edge_id(e_it)==max_edge_of_node[e_it->n1] && g.get_edge_id(e_it)==max_edge_of_node[e_it->n2])
				{
//					std::cout << (*e_it) << std::endl;
					matching.push_back(g.get_edge_id(e_it));
					is_active[e_it->n1] = false;
					is_active[e_it->n2] = false;

//					g.deactivate_edge(e_it);
//					continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
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
	static bool is_maximal_matching(std::list<typename Graph::EdgeIDType> &matchings,
							 	    Graph &g,
							 	    NodeIDType &unmatched_nodes)
	{
		bool result = true;

		std::vector<bool> matched_edge(g.num_edges(), false);
		std::vector<EdgeIDType> incident_edge(g.num_vertices(), std::numeric_limits<EdgeIDType>::max());
		std::vector<bool> matched_node(g.num_vertices(), false);

		for(typename std::list<typename Graph::EdgeIDType>::const_iterator it=matchings.begin(); it!=matchings.end(); it++)
		{
			matched_node[g.get_edge(*it).n1] = true;
			matched_node[g.get_edge(*it).n2] = true;

			incident_edge[g.get_edge(*it).n1] = *it;
			incident_edge[g.get_edge(*it).n2] = *it;

			matched_edge[*it] = true;
		}

		// check if the given matching is a correct matching
		for(typename Graph::edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
		{
			if(matched_edge[g.get_edge_id(e_it)] && (incident_edge[e_it->n1] != g.get_edge_id(e_it) || incident_edge[e_it->n2] != g.get_edge_id(e_it)) )
			{
				result = false;

				std::cerr << "Not a correct matching, " << g.get_edge(e_it) << " is adjacent to another matched edge." << std::endl;
			}
		}

		// check if the given matching is maximal
		for(typename Graph::edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
		{
			if(!matched_edge[g.get_edge_id(e_it)])
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


	static double get_weight(const Graph &g, const std::list<EdgeIDType> &matching)
	{
		double result = 0.;

		for(typename std::list<EdgeIDType>::const_iterator it=matching.begin(); it!=matching.end(); it++)
		{
			result += g.get_edge(*it).weight;
		}

		return result;
	}

};

}

