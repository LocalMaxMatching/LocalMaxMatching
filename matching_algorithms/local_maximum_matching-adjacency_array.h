#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <boost/functional/hash.hpp>

#include <common/hash_functions.h>


namespace dipl
{

template <class Graph>
class LocalMaximumMatching_AdjacencyArray
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;


public:

	static bool better_than_old_edge(const WeightType old_weight, const EdgeIDType old_edge_id,
									 const WeightType new_weight, const EdgeIDType new_edge_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash2(old_edge_id);
		unsigned int new_hash = hash2(new_edge_id);

		return (new_weight>old_weight) ||
			   (new_weight==old_weight && new_hash<old_hash) ||
			   (new_weight==old_weight && new_hash==old_hash && new_edge_id<old_edge_id);

//		return (new_weight>old_weight) ||
//			   (new_weight==old_weight && boost::hash_value(new_edge_id)<boost::hash_value(old_edge_id)) ||
//			   (new_weight==old_weight && boost::hash_value(new_edge_id)==boost::hash_value(old_edge_id) && new_edge_id<old_edge_id);
	}

	/*!
	 * \brief Gets the maximal active edge incident to the given node.
	 *
	 * Sets max_edge to the edge with the largest weight.
	 * If there are several edges with the same max weight then max_edge is set
	 * to the edge with the smallest ID of those edges!!
	 *
	 * Runtime: O(deg_active(node_id))
	 *
	 * @param node_id	The given node.
	 * @param g			The graph to be used.
	 * @param max_edge	A reference to the maximal incident active edge.
	 * @return 			Returns true if a maximal edge was found otherwise false.
	 */
	static bool get_maximal_incident_edge(const NodeIDType node_id,
										   const Graph &g,
										   EdgeIDType &max_edge)
	{
		typedef typename Graph::incident_edges_iterator incident_edges_iterator;

		WeightType current_max_weight = -std::numeric_limits<WeightType>::max();
		bool new_matching_edge = false;

		for(incident_edges_iterator incident_edge_it=g.begin_incident_active_edges(node_id);
			incident_edge_it!=g.end_incident_active_edges(node_id);
			incident_edge_it++)
		{
//			if( g.get_edge(incident_edge_it).weight>current_max_weight ||
//			   (g.get_edge(incident_edge_it).weight==current_max_weight && *incident_edge_it<max_edge))
			if(better_than_old_edge(current_max_weight,                  max_edge,
									g.get_edge(incident_edge_it).weight, *incident_edge_it))
			{
				max_edge = *incident_edge_it;
				current_max_weight = g.get_edge(max_edge).weight;
				new_matching_edge = true;
			}
		}

		return new_matching_edge;
	}


	/*!
	 * \brief Deactivates the given edge and all adjacent edges
	 *
	 * Runtime: O(1) + O(deg_active(end_vertex_1) + O(deg_active(end_vertex_2))
	 *
	 * @param edge_id	The edge to be deactivated
	 * @param g			The graph the edge belongs to
	 */
	static void deactivate_edge_and_adjacent_edges(const EdgeIDType edge_id,
											Graph &g)
	{
		typedef typename Graph::incident_edges_iterator incident_edges_iterator;
		std::list<EdgeIDType> edge_id_list;

		typename Graph::Edge edge = g.get_edge(edge_id);

		// get all edges adjacent to the given edge on end vertex n1
		// we can't directly deactivate those edges
		for(incident_edges_iterator incident_edge_it=g.begin_incident_active_edges(edge.n1); incident_edge_it!=g.end_incident_active_edges(edge.n1); incident_edge_it++)
		{
			edge_id_list.push_back(*incident_edge_it);
		}

		// get all edges adjacent to the given edge on end vertex n2
		// we can't directly deactivate those edges
		for(incident_edges_iterator incident_edge_it=g.begin_incident_active_edges(edge.n2); incident_edge_it!=g.end_incident_active_edges(edge.n2); incident_edge_it++)
		{
			edge_id_list.push_back(*incident_edge_it);
		}

		// deactivate the adjacent edges
		for(typename std::list<EdgeIDType>::iterator edge_it=edge_id_list.begin(); edge_it!=edge_id_list.end(); edge_it++)
		{
			g.deactivate_edge(*edge_it);
		}

		// deactivate the given edge
		g.deactivate_edge(edge_id);

	}


	/*!
	 * \brief Computes the locale maximal edges from the matching-candidates-list
	 * and adds them the matching list.
	 *
	 * Runtime: O(lenght_of(matching_candidates))
	 *
	 * @param matching_candidates	List of candidates for the resulting matching
	 * @param g						The corresponding graph
	 * @param max_edge_of_node		A vector that stores for each end-vertex of matching-candidates
	 * 								the ID of the maximal incident edge. If there's more than one
	 * 								incident maximal edge, the smallest ID of those edges is stored.
	 * @param is_active				A vector that specifies for each vertex if it's active.
	 * @param matching				The list local maximal edges are added to.
	 *
	 * @return						Returns true if an edge has been added to the matching-list.
	 */
	static bool add_local_max_edges_to_matchings(const std::list<EdgeIDType> &matching_candidates,
										  Graph &g,
										  const std::vector<EdgeIDType> &max_edge_of_node,
										  std::vector<bool> &is_active,
										  std::list<EdgeIDType> &matching)
	{
		bool made_a_change = false;

		for(typename std::list<EdgeIDType>::const_iterator candidate_it=matching_candidates.begin(); candidate_it!=matching_candidates.end(); candidate_it++)
		{
			EdgeIDType edge_id = *candidate_it;
			typename Graph::Edge edge = g.get_edge(edge_id);

			// the following condition is ok, because we guaranteed that max_edge_of_node stores the smallest edge-id of incident edges with max weight
			if(max_edge_of_node[edge.n1]==edge_id && max_edge_of_node[edge.n2]==edge_id &&
			   g.is_active(edge_id) /*we don't want to add an edge more than twice to the matching list*/)
			{
				deactivate_edge_and_adjacent_edges(edge_id, g); // important, otherwise we'll check those edges over and over again
				matching.push_back(edge_id);
				is_active[edge.n1] = false;
				is_active[edge.n2] = false;

				made_a_change = true;
			}
		}

		return made_a_change;
	}


	/*!
	 * \brief Removes nodes from the active_nodes-list, that are no longer active.
	 *
	 * Runtime: O(lenght_of(active_nodes))
	 *
	 * @param active_nodes	List of active nodes, that is updated
	 * @param is_active		Vector that specifies which nodes are active
	 */
	static void remove_inactive_nodes(std::list<NodeIDType> &active_nodes,
							   const std::vector<bool> &is_active)
	{
		typename std::list<NodeIDType>::iterator node_it = active_nodes.begin();

		while(node_it != active_nodes.end())
		{
			if(!is_active[*node_it])
			{
				node_it = active_nodes.erase(node_it);
			}
			else
			{
				node_it++;
			}
		}
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
								   std::list<EdgeIDType> &matching,
								   unsigned int &depth)
	{
		typedef typename Graph::incident_edges_iterator incident_edges_iterator;

		std::list<NodeIDType> active_nodes;
		for(NodeIDType i=0; i<g.num_vertices(); i++)
		{
			active_nodes.push_back(i);
		}

		std::vector<bool> is_active(g.num_vertices(), true);

		// stores for each vertex the incident edge id, of the edge with the largest weight
		std::vector<EdgeIDType> max_edge_of_node(g.num_vertices(), 0);

		bool changed = true;
		depth = 0;

		while(changed)
		{
			changed = false;
			depth++;

			std::list<EdgeIDType> matching_candidates;

			// iterate over each active node and check them for local max edges
			for(typename std::list<NodeIDType>::iterator node_it=active_nodes.begin(); node_it!=active_nodes.end(); node_it++)
			{
				EdgeIDType local_max_edge_candidate = 0;

				if(get_maximal_incident_edge(*node_it, g, local_max_edge_candidate))
				{
					matching_candidates.push_back(local_max_edge_candidate);
					max_edge_of_node[*node_it] = local_max_edge_candidate;
				}

			}

			changed = add_local_max_edges_to_matchings(matching_candidates, g, max_edge_of_node, is_active, matching);

			remove_inactive_nodes(active_nodes, is_active);
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
	static bool is_maximal_matching(const std::list<EdgeIDType> &matchings,
							 	    const Graph &g,
							 	    NodeIDType &unmatched_nodes)
	{
		typedef typename Graph::incident_edges_iterator incident_edges_iterator;

		bool result = true;

		std::vector<bool> matched_edge(g.num_edges(), false);
		std::vector<bool> matched_node(g.num_vertices(), false);

		for(typename std::list<EdgeIDType>::const_iterator it=matchings.begin(); it!=matchings.end(); it++)
		{
			matched_node[g.get_edge(*it).n1] = true;
			matched_node[g.get_edge(*it).n2] = true;

			matched_edge[*it] = true;
		}

		// check if the given matching is a correct matching
		for(typename std::list<EdgeIDType>::const_iterator it=matchings.begin(); it!=matchings.end(); it++)
		{
			for(incident_edges_iterator a_it=g.begin_incident_edges(g.get_edge(*it).n1); a_it< g.end_incident_edges(g.get_edge(*it).n1); a_it++)
			{
				if(*a_it!=*it && matched_edge[*a_it])
				{// given matching isn't a correct matching
					result = false;

					std::cout << "Matching not correct, contains edges: " << "(" << g.get_edge(*it).n1 << ", " << g.get_edge(*it).n2 << ", " << g.get_edge(*it).weight << ")"
							<< " and " << "(" << g.get_edge(*a_it).n1 << ", " << g.get_edge(*a_it).n2 << ", " << g.get_edge(*a_it).weight << ")" << std::endl;
				}
			}


			for(incident_edges_iterator a_it=g.begin_incident_edges(g.get_edge(*it).n2); a_it< g.end_incident_edges(g.get_edge(*it).n2); a_it++)
			{
				if(*a_it!=*it && matched_edge[*a_it])
				{// given matching isn't a correct matching
					result = false;

					std::cout << "Matching not correct, contains edges: " << "(" << g.get_edge(*it).n1 << ", " << g.get_edge(*it).n2 << ", " << g.get_edge(*it).weight << ")"
							<< " and " << "(" << g.get_edge(*a_it).n1 << ", " << g.get_edge(*a_it).n2 << ", " << g.get_edge(*a_it).weight << ")" << std::endl;
				}
			}
		}

		// check if the given matching is maximal
		for(EdgeIDType e=0; e<g.num_edges(); e++)
		{
			if(!matched_edge[e])
			{
				if(!matched_node[g.get_edge(e).n1] && !matched_node[g.get_edge(e).n2])
				{
					result = false;

					std::cout << "Matching isn't maximal, vertices " << g.get_edge(e).n1 << "  and  " << g.get_edge(e).n2 << " haven't been matched!" << std::endl;
				}
			}
		}

		unmatched_nodes = g.num_vertices();
		for(NodeIDType i=0; i<g.num_vertices(); i++)
		{
			unmatched_nodes -= matched_node[i];
		}



		return result; // correct maximal matching
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

