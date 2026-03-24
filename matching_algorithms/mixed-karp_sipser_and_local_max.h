#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <boost/functional/hash.hpp>
#include <common/hash_functions.h>


namespace dipl
{

// use AdjacencyArrayGraphWithActiveEdgesVector as the Graph

template <class Graph>
class MixedKarpSipserLocalMax
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;

	typedef typename Graph::Edge Edge;

public:

	/*!
	 * \brief Deactivates the provided edge and adjusts the degree-one management
	 * arrays. Adds new degree-one vertices and removes degree-zero vertices.
	 *
	 * @param g
	 * @param edge_id
	 * @param deg_one_vertices
	 * @param pos_in_deg_one_vertices
	 */
	static void deactivate_edge(Graph &g, EdgeIDType edge_id,
								std::vector<NodeIDType> &deg_one_vertices,
								std::vector<NodeIDType > &pos_in_deg_one_vertices)
	{
		std::list<NodeIDType> modified_nodes;
		typename std::list<EdgeIDType> adjacent_edges;
		typename Graph::Edge current_edge = g.get_edge(edge_id);

		modified_nodes.push_back(current_edge.n1);
		modified_nodes.push_back(current_edge.n2);

		g.deactivate_edge(edge_id);

		for(typename Graph::incident_edges_iterator it=g.begin_incident_active_edges(current_edge.n1);
				it!=g.end_incident_active_edges(current_edge.n1); it++)
		{
			modified_nodes.push_back(g.get_edge(*it).n1);
			modified_nodes.push_back(g.get_edge(*it).n2);
			adjacent_edges.push_back(*it);
		}

		for(typename Graph::incident_edges_iterator it=g.begin_incident_active_edges(current_edge.n2);
						it!=g.end_incident_active_edges(current_edge.n2); it++)
		{
			modified_nodes.push_back(g.get_edge(*it).n1);
			modified_nodes.push_back(g.get_edge(*it).n2);
			adjacent_edges.push_back(*it);
		}

		for(typename std::list<EdgeIDType>::iterator it=adjacent_edges.begin(); it!=adjacent_edges.end(); it++)
		{
			g.deactivate_edge(*it);
		}


		// get new degree-one vertices and remove degree-zero vertices
		for(typename std::list<NodeIDType>::iterator it=modified_nodes.begin(); it!=modified_nodes.end(); it++)
		{
			if(g.active_degree(*it)==1)
			{
				if(pos_in_deg_one_vertices[*it]>=deg_one_vertices.size() || deg_one_vertices[pos_in_deg_one_vertices[*it]]!=*it)
				{
					pos_in_deg_one_vertices[*it] = deg_one_vertices.size();
					deg_one_vertices.push_back(*it);
				}
			}
			else if(g.active_degree(*it)==0)
			{
				if(pos_in_deg_one_vertices[*it]<deg_one_vertices.size() && deg_one_vertices[pos_in_deg_one_vertices[*it]]==*it)
				{
					NodeIDType pos = pos_in_deg_one_vertices[*it];
					deg_one_vertices[pos] = deg_one_vertices.back();
					pos_in_deg_one_vertices[deg_one_vertices[pos]] = pos;
					deg_one_vertices.pop_back();
				}
			}
		}

	}


	static bool better_than_old_edge(const WeightType old_weight, const EdgeIDType old_edge_id,
									 const WeightType new_weight, const EdgeIDType new_edge_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash2(old_edge_id);
		unsigned int new_hash = hash2(new_edge_id);

		return (new_weight>old_weight) ||
			   (new_weight==old_weight && new_hash<old_hash) ||
			   (new_weight==old_weight && new_hash==old_hash && new_edge_id<old_edge_id);
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
	 */
	static void add_local_max_edges_to_matchings(const std::list<EdgeIDType> &matching_candidates,
											  	 Graph &g,
											  	 const std::vector<EdgeIDType> &max_edge_of_node,
											  	 std::vector<bool> &is_active,
											  	 std::list<Edge> &matching,
											  	 std::vector<NodeIDType> &deg_one_vertices,
											  	 std::vector<NodeIDType > &pos_in_deg_one_vertices)
	{
		for(typename std::list<EdgeIDType>::const_iterator candidate_it=matching_candidates.begin(); candidate_it!=matching_candidates.end(); candidate_it++)
		{
			EdgeIDType edge_id = *candidate_it;
			typename Graph::Edge edge = g.get_edge(edge_id);

			// the following condition is ok, because we guaranteed that max_edge_of_node stores the smallest edge-id of incident edges with max weight
			if(max_edge_of_node[edge.n1]==edge_id && max_edge_of_node[edge.n2]==edge_id &&
			   g.is_active(edge_id) /*we don't want to add an edge more than twice to the matching list*/)
			{
				deactivate_edge(g, edge_id, deg_one_vertices, pos_in_deg_one_vertices); // important, otherwise we'll check those edges over and over again
				matching.push_back(g.get_edge(edge_id));
				is_active[edge.n1] = false;
				is_active[edge.n2] = false;
			}
		}

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
	static void compute_approximate_maximum_matching(Graph &g,
									     std::list<Edge> &matching,
									     unsigned int &rounds)
	{
		rounds = 0;
		typedef typename Graph::edge_id_iterator edge_iterator;

		std::vector<NodeIDType> deg_one_vertices;
		std::vector<NodeIDType > pos_in_deg_one_vertices(g.num_vertices(), std::numeric_limits<NodeIDType>::max());

		// stores for each vertex the edge-ID of dominating incident edge
		std::vector<EdgeIDType> max_edge_of_node(g.num_vertices(), 0);


		std::vector<bool> is_active(g.num_vertices(), true);

		std::list<NodeIDType> active_nodes;
		for(NodeIDType i=0; i<g.num_vertices(); i++)
		{
			active_nodes.push_back(i);
		}

		// get degree one vertices
		for(edge_iterator e_it=g.begin_edges(); e_it!=g.end_edges(); e_it++)
		{
			const Edge &current_edge = g.get_edge(e_it);
			if(g.active_degree(current_edge.n1)==1)
			{
				pos_in_deg_one_vertices[current_edge.n1] = deg_one_vertices.size();
				deg_one_vertices.push_back(current_edge.n1);
			}

			if(g.active_degree(current_edge.n2)==1)
			{
				pos_in_deg_one_vertices[current_edge.n2] = deg_one_vertices.size();
				deg_one_vertices.push_back(current_edge.n2);
			}
		}


		while(g.has_active_edges())
		{
			rounds++;

			// karp-sipser: add degree-one vertices to matching
			while(deg_one_vertices.size()>0)
			{ // there's a degree 1 vertex
				NodeIDType rand_node_pos = rand()%deg_one_vertices.size();
				EdgeIDType edge_id = *g.begin_incident_active_edges(deg_one_vertices[rand_node_pos]);

				typename Graph::Edge current_edge = g.get_edge(edge_id);

				matching.push_back(current_edge);

				// also adjusts the degree-one management vectors
				deactivate_edge(g, edge_id, deg_one_vertices, pos_in_deg_one_vertices);
			}


			// use local max algo to reduce the graph

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

			add_local_max_edges_to_matchings(matching_candidates, g, max_edge_of_node, is_active, matching,
											 deg_one_vertices, pos_in_deg_one_vertices);

			remove_inactive_nodes(active_nodes, is_active);

		}

	}


	/*!
	 * This name is confusing since it doesn't compute weighted matchings.
	 * It's just a wrapper function for compute_approximate_maximum_matching.
	 * I've added this function compatibility reasons.
	 *
	 * @param g
	 * @param matching
	 * @param depth
	 */
	static void compute_weighted_matching(Graph &g,
										  std::list<typename Graph::Edge> &matching,
										  unsigned int &depth)
	{
		compute_approximate_maximum_matching(g,
											 matching,
											 depth);
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
	static bool is_maximal_matching(const std::list<Edge> &matchings,
								 	    const Graph &g,
								 	    NodeIDType &unmatched_nodes)
	{
		typedef typename Graph::incident_edges_iterator incident_edges_iterator;

		bool result = true;

		std::vector<bool> matched_node(g.num_vertices(), false);

		for(typename std::list<Edge>::const_iterator it=matchings.begin(); it!=matchings.end(); it++)
		{
			if(matched_node[it->n1])
			{// not a correct matching, 2 matched edges are adjacent
				result = false;

				std::cout << "not a correct matching, 2 matched edges are adjacent at node " << it->n1 << std::endl;
			}

			if(matched_node[it->n2])
			{// not a correct matching, 2 matched edges are adjacent
				result = false;

				std::cout << "not a correct matching, 2 matched edges are adjacent at node " << it->n2 << std::endl;
			}

			matched_node[(*it).n1] = true;
			matched_node[(*it).n2] = true;
		}


		// check if the given matching is maximal
		for(EdgeIDType e=0; e<g.num_edges(); e++)
		{
			if(!matched_node[g.get_edge(e).n1] && !matched_node[g.get_edge(e).n2])
			{
				result = false;

				std::cout << "Matching isn't maximal, end vertices of edge " << g.get_edge(e) << " haven't been matched!" << std::endl;
			}
		}

		unmatched_nodes = g.num_vertices();
		for(NodeIDType i=0; i<g.num_vertices(); i++)
		{
			unmatched_nodes -= matched_node[i];
		}



		return result; // correct maximal matching
	}


	static double get_weight(const Graph &g, const std::list<Edge> &matching)
	{
		double result = 0.;

		for(typename std::list<Edge>::const_iterator it=matching.begin(); it!=matching.end(); it++)
		{
			result += (*it).weight;
		}

		return result;
	}

};

}

