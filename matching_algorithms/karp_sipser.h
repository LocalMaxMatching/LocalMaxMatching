#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <boost/functional/hash.hpp>


namespace dipl
{

template <class Graph>
class KarpSipser
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;
	typedef typename Graph::Edge	   Edge;

public:

	static void deactivate_edge(Graph &g, EdgeIDType edge_id, std::list<NodeIDType> &modified_nodes)
	{
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
									     std::list<typename Graph::Edge> &matching,
									     unsigned int &depth)
	{
		depth = 0;
		typedef typename Graph::edge_id_iterator edge_iterator;

		// stores for each vertex the incident edge id, of the edge with the largest weight
		std::vector<NodeIDType> deg_one_vertices;
		std::vector<NodeIDType > pos_in_deg_one_vertices(g.num_vertices(), -1);

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
			if(deg_one_vertices.size()>0)
			{ // there's a degree 1 vertex
				NodeIDType rand_node_pos = rand()%deg_one_vertices.size();
				EdgeIDType edge_id = *g.begin_incident_active_edges(deg_one_vertices[rand_node_pos]);

				typename Graph::Edge current_edge = g.get_edge(edge_id);

				matching.push_back(current_edge);

				std::list<NodeIDType> modified_nodes;

				deactivate_edge(g, edge_id, modified_nodes);

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
			else
			{
				EdgeIDType edge_id = g.random_active_edge();

				typename Graph::Edge current_edge = g.get_edge(edge_id);

				matching.push_back(current_edge);

				std::list<NodeIDType> modified_nodes;

				deactivate_edge(g, edge_id, modified_nodes);

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

