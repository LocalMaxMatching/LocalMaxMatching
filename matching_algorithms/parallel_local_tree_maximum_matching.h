#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <mpi.h>

#include <common/hash_functions.h>
#include <common/math_funcs.h>
#include <matching_algorithms/parallel_forest_matching.h>


namespace dipl
{

template <class Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType) = dipl::identity<typename Graph::NodeIDType> >
class ParallelLocalTreeMaximumMatching
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;

	typedef typename Graph::Edge	   		Edge;
	typedef typename Graph::edge_iterator 		edge_iterator;
	typedef typename Graph::const_edge_iterator const_edge_iterator;


	/// \todo maybe use weight instead of edge_id, why did i decide to use edge_id?
	/// have to do something like in 2nd loop of compute_local_tree_matchings
	struct Candidate
	{
		EdgeIDType edge_id;
		NodeIDType partner;

		Candidate()
		: edge_id(-1), partner(-1)
		{}

		Candidate(const EdgeIDType edge_id, const NodeIDType partner)
		: edge_id(edge_id), partner(partner)
		{}

		void set(const EdgeIDType edge_id, const NodeIDType partner)
		{
			this->edge_id = edge_id;
			this->partner = partner;
		}
	};

	static Candidate base_candidate;

	struct Message
	{
		NodeIDType src;
		NodeIDType candidate;

		Message()
		: src(-1), candidate(-1)
		{}

		Message(const NodeIDType src, const NodeIDType candidate)
		: src(src), candidate(candidate)
		{}
	};


	static MPI_Datatype CANDIDATE_MESSAGE_TYPE;
	static MPI_Datatype MATCHED_MESSAGE_TYPE;


	static const NodeIDType MATCHED       = -1;//std::numeric_limits<NodeIDType>::max();

	static const int msg_tag = 0;
	static const int matched_nodes_tag = 1;


	// buffers for messages that are used to send candidate requests of ghostnodes
	static std::vector< std::pair< std::vector<Message>, MPI_Request> >    candidate_messages_of_proc;


	// buffers for messages that are used to send information about matched ghostnodes
	static std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > matched_ghostnode_messages_of_proc;

	static long number_of_messages_sent, number_of_messages_received;

public:

	static bool better_than_old_partner(const Graph &g, const NodeIDType src_id,
										const WeightType old_weight, const NodeIDType old_partner_id,
										const WeightType new_weight, const NodeIDType new_partner_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash2(old_partner_id ^ src_id);
		unsigned int new_hash = hash2(new_partner_id ^ src_id);

		return (new_weight>old_weight) ||
			   (new_weight==old_weight && new_hash<old_hash) ||
			   (new_weight==old_weight && new_hash==old_hash && (new_partner_id^src_id)<(old_partner_id^src_id));
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
									     std::vector<Edge> &matching,
									     unsigned int &rounds )
	{
		int procs, proc_id;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

		init();

		dipl::ParallelForestMatching<Graph> forest_matching(g);


		/// \todo use dummy node with ID g.num_all_local_vertices(), allows us to initialize candidates with partner dummy node
		std::vector<bool> matched(g.num_all_local_vertices(), false);

		// let the base candidate point to the dummy edge and dummy node -> don't correspond to actual nodes or edge
		// weight of dummy edge at pos g.num_all_local_edges() is min_val(WeightType)
		base_candidate = Candidate(g.num_all_local_edges(), g.num_all_local_vertices());

		// stores for each vertex the incident edge id, of the edge with the largest weight
		std::vector<Candidate > candidate_of_node(g.num_all_local_vertices(), base_candidate); // initialize with dummy candidate

		// stores for each round a list of nodes that are incident to cross edges that aren't part of
		// parallel tree of this round.
		// Those nodes are ghostnodes on another process and might be matched in this round. Thus we
		// have to inform the other process
		std::vector<NodeIDType> matched_nodes_candidates;

		candidate_messages_of_proc.resize(procs);
		matched_ghostnode_messages_of_proc.resize(procs);

		rounds = 0;

#ifdef MORE_INFORMATION
		double start, end;
#endif

		while(g.has_active_edges())
		{
			rounds++;

#ifdef MORE_INFORMATION
			start = MPI_Wtime();
#endif

			set_maximal_adjacent_nodes(g, candidate_of_node);

#ifdef MORE_INFORMATION
			end = MPI_Wtime();
			forest_matching.outfile_std_out << "proc " << proc_id << ": time to run set_maximal_adjacent_nodes=" << (end-start) << " sec" << std::endl;

			start = MPI_Wtime();
#endif

			compute_local_tree_matchings(g, candidate_of_node, forest_matching, matching, matched, rounds);

#ifdef MORE_INFORMATION
			end = MPI_Wtime();

			forest_matching.outfile_std_out << "proc " << proc_id << ": time to run compute_local_tree_matchings=" << (end-start) << " sec" << std::endl;

			start = MPI_Wtime();
#endif

			deactivate_edges(g, matched, candidate_of_node);

#ifdef MORE_INFORMATION
			end = MPI_Wtime();

			forest_matching.outfile_std_out << "proc " << proc_id << ": time to run deactivate_edges=" << (end-start) << " sec" << std::endl;
#endif
		}

#ifdef MORE_INFORMATION
		forest_matching.outfile_std_out << "proc " << proc_id << ": number_of_mpi_messages_sent=" << (forest_matching.number_of_mpi_messages_sent + number_of_messages_sent) << std::endl;
		forest_matching.outfile_std_out << "proc " << proc_id << ": number_of_mpi_messages_received=" << (forest_matching.number_of_mpi_messages_received + number_of_messages_received) << std::endl;
		forest_matching.outfile_std_out << "proc " << proc_id << ": calls to gid to lid of graph instead of tree=" << (forest_matching.calls_to_g_for_gid_to_lid) << std::endl;
#endif

	}


	static void compute_local_tree_matchings(
						const Graph &g,
						std::vector<Candidate> &candidate_of_node,
						dipl::ParallelForestMatching<Graph> &forest_matching,
						std::vector<Edge> &matching,
						std::vector<bool> &matched,
						int depth)
	{
		int procs, proc_id;
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
		MPI_Comm_size(MPI_COMM_WORLD, &procs);


		std::vector<EdgeIDType> tree_local_edges; // use a counter somewhere for the correct size, if possible
		std::vector<EdgeIDType> tree_cross_edges; // use a counter somewhere for the correct size, if possible
		std::vector<EdgeIDType> tree_matching;

		// get tree edges, those are the incident edges that are maximal
		for(const_edge_iterator e_it=g.begin_active_local_edges(); e_it<g.end_active_local_edges(); e_it++)
		{
			if(g.get_edge_id(e_it)==candidate_of_node[e_it->n1].edge_id)
			{
				tree_local_edges.push_back(candidate_of_node[e_it->n1].edge_id);
			}
			else if(g.get_edge_id(e_it)==candidate_of_node[e_it->n2].edge_id)
			{
				tree_local_edges.push_back(candidate_of_node[e_it->n2].edge_id);
			}
		}

		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(); e_it<g.end_active_local_cross_edges(); e_it++)
		{
			if(g.get_edge_id(e_it)==candidate_of_node[e_it->n1].edge_id)
			{
				tree_cross_edges.push_back(candidate_of_node[e_it->n1].edge_id);
			}
			else if(e_it->n1==candidate_of_node[e_it->n2].partner)
			{
				tree_cross_edges.push_back(g.get_edge_id(e_it));
			}
		}

		// get matched edges
		std::vector<NodeIDType> matched_nodes;

		const EdgeIDType old_matching_size = matching.size();

		forest_matching.get_matching(g, tree_local_edges, tree_cross_edges,
						  	   	   	 matching, matched_nodes);

		set_matched_nodes(g, candidate_of_node, old_matching_size, matching, matched_nodes, matched);
	}


	static void set_matched_nodes(const Graph &g,
								  const std::vector<Candidate> &candidate_of_node,
								  const EdgeIDType old_matching_size,
								  std::vector<Edge> &matching,
								  std::vector<NodeIDType> &matched_nodes,
								  std::vector<bool> &matched)
	{
		// set newly matched nodes
		for(EdgeIDType pos=old_matching_size; pos<matching.size(); pos++)
		{
			matched[matching[pos].n1] = true;
			matched[matching[pos].n2] = true;

			// change local node ids to global ids
			matching[pos].n1 = g.local_vertex_id_to_global_id(matching[pos].n1);
			matching[pos].n2 = g.local_vertex_id_to_global_id(matching[pos].n2);
		}

		for(NodeIDType pos=0; pos<matched_nodes.size(); pos++)
		{
			matched[matched_nodes[pos]] = true;
		}


		// send information about matched nodes incident to non parallel tree cross edges
		send_matched_ghostnodes_incident_to_non_parallel_tree_edges(g, candidate_of_node, matched);

		// receive information about matched ghostnodes
		receive_matched_ghostnodes(g, matched);

		clear_messages(matched_ghostnode_messages_of_proc);
	}



	static void set_maximal_adjacent_nodes(const Graph &g, std::vector<Candidate> &candidate_of_node)
	{

		set_local_candidates(g, candidate_of_node);

		// get the candidates of the ghost nodes
		send_candidates_of_ghostnodes(g, candidate_of_node);
		receive_candidates_of_ghostnodes(g, candidate_of_node);

		clear_messages(candidate_messages_of_proc);
	}


	static void set_local_candidates(const Graph &g, std::vector<Candidate> &candidate_of_node)
	{
		// iterate over each active local edge, to set local max edges
		for(const_edge_iterator e_it=g.begin_active_local_edges(); e_it<g.end_active_local_edges(); e_it++)
		{
//			const NodeIDType adjusted_e_it_n1 = adjust_id(g.local_vertex_id_to_global_id(e_it->n1));
//			const NodeIDType adjusted_e_it_n2 = adjust_id(g.local_vertex_id_to_global_id(e_it->n2));

			/// \todo maybe store the global vertex IDs in the candidate objects?
			// check first endpoint of current edge
//			if(better_than_old_partner(g, adjusted_e_it_n1,
//									   g.get_edge(candidate_of_node[e_it->n1].edge_id).weight,
//									   adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n1].partner)),
//									   e_it->weight,
//									   adjusted_e_it_n2 ) )
			if(dipl::vertex_tie_breaking(g, *e_it,g.get_edge(candidate_of_node[e_it->n1].edge_id)))
			{
				candidate_of_node[e_it->n1].set(g.get_edge_id(e_it), e_it->n2);
			}

			// check second endpoint of current edge
//			if(better_than_old_partner(g, adjusted_e_it_n2,
//									   g.get_edge(candidate_of_node[e_it->n2].edge_id).weight,
//									   adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n2].partner)),
//									   e_it->weight,
//									   adjusted_e_it_n1 ) )
			if(dipl::vertex_tie_breaking(g, *e_it,g.get_edge(candidate_of_node[e_it->n2].edge_id)))
			{
				candidate_of_node[e_it->n2].set(g.get_edge_id(e_it), e_it->n1);
			}

		}


		// iterate over each active cross edge, to set local max edges
		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(); e_it<g.end_active_local_cross_edges(); e_it++)
		{
//			const NodeIDType adjusted_e_it_n1 = adjust_id(g.local_vertex_id_to_global_id(e_it->n1));
//			const NodeIDType adjusted_e_it_n2 = adjust_id(g.local_vertex_id_to_global_id(e_it->n2));

			// we're only interested in the best partner of local nodes, we set the local partners of ghostnodes in the next step
			// n1 is always local and n2 is ghost

			/// \todo maybe store the global vertex IDs in the candidate objects?
			// check first endpoint of current edge
//			if(better_than_old_partner(g, adjusted_e_it_n1,
//									   g.get_edge(candidate_of_node[e_it->n1].edge_id).weight,
//									   adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n1].partner)),
//									   e_it->weight,
//									   adjusted_e_it_n2 ) )
			if(dipl::vertex_tie_breaking(g, *e_it,g.get_edge(candidate_of_node[e_it->n1].edge_id)))
			{
				candidate_of_node[e_it->n1].set(g.get_edge_id(e_it), e_it->n2);
			}
		}
	}


	static void send_candidates_of_ghostnodes(const Graph &g, std::vector<Candidate> &candidate_of_node)
	{
		// if a local node wants to match with a ghost node then send this information
		// to the process of this ghost node!!
		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(); e_it<g.end_active_local_cross_edges(); e_it++)
		{
			const Edge &e = *e_it;
			// n1 is always local and n2 is ghost
			if(e.n2==candidate_of_node[e.n1].partner)
			{// we might send unnecessary messages if there are multi-edges with the same weight
				candidate_messages_of_proc[g.get_proc_of_ghost_vertex(e.n2)].first.push_back(Message(g.local_vertex_id_to_global_id(e.n1),
																								 	 g.local_vertex_id_to_global_id(e.n2)));
			}
		}


		send_messages(g, candidate_messages_of_proc);
	}


	/// \todo implement IRecv version of this function
	static void receive_candidates_of_ghostnodes(const Graph &g, std::vector<Candidate> &candidate_of_node)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Status msg_status;


		// receive matched messages
		for(int p=0; p<procs; p++)
		{
			if(!g.is_active_partner(p))
			{// if we send a message to p, we also expect a response from this process, even if we no longer have any cross-edges this process
				continue;
			}

			// only receive messages from active partners

			/// \todo maybe use ANY_SOURCE instead (don't wait for a specific message). In this case we need a loop the iterates over all "incoming messages"
			MPI_Probe(p, msg_tag, MPI_COMM_WORLD, &msg_status);

			int msg_count;
			MPI_Get_count(&msg_status, CANDIDATE_MESSAGE_TYPE, &msg_count);//g.get_active_cross_edges_count_of_proc(p);

			Message *messages = new Message[msg_count];

			MPI_Recv(messages, msg_count, CANDIDATE_MESSAGE_TYPE, p, msg_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			number_of_messages_received++;


			for(int i=0; i<msg_count; i++)
			{
				const NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(messages[i].src); /// \todo not nice, get_local_ghost_id is expensive
				const NodeIDType local_id		= g.global_vertex_id_to_local_id_of_local_vertex(messages[i].candidate);

				// the weight isn't interesting in this case
				candidate_of_node[local_ghost_id].set(-1, local_id);
			}

			delete[] messages;
		}
	}



	static void send_matched_ghostnodes_incident_to_non_parallel_tree_edges(const Graph &g,
																			const std::vector<Candidate> &candidate_of_node,
																			std::vector<bool> &matched)
	{
		// if a local node wants to match with a ghost node then send this information
		// to the process of this ghost node!!
		/// \todo maybe it's fast to use a list candidate nodes (well pairs, or the edge itself)
		/// could be filled in send_candidates_of_ghostnodes -> we would reduce the number of runs through
		/// active cross edges
		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(); e_it<g.end_active_local_cross_edges(); e_it++)
		{
			const Edge &e = *e_it;
			// n1 is always local and n2 is ghost
			if(e.n2!=candidate_of_node[e.n1].partner)
			{
				// n1 is a ghostnode on the other process
				if(matched[e.n1])
				{
					matched_ghostnode_messages_of_proc[g.get_proc_of_ghost_vertex(e.n2)].first.push_back(g.local_vertex_id_to_global_id(e.n1));
				}
			}
		}


		send_messages(g, matched_ghostnode_messages_of_proc);
	}


	/// \todo implement IRecv version of this function
	static void receive_matched_ghostnodes(const Graph &g,
										   std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Status msg_status;


		// receive matched messages
		for(int p=0; p<procs; p++)
		{
			if(!g.is_active_partner(p))
			{// if we send a message to p, we also expect a response from this process, even if we no longer have any cross-edges this process
				continue;
			}

			// only receive messages from active partners

			/// \todo maybe use ANY_SOURCE instead (don't wait for a specific message). In this case we need a loop the iterates over all "incoming messages"
			MPI_Probe(p, matched_nodes_tag, MPI_COMM_WORLD, &msg_status);

			int msg_count;
			MPI_Get_count(&msg_status, MATCHED_MESSAGE_TYPE, &msg_count);//g.get_active_cross_edges_count_of_proc(p);

			NodeIDType *messages = new NodeIDType[msg_count];

			MPI_Recv(messages, msg_count, MATCHED_MESSAGE_TYPE, p, matched_nodes_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			number_of_messages_received++;

			for(int i=0; i<msg_count; i++)
			{
				const NodeIDType matched_ghost_node = g.global_vertex_id_to_local_id_of_ghost_vertex(messages[i]); /// \todo not nice, get_local_ghost_id is expensive

				matched[matched_ghost_node] = true;
			}

			delete[] messages;
		}
	}



	static void send_messages(const Graph &g, std::vector< std::pair< std::vector<Message>, MPI_Request> > &messages_of_proc)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		//for each active partner: send messages;
		/// \todo think about using a list with active partners
		for(int p=0; p<procs; p++)
		{
			if(!messages_of_proc[p].first.empty())
			{// check for messages_of_proc[p].size() == g.get_active_cross_edges_of_proc(p)??

				MPI_Isend(&messages_of_proc[p].first[0], messages_of_proc[p].first.size(), CANDIDATE_MESSAGE_TYPE, p, msg_tag, MPI_COMM_WORLD, &messages_of_proc[p].second);
				number_of_messages_sent++;
			}
			else if(g.is_active_partner(p))
			{
				// we have to send empty message to active partner, because they expect a message from each active partner
				// even if this message is empty!!
				MPI_Isend(0, 0, CANDIDATE_MESSAGE_TYPE, p, msg_tag, MPI_COMM_WORLD, &messages_of_proc[p].second);
				number_of_messages_sent++;
			}
		}
	}


	static void send_messages(const Graph &g, std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > &messages_of_proc)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		//for each active partner: send messages;
		/// \todo think about using a list with active partners
		for(int p=0; p<procs; p++)
		{
			if(!messages_of_proc[p].first.empty())
			{// check for messages_of_proc[p].size() == g.get_active_cross_edges_of_proc(p)??

				MPI_Isend(&messages_of_proc[p].first[0], messages_of_proc[p].first.size(), MATCHED_MESSAGE_TYPE, p, matched_nodes_tag, MPI_COMM_WORLD, &messages_of_proc[p].second);
				number_of_messages_sent++;
			}
			else if(g.is_active_partner(p))
			{
				// we have to send empty message to active partner, because they expect a message from each active partner
				// even if this message is empty!!
				MPI_Isend(0, 0, MATCHED_MESSAGE_TYPE, p, matched_nodes_tag, MPI_COMM_WORLD, &messages_of_proc[p].second);
				number_of_messages_sent++;
			}
		}
	}



	// forces some kind of synchronization between communicating processes
	/*!
	 *
	 * @param messages_of_proc
	 */
	template <typename MTYPE>
	static void clear_messages(std::vector< std::pair< std::vector<MTYPE>, MPI_Request> > &messages_of_proc)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		// clear all message vectors
		// we might have send a few more messages but all of them were empty, thus we don't have to clear
		// messages_of_proc[p].first, because it's already empty.
		for(int p=0; p<procs; p++)
		{
			if(!messages_of_proc[p].first.empty())
			{
				MPI_Wait(&messages_of_proc[p].second, MPI_STATUS_IGNORE);
				messages_of_proc[p].first.clear(); // probably doesn't free any allocated memory, but in this case that's not too bad.
												   // I'm not totally sure about the definition of clear
			}
		}
	}


	// deactivates edges incident to matched nodes
	static void deactivate_edges(Graph &g, const std::vector<bool> &matched, std::vector<Candidate > &candidate_of_node)
	{
		Candidate base_candidate = Candidate(g.num_all_local_edges(), g.num_all_local_vertices());


		// deactivate local edges that are incident to matched nodes
		edge_iterator e_it=g.begin_active_local_edges();
		while(e_it<g.end_active_local_edges())
		{
			// check first endpoint of current edge
			if(matched[e_it->n1] || matched[e_it->n2])
			{
				g.deactivate_local_edge(e_it);
				continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
			}

			// reset the candidate for the nodes that aren't matched
			// we might miss a few nodes, but those are unmatched nodes,
			// that are no longer incident to an active edge
			candidate_of_node[e_it->n1] = base_candidate;
			candidate_of_node[e_it->n2] = base_candidate;

			e_it++;
		}

		// deactivate cross edges that are incident to matched nodes
		e_it=g.begin_active_local_cross_edges();
		while(e_it<g.end_active_local_cross_edges())
		{
			// check first endpoint of current edge
			if(matched[e_it->n1] || matched[e_it->n2])
			{
				g.deactivate_local_cross_edge(e_it);
				continue ; // we just set an edge to inactive, thus the current iterator position references a new active edge
			}

			// reset the candidate for the nodes that aren't matched
			// we might miss a few nodes, but those are unmatched nodes,
			// that are no longer incident to an active edge
			candidate_of_node[e_it->n1] = base_candidate;
			candidate_of_node[e_it->n2] = base_candidate;

			e_it++;
		}
	}


	static void init()
	{
		number_of_messages_sent = 0;
		number_of_messages_received = 0;
		MPI_Type_contiguous(sizeof(Message), MPI_BYTE, &CANDIDATE_MESSAGE_TYPE);
		MPI_Type_commit(&CANDIDATE_MESSAGE_TYPE);

		MPI_Type_contiguous(sizeof(NodeIDType), MPI_BYTE, &MATCHED_MESSAGE_TYPE);
		MPI_Type_commit(&MATCHED_MESSAGE_TYPE);
	}

	static double get_weight(const Graph &g, const std::vector<Edge> &matching)
	{
		double local_result = 0.;

		for(typename std::vector<Edge>::const_iterator it=matching.begin(); it!=matching.end(); it++)
		{
			local_result += it->weight;
		}

		double global_result;

		MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		return global_result;
	}

};


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
MPI_Datatype ParallelLocalTreeMaximumMatching<Graph, adjust_id>::CANDIDATE_MESSAGE_TYPE;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
MPI_Datatype ParallelLocalTreeMaximumMatching<Graph, adjust_id>::MATCHED_MESSAGE_TYPE;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::pair< std::vector<typename ParallelLocalTreeMaximumMatching<Graph, adjust_id>::Message>, MPI_Request> >
ParallelLocalTreeMaximumMatching<Graph, adjust_id>::candidate_messages_of_proc;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::pair< std::vector<typename Graph::NodeIDType>, MPI_Request> >
ParallelLocalTreeMaximumMatching<Graph, adjust_id>::matched_ghostnode_messages_of_proc;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
long
ParallelLocalTreeMaximumMatching<Graph, adjust_id>::number_of_messages_sent;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
long
ParallelLocalTreeMaximumMatching<Graph, adjust_id>::number_of_messages_received;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
typename ParallelLocalTreeMaximumMatching<Graph, adjust_id>::Candidate
ParallelLocalTreeMaximumMatching<Graph, adjust_id>::base_candidate;

}

