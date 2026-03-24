#ifndef _PARALLEL_LOCAL_MAXIMUM_MATCHING_H_
#define _PARALLEL_LOCAL_MAXIMUM_MATCHING_H_

#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <utility>

#include <mpi.h>

#include <boost/functional/hash.hpp>

#include <common/hash_functions.h>
#include <common/math_funcs.h>

namespace dipl
{





template <class Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType) = dipl::identity<typename Graph::NodeIDType> >
class ParallelLocalMaximumMatching
{
	typedef typename Graph::WeightType WeightType;
	typedef typename Graph::NodeIDType NodeIDType;
	typedef typename Graph::EdgeIDType EdgeIDType;
	typedef typename Graph::Edge	   Edge;
	typedef typename Graph::edge_iterator edge_iterator;
	typedef typename Graph::const_edge_iterator const_edge_iterator;


	struct Candidate
	{
		WeightType weight;
		NodeIDType partner;

		Candidate()
		: weight(min_val<WeightType>()),
		  partner(0) /// \todo this only works as long min_val<WeightType>() is smaller then any weight that might occur in graphs
		{}

		Candidate(const WeightType weight, const NodeIDType partner)
		: weight(weight), partner(partner)
		{}

		void set(const WeightType weight, const NodeIDType partner)
		{
			this->weight = weight;
			this->partner = partner;
		}
	};

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

	static const unsigned int msg_tag_candidates = 0;
	static const unsigned int msg_tag_ghost_vertices = 0x40000000;

	static unsigned int msgs_to_receive;
	static unsigned int round;

	// MSG-Buffers for IRecv-Version
	static std::vector<MPI_Request> recv_requests;
	static std::vector< std::vector<Message> > candidate_recv_messages;
	static std::vector< std::vector<NodeIDType> > matched_recv_messages;

	// buffers for messages that are used to send candidate requests of ghostnodes
	static std::vector< std::pair< std::vector<Message>, MPI_Request> >    candidate_messages_of_proc;
	// buffers for messages that are used to inform other processes about matched ghostnodes
	static std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > matched_messages_of_proc;

#ifdef MORE_INFORMATION
public:
	static double communication_time;
	static long long number_of_sent_messages;
	static long long sent_bytes;

	static unsigned int communication_factor;
#endif


public:




	static bool better_than_old_partner(const Graph &g, const NodeIDType src_id,
										const WeightType old_weight, const NodeIDType old_partner_id,
										const WeightType new_weight, const NodeIDType new_partner_id)
	{
		/// \todo think about the type of the hash, also the return type
		unsigned int old_hash = hash(old_partner_id ^ src_id);
		unsigned int new_hash = hash(new_partner_id ^ src_id);

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
									     unsigned int &rounds
									    )
	{
		init(); // e.g. inits MESSAGE_TYPE

		int procs, proc_id;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

//		set_recv_buffers(g);

		candidate_messages_of_proc.resize(procs, make_pair(std::vector<Message>(), MPI_Request()));
		matched_messages_of_proc.resize(procs, make_pair(std::vector<NodeIDType>(), MPI_Request()));

		/// \todo use dummy node with ID g.num_all_local_vertices(), allows us to initialize candidates with partner dummy node
		std::vector<bool> matched(g.num_all_local_vertices(), false);

		// stores for each vertex the incident edge id, of the edge with the largest weight
		std::vector<Candidate > candidate_of_node(g.num_all_local_vertices(), Candidate()); // initialize with dummy candidate

		rounds = 0;

		while(g.has_active_edges() /*&& rounds<4*/)
		{
			set_local_maximal_candidate_of_nodes(g, candidate_of_node, matched);

			set_matched_local_nodes_and_add_edge_to_matching(g, candidate_of_node, matched, matching);
			set_matched_ghost_nodes_and_add_edge_to_matching(g, candidate_of_node, matched, matching);

			deactivate_edges_incident_to_matched_nodes(g, matched, candidate_of_node);

			rounds++;
			round++;
			round %= msg_tag_ghost_vertices;// we don't want negative tags!!
		}


//		clear_recv_buffers();
	}


	// only process 0 gets the correct global result
	static double get_weight(const Graph &g, std::vector<Edge> &matching)
	{
		double local_result = 0.;
		double result = 0.;

		for(typename std::vector<Edge>::iterator it=matching.begin(); it!=matching.end(); it++)
		{
			local_result += it->weight;
		}

		MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		return result;
	}



private:

	static void set_local_maximal_candidate_of_nodes(Graph &g, std::vector<Candidate> &candidate_of_node, const std::vector<bool> &matched)
	{
		// iterate over each active local edge, to set local max edges
		for(edge_iterator e_it=g.begin_active_local_edges(),
						  it_end = g.end_active_local_edges();
			e_it<it_end; e_it++)
		{
			const NodeIDType adjusted_e_it_n1 = adjust_id(g.local_vertex_id_to_global_id(e_it->n1));
			const NodeIDType adjusted_e_it_n2 = adjust_id(g.local_vertex_id_to_global_id(e_it->n2));

			/// \todo maybe store the global vertex IDs in the candidate objects?
			// check first endpoint of current edge
			if(better_than_old_partner(g, adjusted_e_it_n1,
									   candidate_of_node[e_it->n1].weight, adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n1].partner)),
									   e_it->weight,                       adjusted_e_it_n2 ) )
			{
				candidate_of_node[e_it->n1].set(e_it->weight, e_it->n2);
			}

			// check second endpoint of current edge
			if(better_than_old_partner(g, adjusted_e_it_n2,
									   candidate_of_node[e_it->n2].weight, adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n2].partner)),
									   e_it->weight,                       adjusted_e_it_n1 ) )
			{
				candidate_of_node[e_it->n2].set(e_it->weight, e_it->n1);
			}

		}


		// iterate over each active local ghost edge, to set local max edges
		for(edge_iterator e_it=g.begin_active_local_cross_edges(),
						  it_end=g.end_active_local_cross_edges();
			e_it<it_end; e_it++)
		{// we have to set the candidates of both vertices -> we need this information to send the correct information

			const NodeIDType adjusted_e_it_n1 = adjust_id(g.local_vertex_id_to_global_id(e_it->n1));
			const NodeIDType adjusted_e_it_n2 = adjust_id(g.local_vertex_id_to_global_id(e_it->n2));

			// check first endpoint of current edge
			if(better_than_old_partner(g, adjusted_e_it_n1,
									   candidate_of_node[e_it->n1].weight, adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n1].partner)),
									   e_it->weight,                       adjusted_e_it_n2 ) )
			{
				candidate_of_node[e_it->n1].set(e_it->weight, e_it->n2);
			}

			// check second endpoint of current edge
			if(better_than_old_partner(g, adjusted_e_it_n2,
									   candidate_of_node[e_it->n2].weight, adjust_id(g.local_vertex_id_to_global_id(candidate_of_node[e_it->n2].partner)),
									   e_it->weight,                       adjusted_e_it_n1 ) )
			{
				candidate_of_node[e_it->n2].set(e_it->weight, e_it->n1);
			}

		}

	}


	static void set_matched_local_nodes_and_add_edge_to_matching(Graph &g,
																 std::vector<Candidate> &candidate_of_node,
																 std::vector<bool> &matched,
																 std::vector<Edge> &matching)
	{
		// add local matching edges to matchings list and set them inactive, also set the incident nodes inactive
		//
		// We have to check for local dominant edges before we handle dominant cross-edges. This makes it easy
		// to tell partner processes about locally matched nodes in the same round. Thus after each round every
		// process knows the correct amount of active cross edges.
		for(edge_iterator e_it=g.begin_active_local_edges(),
						  it_end=g.end_active_local_edges();
			e_it<it_end; e_it++)
		{
			// check first endpoint of current edge
			if(e_it->n2==candidate_of_node[e_it->n1].partner && e_it->n1==candidate_of_node[e_it->n2].partner &&
			   e_it->weight==candidate_of_node[e_it->n1].weight && !matched[e_it->n1] && !matched[e_it->n2])
			{ // have to check for matched, because of multi-graphs

				matching.push_back(Edge(g.local_vertex_id_to_global_id(e_it->n1),
										g.local_vertex_id_to_global_id(e_it->n2),
										e_it->weight));

				matched[e_it->n1] = true;
				matched[e_it->n2] = true;
			}
		}
	}


	static void set_matched_ghost_nodes_and_add_edge_to_matching(Graph &g,
																 std::vector<Candidate> &candidate_of_node,
																 std::vector<bool> &matched,
																 std::vector<Edge> &matching)
	{
		/*! \todo clear the message buffer at the end of each round, they cause some kind of
		 * synchronization, so there might less time waiting when doing other stuff
		 * before clearing those messages.
		 */
		get_ghost_candidates(g, candidate_of_node, matched, candidate_messages_of_proc);

#ifdef MORE_INFORMATION
		double start, end;
		start = MPI_Wtime();
#endif
		send_messages(g, candidate_messages_of_proc);

		receive_ghost_candidates(g, candidate_of_node, matched);
//		receive_ghost_candidates_irecv(g, candidate_of_node, matched);

//		exchange_candidates_all_to_all(g, candidate_messages_of_proc,
//									   candidate_of_node, matched);

		clear_messages(candidate_messages_of_proc);
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		// clear all message vectors
		// we might have send a few more messages but all of them were empty, thus we don't have to clear
		// messages_of_proc[p].first, because it's already empty.
		for(int p=0; p<procs; p++)
		{
			candidate_messages_of_proc[p].first.clear(); // probably doesn't free any allocated memory, but in this case that's not too bad.
												   // I'm not totally sure about the definition of clear
		}
#ifdef MORE_INFORMATION
		end = MPI_Wtime();
		communication_time += end-start;
#endif


		add_matched_cross_edges_to_matching_and_prepare_messages_for_ghost_nodes_adjacent_to_matched_local_nodes
														(	g,
															candidate_of_node,
															matched,
															matching,
															matched_messages_of_proc
														);

#ifdef MORE_INFORMATION
		start = MPI_Wtime();
#endif
//		exchange_ghost_vertices_all_to_all(g, matched_messages_of_proc, matched);

		send_messages(g, matched_messages_of_proc);

		receive_information_about_matched_ghost_nodes(g, matched);
//		receive_information_about_matched_ghost_nodes_irecv(g, matched);

		clear_messages(matched_messages_of_proc);
		for(int p=0; p<procs; p++)
		{
			matched_messages_of_proc[p].first.clear(); // probably doesn't free any allocated memory, but in this case that's not too bad.
												   // I'm not totally sure about the definition of clear
		}
#ifdef MORE_INFORMATION
		end = MPI_Wtime();
		communication_time += end-start;
#endif
	}


	static void get_ghost_candidates(const Graph &g,
									  const std::vector<Candidate> &candidate_of_node,
									  const std::vector<bool> &matched,
									  std::vector< std::pair< std::vector<Message>, MPI_Request> > &messages_of_proc)
	{
		int procs, proc_id;
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		// iterate ghost edges and send appropriate messages
		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(),
						  it_end=g.end_active_local_cross_edges();
			e_it<it_end; e_it++)
		{
			// n1 is always local -> n2 is ghost
			if(!matched[e_it->n1] &&
				e_it->n2==candidate_of_node[e_it->n1].partner &&
				e_it->n1==candidate_of_node[e_it->n2].partner &&
				e_it->weight==candidate_of_node[e_it->n1].weight) /// \todo weight comparison only necessary for multi-graphs
			{// we might send unnecessary messages if there are multi-edges with the same weight

#ifdef MORE_INFORMATION
				for(unsigned int i=0; i<communication_factor; i++)
#endif
				messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n2)]
				                 .first.push_back(Message(g.local_vertex_id_to_global_id(e_it->n1),
				                		 	 	 	 	  g.local_vertex_id_to_global_id(e_it->n2)));
			}

		}

	}

	static void exchange_candidates_all_to_all(const Graph &g,
									  	  	   const std::vector< std::pair< std::vector<Message>, MPI_Request> > &messages_of_proc,
									  	  	   std::vector<Candidate> &candidate_of_node,
									  	  	   std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		std::vector<int> counts(procs);
		std::vector<int> displs(procs);

		int total_number_of_messages = 0;

		for(int p=0; p<procs; p++)
		{
			int cross_edges_to_p = g.get_active_cross_edges_count_of_proc(p);
#ifdef MORE_INFORMATION
			cross_edges_to_p *= communication_factor;
#endif
			counts[p] = cross_edges_to_p;
			displs[p] = total_number_of_messages;
			total_number_of_messages += cross_edges_to_p;
		}

		std::vector<Message> sendbuf(total_number_of_messages,
									 Message(g.dummy_edge.n1,g.dummy_edge.n2));
		std::vector<Message> recvbuf(total_number_of_messages);

		//add messages to sendbuf
		for(int p=0; p<procs; p++)
		{
			for(unsigned int i=0;
				i<messages_of_proc[p].first.size();
				i++)
			{
				sendbuf[i+displs[p]] = messages_of_proc[p].first[i];
			}
		}

		MPI_Alltoallv(&sendbuf[0], &counts[0], &displs[0], CANDIDATE_MESSAGE_TYPE,
					  &recvbuf[0], &counts[0], &displs[0], CANDIDATE_MESSAGE_TYPE,
					  MPI_COMM_WORLD);


		// set received candidates
		for(int p=0; p<procs; p++)
		{
			unsigned int i=displs[p];
			unsigned int upper_bound = i+counts[p];

			while(i<upper_bound && recvbuf[i].src!=g.dummy_edge.n1)
			{
				Message msg = recvbuf[i];
				NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(msg.src); /// \todo not nice, get_local_ghost_id is expensive
				NodeIDType local_id		  = g.global_vertex_id_to_local_id(msg.candidate);

				if(candidate_of_node[local_id].partner == local_ghost_id)
				{
					/// \todo I don't think that we have to adjust candidate_of_node, because
					/// the candidates have been set in an earlier step, otherwise we wouldn't
					/// have received this message or the ghostnode wouldn't be the partner
					matched[local_id] = true;
					matched[local_ghost_id] = true;
				}

				i++;
			}
		}

	}


	static void send_messages(const Graph &g, std::vector< std::pair< std::vector<Message>, MPI_Request> > &messages_of_proc)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		msgs_to_receive = 0;

		int tag = msg_tag_candidates^round;

		//for each active partner: send messages;
		/// \todo think about using a list with active partners
		for(int p=0; p<procs; p++)
		{
			if(!messages_of_proc[p].first.empty())
			{// check for messages_of_proc[p].size() == g.get_active_cross_edges_of_proc(p)??

#ifdef MORE_INFORMATION
				number_of_sent_messages++;
				sent_bytes += messages_of_proc[p].first.size()*sizeof(Message);
#endif
				MPI_Isend(&messages_of_proc[p].first[0],
						   messages_of_proc[p].first.size(),
						   CANDIDATE_MESSAGE_TYPE, p, tag,
						   MPI_COMM_WORLD, &messages_of_proc[p].second);

				msgs_to_receive++;
			}
			else if(g.is_active_partner(p))
			{
#ifdef MORE_INFORMATION
				number_of_sent_messages++;
#endif
				// we have to send empty message to active partner, because they expect a message from each active partner
				// even if this message is empty!!
				MPI_Isend(0, 0, CANDIDATE_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD, &messages_of_proc[p].second);

				msgs_to_receive++;
			}
		}
	}


	static void send_messages(const Graph &g, std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > &messages_of_proc)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		msgs_to_receive = 0;

		int tag = msg_tag_ghost_vertices^round;

		//for each active partner: send messages;
		/// \todo think about using a list with active partners
		for(int p=0; p<procs; p++)
		{
			if(!messages_of_proc[p].first.empty())
			{// check for messages_of_proc[p].size() == g.get_active_cross_edges_of_proc(p)??
#ifdef MORE_INFORMATION
				number_of_sent_messages++;
				sent_bytes += messages_of_proc[p].first.size()*sizeof(NodeIDType);
#endif
				MPI_Isend(&messages_of_proc[p].first[0], messages_of_proc[p].first.size(),
						MATCHED_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD,
						&messages_of_proc[p].second);

				msgs_to_receive++;
			}
			else if(g.is_active_partner(p))
			{
#ifdef MORE_INFORMATION
				number_of_sent_messages++;
#endif
				// we have to send empty message to active partner, because they expect a message from each active partner
				// even if this message is empty!!
				MPI_Isend(0, 0, MATCHED_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD, &messages_of_proc[p].second);

				msgs_to_receive++;
			}
		}
	}

	static void exchange_ghost_vertices_all_to_all(const Graph &g,
									  	  	  	     const std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > &messages_of_proc,
									  	  	  	     std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		std::vector<int> counts(procs);
		std::vector<int> displs(procs);

		int total_number_of_messages = 0;

		for(int p=0; p<procs; p++)
		{
			int cross_edges_to_p = g.get_active_cross_edges_count_of_proc(p);
#ifdef MORE_INFORMATION
			cross_edges_to_p *= communication_factor;
#endif
			counts[p] = cross_edges_to_p;
			displs[p] = total_number_of_messages;
			total_number_of_messages += cross_edges_to_p;
		}

		std::vector<NodeIDType> sendbuf(total_number_of_messages, g.dummy_vertex);
		std::vector<NodeIDType> recvbuf(total_number_of_messages);

		//add messages to sendbuf
		for(int p=0; p<procs; p++)
		{
			for(unsigned int i=0;
				i<messages_of_proc[p].first.size();
				i++)
			{
				sendbuf[i+displs[p]] = messages_of_proc[p].first[i];
			}
		}

		MPI_Alltoallv(&sendbuf[0], &counts[0], &displs[0], MATCHED_MESSAGE_TYPE,
					  &recvbuf[0], &counts[0], &displs[0], MATCHED_MESSAGE_TYPE,
					  MPI_COMM_WORLD);


		// set received candidates
		for(int p=0; p<procs; p++)
		{
			unsigned int i=displs[p];
			unsigned int upper_bound = i+counts[p];

			while(i<upper_bound && recvbuf[i]!=g.dummy_vertex)
			{
				const NodeIDType global_ghost_id = recvbuf[i];
				const NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(global_ghost_id); /// \todo not nice, get_local_ghost_id is expensive

				matched[local_ghost_id] = true;

				i++;
			}
		}

	}


	static void receive_ghost_candidates(const Graph &g, const std::vector<Candidate> &candidate_of_node, std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Status msg_status;

		int tag = msg_tag_candidates^round;

		// receive messages
//		for(int p=0; p<procs; p++)
//		{
//			if(!g.is_active_partner(p))
//			{
//				continue;
//			}
		while(msgs_to_receive>0)
		{
			MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &msg_status);

			int msg_count;
			MPI_Get_count(&msg_status, CANDIDATE_MESSAGE_TYPE, &msg_count);

			int src = msg_status.MPI_SOURCE;

			Message *messages = new Message[msg_count];

			MPI_Recv(messages, msg_count, CANDIDATE_MESSAGE_TYPE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			msgs_to_receive--;

			for(int i=0; i<msg_count; i++)
			{
				Message msg = messages[i];
				NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(msg.src); /// \todo not nice, get_local_ghost_id is expensive
				NodeIDType local_id		  = g.global_vertex_id_to_local_id(msg.candidate);

				if(candidate_of_node[local_id].partner == local_ghost_id)
				{
					/// \todo I don't think that we have to adjust candidate_of_node, because
					/// the candidates have been set in an earlier step, otherwise we wouldn't
					/// have received this message or the ghostnode wouldn't be the partner
					matched[local_id] = true;
					matched[local_ghost_id] = true;
				}

			}

			delete[] messages;
		}

	}


	static void receive_ghost_candidates_irecv(const Graph &g, const std::vector<Candidate> &candidate_of_node, std::vector<bool> &matched)
	{
		int procs, proc_id;
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		/// \todo would be nice to set the size of recv_requests to number of active partners

//		std::vector<MPI_Request> recv_requests(g.get_active_partner_count(), MPI_REQUEST_NULL);

		int recv_count = 0;

		for(int p=0; p<procs; p++)
		{
			if(!g.is_active_partner(p))
			{
				continue;
			}

			MPI_Irecv(&candidate_recv_messages[p][0], candidate_recv_messages[p].size(), CANDIDATE_MESSAGE_TYPE,
					   p, msg_tag_candidates, MPI_COMM_WORLD, &recv_requests[/*recv_count*/p]);
			recv_count++;
		}


		int recvd_msgs;

//		std::vector<int> positions(recv_messages.size());
		std::vector<MPI_Status> statuses(candidate_recv_messages.size(), MPI_Status());

		while(recv_count>0)
		{
//			MPI_Waitsome(recv_requests.size(), &recv_requests[0], &recvd_msgs, &positions[0], &statuses[0]);

			MPI_Waitall(recv_requests.size(), &recv_requests[0], &statuses[0]);
			recvd_msgs = recv_requests.size();


			recv_count -= recvd_msgs;

			for(int m=0; m<recvd_msgs; m++)
			{
				int src_proc = statuses[m].MPI_SOURCE;
				int msgs_of_proc;

				MPI_Get_count(&statuses[m], CANDIDATE_MESSAGE_TYPE, &msgs_of_proc);

//				 adjust msgs_of_proc, because we counted bytes
//				msgs_of_proc /= sizeof(Message);


				for(int i=0; i<msgs_of_proc; i++)
				{
					Message msg = candidate_recv_messages[src_proc][i];
					NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(msg.src); /// \todo not nice, get_local_ghost_id is expensive
					NodeIDType local_id		  = g.global_vertex_id_to_local_id(msg.candidate);

					if(candidate_of_node[local_id].partner == local_ghost_id)
					{
						/// \todo I don't think that we have to adjust candidate_of_node,
						/// because the candidates have been set in an earlier step, otherwise
						/// we wouldn't have received this message or the ghostnode wouldn't be the partner
						matched[local_id] = true;
						matched[local_ghost_id] = true;
					}

				}
			}

		}


	}


	static void set_recv_buffers(const Graph &g)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		candidate_recv_messages.resize(procs);
		matched_recv_messages.resize(procs);
		recv_requests.resize(procs, MPI_REQUEST_NULL);

		for(int p=0; p<procs; p++)
		{
			candidate_recv_messages[p].resize(g.get_active_cross_edges_count_of_proc(p));
			matched_recv_messages[p].resize(g.get_active_cross_edges_count_of_proc(p));
		}

	}

	static void clear_recv_buffers()
	{
		candidate_recv_messages.clear();
		matched_recv_messages.clear();
		recv_requests.clear();
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


	static void deactivate_local_edges_incident_to_matched_nodes(Graph &g, const std::vector<bool> &matched, std::vector<Candidate > &candidate_of_node)
	{
		Candidate base_candidate = Candidate();
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
	}


	static void add_matched_cross_edges_to_matching_and_prepare_messages_for_ghost_nodes_adjacent_to_matched_local_nodes
												(
													const Graph &g,
													const std::vector<Candidate> &candidate_of_node,
													const std::vector<bool> &matched,
													std::vector<Edge> &matching,
													std::vector< std::pair< std::vector<NodeIDType>, MPI_Request> > &messages_of_proc
												)
	{
		for(const_edge_iterator e_it=g.begin_active_local_cross_edges(),
						  it_end=g.end_active_local_cross_edges();
			e_it!=it_end; e_it++)
		{
			// check first endpoint of current edge
			if(matched[e_it->n1] || matched[e_it->n2])
			{
				/// \todo because of this condition no multi-graphs are supported which have multi-edges with the same weight
				if(matched[e_it->n1] && matched[e_it->n2] && candidate_of_node[e_it->n1].partner==e_it->n2 && candidate_of_node[e_it->n2].partner==e_it->n1)
				{/// \todo only add matched cross edges on one process, maybe on the one with the smaller proc-ID (or use the node ID to decide)

					if( (g.is_local(e_it->n1) && g.local_vertex_id_to_global_id(e_it->n1)<g.local_vertex_id_to_global_id(e_it->n2)) ||
						(g.is_local(e_it->n2) && g.local_vertex_id_to_global_id(e_it->n2)<g.local_vertex_id_to_global_id(e_it->n1))
					  )
					{ // only add the edge if the global id of the local endpoint is smaller than the global id of the ghost endpoint
						matching.push_back(Edge(g.local_vertex_id_to_global_id(e_it->n1),
												g.local_vertex_id_to_global_id(e_it->n2),
												e_it->weight));
					}
				}
				else if(matched[e_it->n1] && g.is_local(e_it->n1))
				{
#ifdef MORE_INFORMATION
					for(unsigned int i=0; i<communication_factor; i++)
#endif
					messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n2)].first.push_back(g.local_vertex_id_to_global_id(e_it->n1));
				}
				else if(matched[e_it->n2] && g.is_local(e_it->n2))
				{
#ifdef MORE_INFORMATION
					for(unsigned int i=0; i<communication_factor; i++)
#endif
					messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n1)].first.push_back(g.local_vertex_id_to_global_id(e_it->n2));
				}

			}

		}
	}


	static void receive_information_about_matched_ghost_nodes(const Graph &g, std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		MPI_Status msg_status;

		int tag = msg_tag_ghost_vertices^round;

		// receive matched messages
//		for(int p=0; p<procs; p++)
//		{
//			if(!g.is_active_partner(p))
//			{// if we send a message to p, we also expect a response from this process, even if we no longer have any cross-edges this process
//				continue;
//			}
		while(msgs_to_receive>0)
		{
			/// \todo maybe use ANY_SOURCE instead (don't wait for a specific message). In this case we need a loop the iterates over all "incoming messages"
			MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &msg_status);

			int msg_count;
			MPI_Get_count(&msg_status, MATCHED_MESSAGE_TYPE, &msg_count);//g.get_active_cross_edges_count_of_proc(p);

//			msg_count /= sizeof(Message); // adjust msg count, because we only counted the number of Bytes, well maybe not the nicest thing to do

			NodeIDType *messages = new NodeIDType[msg_count];

			int src = msg_status.MPI_SOURCE;

			MPI_Recv(messages, msg_count, MATCHED_MESSAGE_TYPE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			msgs_to_receive--;

			for(int i=0; i<msg_count; i++)
			{
				const NodeIDType global_ghost_id = messages[i];
				const NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(global_ghost_id); /// \todo not nice, get_local_ghost_id is expensive

				matched[local_ghost_id] = true;
			}

			delete[] messages;
		}
	}


	static void receive_information_about_matched_ghost_nodes_irecv(const Graph &g, std::vector<bool> &matched)
	{
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);

		/// \todo would be nice to set the size of recv_requests to number of active partners

//		std::vector<MPI_Request> recv_requests(g.get_active_partner_count(), MPI_REQUEST_NULL);

		int recv_count = 0;

		for(int p=0; p<procs; p++)
		{
			if(g.is_active_partner(p))
			{
				MPI_Irecv(&matched_recv_messages[p][0], matched_recv_messages[p].size(), MATCHED_MESSAGE_TYPE, p, msg_tag_ghost_vertices, MPI_COMM_WORLD, &recv_requests[/*recv_count*/p]);
				recv_count++;
			}
		}


		int recvd_msgs;

//		std::vector<int> positions(matched_recv_messages.size());
		std::vector<MPI_Status> statuses(matched_recv_messages.size(), MPI_Status());

		while(recv_count>0)
		{
//			MPI_Waitsome(recv_requests.size(), &recv_requests[0], &recvd_msgs, &positions[0], &statuses[0]);

			MPI_Waitall(recv_requests.size(), &recv_requests[0], &statuses[0]);
			recvd_msgs = recv_requests.size();

			recv_count -= recvd_msgs;

			for(int m=0; m<recvd_msgs; m++)
			{
				int src_proc = statuses[m].MPI_SOURCE;
				int msgs_of_proc;

				MPI_Get_count(&statuses[m], MATCHED_MESSAGE_TYPE, &msgs_of_proc);

//				// adjust msgs_of_proc, because we counted bytes
//				msgs_of_proc /= sizeof(Message);


				for(int i=0; i<msgs_of_proc; i++)
				{
					const NodeIDType global_ghost_id = matched_recv_messages[src_proc][i];

					const NodeIDType local_ghost_id = g.global_vertex_id_to_local_id_of_ghost_vertex(global_ghost_id); /// \todo not nice, get_local_ghost_id is expensive

					matched[local_ghost_id] = true;
				}
			}
		}

	}


	static void deactivate_cross_edges_incident_to_matched_nodes(Graph &g, const std::vector<bool> &matched, std::vector<Candidate > &candidate_of_node)
	{
		Candidate base_candidate = Candidate();

		// deactivate all remaining cross edges that are incident to a matched ghost-vertex
		edge_iterator e_it=g.begin_active_local_cross_edges();
		while(e_it!=g.end_active_local_cross_edges())
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

	static inline void deactivate_edges_incident_to_matched_nodes(Graph &g, const std::vector<bool> &matched, std::vector<Candidate > &candidate_of_node)
	{
		deactivate_local_edges_incident_to_matched_nodes(g, matched, candidate_of_node);
		deactivate_cross_edges_incident_to_matched_nodes(g, matched, candidate_of_node);
	}


	static void init()
	{
		MPI_Type_contiguous(sizeof(Message), MPI_BYTE, &CANDIDATE_MESSAGE_TYPE);
		MPI_Type_commit(&CANDIDATE_MESSAGE_TYPE);

		MPI_Type_contiguous(sizeof(NodeIDType), MPI_BYTE, &MATCHED_MESSAGE_TYPE);
		MPI_Type_commit(&MATCHED_MESSAGE_TYPE);
	}

};


// Define within in file scope -> static members aren't automatically defined, you have to "use" them.
// http://msdn.microsoft.com/en-us/library/b1b5y48f.aspx
template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector<MPI_Request> ParallelLocalMaximumMatching<Graph, adjust_id>::recv_requests;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::vector<typename ParallelLocalMaximumMatching<Graph, adjust_id>::Message> >
	ParallelLocalMaximumMatching<Graph, adjust_id>::candidate_recv_messages;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::vector<typename ParallelLocalMaximumMatching<Graph, adjust_id>::NodeIDType> >
	ParallelLocalMaximumMatching<Graph, adjust_id>::matched_recv_messages;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
MPI_Datatype ParallelLocalMaximumMatching<Graph, adjust_id>::CANDIDATE_MESSAGE_TYPE;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
MPI_Datatype ParallelLocalMaximumMatching<Graph, adjust_id>::MATCHED_MESSAGE_TYPE;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::pair< std::vector<typename ParallelLocalMaximumMatching<Graph, adjust_id>::Message>, MPI_Request> >
	ParallelLocalMaximumMatching<Graph, adjust_id>::candidate_messages_of_proc;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
std::vector< std::pair< std::vector<typename ParallelLocalMaximumMatching<Graph, adjust_id>::NodeIDType>, MPI_Request> >
	ParallelLocalMaximumMatching<Graph, adjust_id>::matched_messages_of_proc;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
unsigned int ParallelLocalMaximumMatching<Graph, adjust_id>::msgs_to_receive = 0;


template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
unsigned int ParallelLocalMaximumMatching<Graph, adjust_id>::round = 0;


#ifdef MORE_INFORMATION
template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
double ParallelLocalMaximumMatching<Graph, adjust_id>::communication_time = 0.;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
long long ParallelLocalMaximumMatching<Graph, adjust_id>::number_of_sent_messages = 0;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
long long ParallelLocalMaximumMatching<Graph, adjust_id>::sent_bytes = 0;

template <typename Graph, typename Graph::NodeIDType (*adjust_id)(const typename Graph::NodeIDType)>
unsigned int ParallelLocalMaximumMatching<Graph, adjust_id>::communication_factor = 1;
#endif


}


#endif


