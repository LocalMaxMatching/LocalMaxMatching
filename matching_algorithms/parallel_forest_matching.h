#ifndef _TREE_MATCHING_H_
#define _TREE_MATCHING_H_

#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <math.h>

#include <fstream>

#include <mpi.h>

#include <graphs/parallel_forest.h>
#include <common/union_find.h>
#include <common/hash_functions.h>
#include <common/math_funcs.h>


namespace dipl
{


template <class Graph>
class ParallelForestMatching
{
	typedef typename Graph::Edge Edge;
	typedef typename Graph::NodeIDType	NodeIDType;
	typedef typename Graph::NodeIDType	EdgeIDType;
	typedef typename Graph::WeightType	WeightType;
	typedef AdjacencyList<NodeIDType, EdgeIDType> Forest;


	struct Partner
	{
		int id; // Pocess ID
		NodeIDType local_node, ghost_node;
		bool is_active_src; // states if we still receive messages from this process
		bool is_active_tgt; // states if we still send messages to this process
		unsigned int msgs_from_other_partners; // used for: do we receive messages from c\p? if > 0: we receive messages from other partners from the current component

		Partner()
		: id(0), local_node(0), ghost_node(0), is_active_src(true), is_active_tgt(true)
		{}


		Partner(const int id, const NodeIDType ln, const NodeIDType gn)
		: id(id), local_node(ln), ghost_node(gn), is_active_src(true), is_active_tgt(true)
		{}

		friend std::ostream& operator<< (std::ostream &out, const Partner &p)
		{
			out << "(" << p.id << ", (" << p.local_node << ", "
					<< p.ghost_node << "), " << p.is_active_src << ", "
					<< p.is_active_tgt << ", " << p.msgs_from_other_partners<< ")";
			return out;
		}

	};



	struct ComponentMSG
	{
		NodeIDType sender_local_node;
		NodeIDType sender_ghost_node;

		/// positive -> candidate value
		/// negative -> no more messages
		long msg_info;

		friend std::ostream& operator<< (std::ostream &out, const ComponentMSG &msg)
		{
			out << "sender_local_node: " << msg.sender_local_node << ", sender_ghost_node: " << msg.sender_ghost_node << ", msg_info: " << msg.msg_info;
			return out;
		}
	};

	MPI_Datatype ComponentMSGType;


	struct SubtreeResultMSG
	{
		NodeIDType subtree_root; // is ghost node on receiving side - local on sending side
		NodeIDType hint_node; // is local node on receiving side that is adjacent to subtree_root

		WeightType incoming_matched;
		WeightType incoming_not_matched;

		bool uses_outgoing_edge;


		friend std::ostream& operator<< (std::ostream &out, const SubtreeResultMSG &msg)
		{
			out << "subtree_root=" << msg.subtree_root << ", hin_node=" << msg.hint_node
					<< ", weight incoming matched=" << msg.incoming_matched
					<< ", weight incoming not matched=" << msg.incoming_not_matched
					<< ", uses_outgoing_edge=" << msg.uses_outgoing_edge;

			return out;
		}
	};

	MPI_Datatype SubtreeMSGType;

	struct TopDownMSG
	{
		NodeIDType tgt_node; // is local on receiving side and ghost on sending side
		unsigned char msg_info;

		static const unsigned char use_incoming_edge = 0;
		static const unsigned char parent_matched = 1;
		static const unsigned char parent_not_matched = 2;
	};

	MPI_Datatype TopDownMSGType;


	struct ConnectionComponent
	{
		typedef typename std::list<Partner>::iterator iterator;
		typedef typename std::list<Partner>::const_iterator const_iterator;

		std::list< Partner > partner_set;

		NodeIDType root_candidate; // id of the root node

		NodeIDType root_edge_local_node;
		NodeIDType root_edge_ghost_node;
//		Partner *local_root;

		ConnectionComponent()
		: root_candidate(std::numeric_limits<NodeIDType>::max()),
		  root_edge_local_node(-1),
		  root_edge_ghost_node(-1)
//		  local_root(NULL)
		{}


		iterator begin()
		{
			return partner_set.begin();
		}

		const_iterator begin() const
		{
			return partner_set.begin();
		}

		iterator end()
		{
			return partner_set.end();
		}

		const_iterator end() const
		{
			return partner_set.end();
		}

		void push_back(const Partner &e)
		{
			partner_set.push_back(e);
		}

		Partner& back()
		{
			return partner_set.back();
		}

		const Partner& back() const
		{
			return partner_set.back();
		}


		unsigned int size() const
		{
			return partner_set.size();
		}


		bool empty() const
		{
			return partner_set.empty();
		}


		void clear()
		{
			partner_set.clear();
			root_candidate = std::numeric_limits<long>::max();
			root_edge_local_node = -1;
			root_edge_ghost_node = -1;
//			local_root = NULL;
		}


		template <typename T>
		friend std::ostream& operator<< (std::ostream &out, const ConnectionComponent &c);
	};

	struct BorderComponent
	{
		typedef typename std::list<Partner>::iterator iterator;
		typedef typename std::list<Partner>::const_iterator const_iterator;

//		std::list< Partner > partner_set;

		Edge root_edge;
		bool is_global_max;

		BorderComponent()
		: root_edge(0,0,dipl::min_val<double>()), is_global_max(false)
		{}

//
//		iterator begin()
//		{
//			return partner_set.begin();
//		}
//
//		const_iterator begin() const
//		{
//			return partner_set.begin();
//		}
//
//		iterator end()
//		{
//			return partner_set.end();
//		}
//
//		const_iterator end() const
//		{
//			return partner_set.end();
//		}
//
//		void push_back(const Partner &e)
//		{
//			partner_set.push_back(e);
//		}
//
//		Partner& back()
//		{
//			return partner_set.back();
//		}
//
//		const Partner& back() const
//		{
//			return partner_set.back();
//		}
//
//
//		unsigned int size() const
//		{
//			return partner_set.size();
//		}
//
//
//		bool empty() const
//		{
//			return partner_set.empty();
//		}


		void clear()
		{
//			partner_set.clear();
			root_edge = Edge(0,0,dipl::min_val<double>());
		}


//		template <typename T>
//		friend std::ostream& operator<< (std::ostream &out, const BorderComponent &c);
	};

	struct BorderComponentMSG
	{
		NodeIDType sender_local;
		NodeIDType sender_ghost;
	};

	MPI_Datatype BorderComponentMSGType;

	typedef std::list<ConnectionComponent>	ConnectionComponents;

	typedef std::list<BorderComponent>	BorderComponents;

	template <typename T>
	friend std::ostream& operator<< (std::ostream &out, const typename ParallelForestMatching<T>::ConnectionComponents &c);

private:
	const Graph &g;

	// UnionFind-datastructure representing sets of locally connected components
	dipl::UnionFind<NodeIDType> uf; // don't forget to reset this values after each call of get_matching

//	std::vector<ConnectionComponent*> connection_component_ptr_of_vertex;
	std::vector<BorderComponent*> connection_component_ptr_of_vertex;

	NodeIDType num_local_vertices, num_ghost_vertices, log_num_ghost_vertices;

	Forest forest;

	std::vector<NodeIDType> pending_msgs_of_root;

	/*!
	 * Stores for each node the two ratings for the subtree starting at this node.
	 * The first rating-value specifies the rating if the incoming edge to this node is
	 * matched and the second rating value specifies the rating if the incoming edge
	 * isn't matched.
	 * We also store the ID of the child of the outgoing edge that is matched if the
	 * incoming isn't matched (-1 if no outgoing edge is matched???)
	 */
	std::vector<std::pair<std::pair<typename Graph::WeightType, typename Graph::WeightType>, typename Graph::NodeIDType> > subtree_table;


	/// pointer to the vector storing the root nodes of the local trees
	std::vector<NodeIDType> *roots_of_local_trees_ptr;

	/// index of the next local tree to be considered.
	/// matchings for local trees at positions < next_local_tree have already
	/// been computed. Initially has to be 0
	NodeIDType	next_local_tree_pos;

#ifdef MORE_INFORMATION
public:
	struct ParallelTreeInfo
	{
		NodeIDType subtree_root; // is ghost node on receiving side - local on sending side
		NodeIDType hint_node; // is local node on receiving side that is adjacent to subtree_root

		NodeIDType depth;
		NodeIDType number_of_nodes;

		NodeIDType component_depth;
		NodeIDType number_of_components;
	};

	MPI_Datatype ParallelTreeInfoMSGType;

	std::vector<ParallelTreeInfo> subtree_info;
	std::vector< std::list< std::vector<ParallelTreeInfo> > > proc_subtree_info_buffers;

	std::ofstream outfile_parallel_tree_info;
	std::ofstream outfile_std_out;

	static const int parallel_tree_info_tag = 11;
private:
#endif


	int procs, proc_id;

	std::vector< std::list<MPI_Request> > proc_requests;
	std::vector< std::list< std::vector<SubtreeResultMSG> > > proc_bottom_up_buffers;
	std::vector< std::list< std::vector<TopDownMSG> > > proc_top_down_buffers;

	static const int bottom_up_tag = 8; // need a better mechanism for defining unique msg tags over the whole project.
	static const int top_down_tag = 9;

	static const int decide_on_root_node_tag = 10;

	std::vector<bool> proc_is_active_src;

public:
	long number_of_mpi_messages_sent, number_of_mpi_messages_received;
	long calls_to_g_for_gid_to_lid;

private:

	template <typename T>
	void wait_for_message_completion_and_reset_buffers(std::vector< std::list<MPI_Request> > &proc_requests,
													   std::vector< std::list< std::vector<T> > > &buffers)
	{
		for(unsigned int p=0; p<proc_requests.size(); p++)
		{
			for(std::list<MPI_Request>::iterator req_it=proc_requests[p].begin(),
												 req_it_end=proc_requests[p].end();
					req_it!=req_it_end; req_it++)
			{
				MPI_Wait(&(*req_it), MPI_STATUS_IGNORE);
			}

			proc_requests[p].clear();
			buffers[p].clear();
		}
	}



	/// \todo unfortunately this functions isn't linear time, because we use a union find
	/// structure - but it's almost linear. the uf.find(...) operation isn't constant.
	///
	/// we could use a similar technique as in build_forest to find the
	/// connected components in linear time
	void get_connection_components(const std::vector<EdgeIDType> &local_edges,
		  	  	  	  	  	  	   const std::vector<EdgeIDType> &cross_edges,
								   ConnectionComponents &connection_components,
								   std::vector<NodeIDType> &roots_of_local_trees)
	{

		for(typename std::vector<EdgeIDType>::const_iterator eit=local_edges.begin(),
															 eit_end=local_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);
			uf.combine(e.n1, e.n2);
		}


		for(typename std::vector<EdgeIDType>::const_iterator eit=cross_edges.begin(),
															 eit_end=cross_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			NodeIDType rep = uf.find(e.n1);

			ConnectionComponent *c_ptr = connection_component_ptr_of_vertex[rep];

			if(c_ptr==NULL)
			{// haven't added a component yet for this node
				connection_components.push_back(ConnectionComponent());
				c_ptr = &connection_components.back();
				connection_component_ptr_of_vertex[rep] = c_ptr;
				uf.set_visited(rep, true); //! < \todo do i really need this? i think we need it at least for the roots of local trees, otherwise we might add a root twice
			}

			// let the local node point to it's connection component
			// the connection component is unique for each local node
			connection_component_ptr_of_vertex[e.n1] = c_ptr;

			// n2 is always the ghost node
			c_ptr->push_back(Partner(g.get_proc_of_ghost_vertex(e.n2),
									 e.n1,
									 e.n2));

		}


		// add the local trees - we've already added the parallel trees, thus all the remaining
		// (not visited) components must be local trees
		for(typename std::vector<EdgeIDType>::const_iterator eit=local_edges.begin(),
															 eit_end=local_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			NodeIDType rep = uf.find(e.n1); // == uf.find(eit->n2) because they're both in the same connected component
			if( ! uf.get_visited(rep) )
			{
				roots_of_local_trees.push_back(rep);
				uf.set_visited(rep, true);
			}
		}

	}


	void get_border_components(const std::vector<EdgeIDType> &local_edges,
		  	  	  	  	  	  	   const std::vector<EdgeIDType> &cross_edges,
								   BorderComponents &border_components,
								   std::vector<NodeIDType> &roots_of_local_trees)
	{

		for(typename std::vector<EdgeIDType>::const_iterator eit=local_edges.begin(),
															 eit_end=local_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);
			uf.combine(e.n1, e.n2);
		}


		for(typename std::vector<EdgeIDType>::const_iterator eit=cross_edges.begin(),
															 eit_end=cross_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			NodeIDType rep = uf.find(e.n1);

			BorderComponent *b_ptr = connection_component_ptr_of_vertex[rep];

			if(b_ptr==NULL)
			{// haven't added a component yet for this node
				border_components.push_back(BorderComponent());
				b_ptr = &border_components.back();
				connection_component_ptr_of_vertex[rep] = b_ptr;
				uf.set_visited(rep, true); //! < \todo do i really need this? i think we need it at least for the roots of local trees, otherwise we might add a root twice
			}

			// let the local node point to it's connection component
			// the connection component is unique for each local node
			connection_component_ptr_of_vertex[e.n1] = b_ptr;

			// update heaviest cross edge, the heaviest cross edge is the root of each
			// local tree, and the heaviest cross edge of a parallel tree is its root
			if(dipl::vertex_tie_breaking(g,e, b_ptr->root_edge))
			{
				b_ptr->root_edge = e;
			}

//			// n2 is always the ghost node
//			b_ptr->push_back(Partner(g.get_proc_of_ghost_vertex(e.n2),
//									 e.n1,
//									 e.n2));

		}


		// add the local trees - we've already added the parallel trees, thus all the remaining
		// (not visited) components must be local trees
		for(typename std::vector<EdgeIDType>::const_iterator eit=local_edges.begin(),
															 eit_end=local_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			NodeIDType rep = uf.find(e.n1); // == uf.find(eit->n2) because they're both in the same connected component
			if( ! uf.get_visited(rep) )
			{
				roots_of_local_trees.push_back(rep);
				uf.set_visited(rep, true);
			}
		}

	}

	void decide_on_global_roots(const Graph &g,
								BorderComponents &bcs)
	{
		std::vector<std::vector<BorderComponentMSG> > msg_for_proc(procs);
		std::vector<MPI_Request> req_of_msg_for_proc(procs, MPI_REQUEST_NULL);

		EdgeIDType number_of_recv_msgs = 0;

		for(typename BorderComponents::iterator b_it=bcs.begin(),
												b_it_end=bcs.end();
			b_it!=b_it_end; b_it++)
		{
				BorderComponentMSG msg;
				msg.sender_local = g.local_vertex_id_to_global_id(b_it->root_edge.n1);
				msg.sender_ghost = g.local_vertex_id_to_global_id(b_it->root_edge.n2);

				int tgt_proc = g.get_proc_of_ghost_vertex(b_it->root_edge.n2);

				msg_for_proc[tgt_proc].push_back(msg);


		}

		for(int p=0; p<procs; p++)
		{
			if(g.is_active_partner(p))
			{// we allow empty messages!! but only to active partners
				MPI_Issend(&msg_for_proc[p][0], msg_for_proc[p].size(),
						   BorderComponentMSGType, p, decide_on_root_node_tag,
						   MPI_COMM_WORLD, &req_of_msg_for_proc[p]);
				number_of_recv_msgs++;
			}
		}

		unsigned int m;
		// receive all msgs
		for(m=0; m<number_of_recv_msgs; )
		{
			MPI_Status status;

			MPI_Probe(MPI_ANY_SOURCE, decide_on_root_node_tag, MPI_COMM_WORLD, &status);

			int src_proc = status.MPI_SOURCE;
			int count;
			MPI_Get_count(&status, BorderComponentMSGType, &count);

			// adjust the number of received messages
			m++; // one message from each active partner

			std::vector<BorderComponentMSG> msgs(count);
			MPI_Recv(&msgs[0], count, ComponentMSGType, src_proc, decide_on_root_node_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(int i=0; i<count; i++)
			{
				BorderComponentMSG &msg = msgs[i];

				NodeIDType local_node = g.global_vertex_id_to_local_id_of_local_vertex(msg.sender_ghost);
				NodeIDType ghost_node = g.global_vertex_id_to_local_id_of_ghost_vertex(msg.sender_local);

				// only set the corresponding bordercomponent for the local vertices,
				// it's not necessarily unique for ghost vertices
				BorderComponent &b = *(connection_component_ptr_of_vertex[local_node]);

				if(b.root_edge.n1==local_node && b.root_edge.n2==ghost_node)
				{// that's the heaviest edge of the corresponding parallel tree,
				 // the larger vertex is the global root
					if(g.local_vertex_id_to_global_id(local_node)>g.local_vertex_id_to_global_id(ghost_node))
					{// this border component is the global root of the corresponding parallel tree
						b.is_global_max = true;
					}
				}
			}

		}

		MPI_Waitall(procs, &req_of_msg_for_proc[0], MPI_STATUSES_IGNORE);
	}


	/*!
	 * \brief Like MPI_Probe, but it probes for any src that is active.
	 *
	 * If there's a waiting message from an active src, the function sets this
	 * src to inactive and returns true. Information about the incoming message
	 * is stored in status.
	 *
	 * This function returns false if all sources in active_src are inactive
	 * (i.e. they're set to false):
	 *
	 * @param active_src	Bool vector that specifies for each process in comm
	 * 						if it is active (i.e. set to true)
	 * @param tag			The tag for incoming messages
	 * @param comm			Communicator
	 * @param status		Status object, returns information about the incoming
	 * 						message.
	 *
	 * @return				Returns true if there's an incoming message and false
	 * 						if all sources are inactive
	 */
	bool probe_some(std::vector<bool> &active_src, int tag, MPI_Comm comm, MPI_Status *status)
	{
		int flag;
		bool active_sources = true;
		while(active_sources)
		{
			active_sources = false;

			for(int p=0; p<(int)active_src.size(); p++)
			{
				if(active_src[p])
				{
					MPI_Iprobe(p, tag, comm, &flag, status);

					if(flag)
					{
						active_src[p] = false;
						return true;
					}

					active_sources = true; //set active sources to true if there's at least one active src
				}
			}
		}

		return false;
	}


	void decide_on_root_process(const Graph &g,
								ConnectionComponents &cc)
	{
		/// \todo remove unnecessary Partners from components
		/// those are the partner we don't send any messages to and from which we also don't receive any messages

		/// \todo maybe too big - only send to active partners
		std::vector<std::vector<ComponentMSG> > msg_for_proc(procs);
		std::vector<MPI_Request> req_of_msg_for_proc(procs, MPI_REQUEST_NULL);

		unsigned int number_of_send_msgs = 0;
		unsigned int number_of_recv_msgs = 0;

		for(typename ConnectionComponents::iterator c_it=cc.begin(),
													c_it_end=cc.end();
			c_it!=c_it_end; c_it++)
		{
			for(typename ConnectionComponent::iterator p_it=c_it->begin(),
													   p_it_end=c_it->end();
				p_it!=p_it_end; p_it++)
			{
				Partner &p = *p_it;

				// count the number of msgs we'll receive by all partners
				number_of_recv_msgs += p.is_active_src;

				proc_is_active_src[p.id] = proc_is_active_src[p.id] || p.is_active_src; // specifies if this proc is an active source


				// send messages to the partner
				if(p.msgs_from_other_partners==0)
				{
					if(p.is_active_tgt)
					{
						p.is_active_tgt = false;

						ComponentMSG msg;

						// send msg: no more messages
						msg.msg_info = -1;
						// always send global vertex ids
						msg.sender_ghost_node = g.local_vertex_id_to_global_id(p.ghost_node);
						msg.sender_local_node = g.local_vertex_id_to_global_id(p.local_node);

						msg_for_proc[p.id].push_back(msg);

						number_of_send_msgs++;
					}
				}
				else
				{
					ComponentMSG msg;

					msg.msg_info = c_it->root_candidate;
					// always send global vertex ids
					msg.sender_ghost_node = g.local_vertex_id_to_global_id(p.ghost_node);
					msg.sender_local_node = g.local_vertex_id_to_global_id(p.local_node);

					msg_for_proc[p.id].push_back(msg);

					number_of_send_msgs++;
				}
			}
		}

		for(int p=0; p<procs; p++)
		{
			if(!msg_for_proc[p].empty())
			{
				/// \todo think about using MPI_Isend -> might have a smaller overhead for small messages.
				///	we might need a different tag (counter?)
				MPI_Issend(&msg_for_proc[p][0], msg_for_proc[p].size(), ComponentMSGType, p, decide_on_root_node_tag, MPI_COMM_WORLD, &req_of_msg_for_proc[p]);
				number_of_mpi_messages_sent++;
			}
		}

		unsigned int m;
		// receive all msgs
		for(m=0; m<number_of_recv_msgs; )
		{
			MPI_Status status;

//			MPI_Probe(MPI_ANY_SOURCE, decide_on_root_node_tag, MPI_COMM_WORLD, &status);
			probe_some(proc_is_active_src, decide_on_root_node_tag, MPI_COMM_WORLD, &status);

			int src_proc = status.MPI_SOURCE;
			int count;
			MPI_Get_count(&status, ComponentMSGType, &count);

			// adjust the number of received messages
			m += count;

			std::vector<ComponentMSG> msgs(count);
			MPI_Recv(&msgs[0], count, ComponentMSGType, src_proc, decide_on_root_node_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			number_of_mpi_messages_received++;

			for(int i=0; i<count; i++)
			{
				ComponentMSG &msg = msgs[i];

				NodeIDType local_node = g.global_vertex_id_to_local_id_of_local_vertex(msg.sender_ghost_node);
				NodeIDType local_ghost_node = g.global_vertex_id_to_local_id_of_ghost_vertex(msg.sender_local_node);

				ConnectionComponent *c_ptr = connection_component_ptr_of_vertex[local_node];

				if(msg.msg_info<0)
				{// we won't receive any more msgs from this source

					// adjust the msgs_from_other_partners of the other partners of the current component
					/// \todo think of something faster/better to check the other partners if they're still receiving
					/// messages
					for(typename ConnectionComponent::iterator p_it=c_ptr->begin(),
															   p_it_end=c_ptr->end();
						p_it!=p_it_end; p_it++)
					{
						Partner &p = *p_it;

						if(p.local_node==local_node && p.ghost_node==local_ghost_node)
						{
							p.is_active_src = false;
						}
						else
						{
							p.msgs_from_other_partners--;
						}
					}
				}
				else
				{
					// msg_info is of Type long -> to call the correct hash-function
					// we need a cast to NodeIDType
					NodeIDType hash_msg_info = hash((NodeIDType)msg.msg_info);
					NodeIDType hash_candidate = hash(c_ptr->root_candidate);
					// check for better root
					if( (hash_msg_info < hash_candidate) ||
					   ((hash_msg_info==hash_candidate) && (msg.msg_info<c_ptr->root_candidate)))
					{
						c_ptr->root_candidate = msg.msg_info;
						c_ptr->root_edge_local_node = local_node;
						c_ptr->root_edge_ghost_node = local_ghost_node;
					}
				}

			}

		}

		/// \todo is it possible to get rid of this synchronization
		MPI_Waitall(procs, &req_of_msg_for_proc[0], MPI_STATUSES_IGNORE);

		if(number_of_recv_msgs==0 && number_of_send_msgs==0)
		{
			return ;
		}
		else
		{
			// wait for completion of all send Messages
//			MPI_Waitall(number_of_send_msgs, &send_req_buffer[0], MPI_STATUSES_IGNORE);
			return decide_on_root_process(g, cc);
		}

	}






	void get_structure_of_parallel_trees(const std::vector<EdgeIDType> &local_edges,
										 const std::vector<EdgeIDType> &cross_edges,
										 ConnectionComponents &connection_components,
										 std::vector<NodeIDType> &roots_of_local_trees)
	{
#ifdef MORE_INFORMATION
		double start, end;

		start = MPI_Wtime();
#endif

		get_connection_components(local_edges, cross_edges, connection_components, roots_of_local_trees);

#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run get_connection_components=" << (end-start) << " sec" << std::endl;
#endif

		// set local candidates
//		unsigned int max_send_msgs = 0;

		// set the root candidates

		// set the initial candidates - just use the first partner of each component
		// to set the initial candidate
		for(typename ConnectionComponents::iterator c_it=connection_components.begin(),
													c_it_end=connection_components.end();
			c_it!=c_it_end; c_it++)
		{/// \todo is it possible that there are empty connection components? i don't think so
			typename ConnectionComponent::iterator p_it=c_it->begin();
			c_it->root_candidate = g.local_vertex_id_to_global_id(p_it->local_node);
			c_it->root_edge_local_node = p_it->local_node;
			c_it->root_edge_ghost_node = p_it->ghost_node;
		}

		// find the correct root candidates
		for(typename ConnectionComponents::iterator c_it=connection_components.begin(),
													c_it_end=connection_components.end();
			c_it!=c_it_end; c_it++)
		{
			for(typename ConnectionComponent::iterator p_it=c_it->begin(),
													   p_it_end=c_it->end();
				p_it!=p_it_end; p_it++)
			{
				p_it->msgs_from_other_partners = c_it->size()-1;

				NodeIDType hash_local_node = hash(g.local_vertex_id_to_global_id(p_it->local_node));
				NodeIDType hash_candidate = hash(c_it->root_candidate);

				// using the hash function to distribute the roots of the parallel trees more evenly
				// over the processes
				if(  (hash_local_node < hash_candidate) ||
					((hash_local_node==hash_candidate) && (g.local_vertex_id_to_global_id(p_it->local_node)<c_it->root_candidate)))
				{
					c_it->root_candidate = g.local_vertex_id_to_global_id(p_it->local_node);
					c_it->root_edge_local_node = p_it->local_node;
					c_it->root_edge_ghost_node = p_it->ghost_node;
				}

				NodeIDType hash_ghost_node = hash(g.local_vertex_id_to_global_id(p_it->ghost_node));
				hash_candidate = hash(c_it->root_candidate);

				if(  (hash_ghost_node < hash_candidate) ||
					((hash_ghost_node==hash_candidate) && (g.local_vertex_id_to_global_id(p_it->ghost_node)<c_it->root_candidate)))
				{
					c_it->root_candidate = g.local_vertex_id_to_global_id(p_it->ghost_node);
					c_it->root_edge_local_node = p_it->local_node;
					c_it->root_edge_ghost_node = p_it->ghost_node;
				}
			}

//			max_send_msgs += c_it->size();
		}

#ifdef MORE_INFORMATION
		start = MPI_Wtime();
#endif

		decide_on_root_process(g, connection_components);

#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run decide_on_root_process=" << (end-start) << std::endl;//" sec,   probe time: " << time << ", loop_time: " << loop_time << ", depth: " << depth << ", msgs_sent: " << msgs_sent << ",  cc.size: " << connection_components.size() << std::endl;
#endif

	}






	///! \todo move this function to parallel_forest.h
	void build_forest(const Graph &g,
					  const std::vector<EdgeIDType> &local_edges,
					  const std::vector<EdgeIDType> &cross_edges,
					  const std::vector<NodeIDType> &roots_of_local_trees,
					        std::vector<NodeIDType> &roots_of_parallel_trees,
					  const BorderComponents &border_components)
	{
	//	std::vector<bool> visited(tree.num_vertices(), false);

	/// \todo build the forest using a union-find like structure?? - could be used to identify connected components (trees).
	/// but I assume not directly to create the trees. because when combining two subtrees the nodes that are used to link
	/// those subtrees are most likely not the representatives (i.e. the "directions" of the trees is wrong)

		for(typename std::vector<EdgeIDType>::const_iterator e_it=local_edges.begin(),
															 e_it_end=local_edges.end();
			e_it!=e_it_end; e_it++)
		{
			const Edge &e = g.get_edge(*e_it);

			forest.add_incident_edge(e.n1, e.n2, *e_it);
			forest.add_incident_edge(e.n2, e.n1, *e_it);

			forest.set_visited(e.n1, false);
			forest.set_visited(e.n2, false);
		}

		for(typename std::vector<EdgeIDType>::const_iterator e_it=cross_edges.begin(),
															 e_it_end=cross_edges.end();
			e_it!=e_it_end; e_it++)
		{
			const Edge &e = g.get_edge(*e_it);

			forest.add_incident_edge(e.n1, e.n2, *e_it);
			forest.add_incident_edge(e.n2, e.n1, *e_it);

			forest.set_visited(e.n1, false);
			forest.set_visited(e.n2, false);
		}


		// build local trees
		for(typename std::vector<NodeIDType>::const_iterator r_it=roots_of_local_trees.begin(),
															 r_it_end=roots_of_local_trees.end();
			r_it!=r_it_end; r_it++)
		{
			// the root doesn't have a parent node, thus use the roots ID as the parent node.
			// that's not problematic because we assume that there aren't any self loops
			make_tree(g,*r_it, *r_it, forest);
		}

		// build parallel trees
		for(typename BorderComponents::const_iterator b_it=border_components.begin(),
													  b_it_end=border_components.end();
			b_it!=b_it_end; b_it++)
		{
			if(b_it->is_global_max)
			{// this is the root part of the global tree - this trees root starts at a local node

				// the root doesn't have a parent node, thus use the roots ID as the parent node.
				// that's not problematic because we assume that there aren't any self loops
				make_tree(g, b_it->root_edge.n1, b_it->root_edge.n1, forest);

				roots_of_parallel_trees.push_back(b_it->root_edge.n1);
			}
			else
			{// this is not the root of the global tree, thus we use the ghost_node (connecting to the parent tree)
							 // as the local root

				// ghost_node is the root, but we start building the tree beginning at local_root (the child node).
				// the other outgoing edges of ghost_node will be added later, this way we don't move up the
				// tree if ghost_node is not a pure root
				make_tree(g, b_it->root_edge.n1, b_it->root_edge.n2, forest);

				// the name visited actually means: is_root?
				// we only want to add a root node once to the list of roots
				if(!forest.visited(b_it->root_edge.n2))
				{
					roots_of_parallel_trees.push_back(b_it->root_edge.n2);
					forest.set_visited(b_it->root_edge.n2, true);
				}
			}
		}
//		for(typename ConnectionComponents::const_iterator c_it=connection_components.begin(),
//														  c_it_end=connection_components.end();
//			c_it!=c_it_end; c_it++)
//		{
//			if(g.local_vertex_id_to_global_id(c_it->root_edge_local_node) == c_it->root_candidate)
//			{// this is the root part of the global tree - this trees root starts at a local node
//
//				// the root doesn't have a parent node, thus use the roots ID as the parent node.
//				// that's not problematic because we assume that there aren't any self loops
//				make_tree(g, c_it->root_edge_local_node, c_it->root_edge_local_node, forest);
//
//				// this actually set the status to is_root of local node
//				forest.set_visited(c_it->root_edge_local_node, true);
//				roots_of_parallel_trees.push_back(c_it->root_edge_local_node);
//			}
//			else
//			{// this is not the root of the global tree, thus we use the ghost_node (connecting to the parent tree)
//			 // as the local root
//
//				// ghost_node is the root, but we start building the tree beginning at local_root (the child node).
//				// the other outgoing edges of ghost_node will be added later, this way we don't move up the
//				// tree if ghost_node is not a pure root
//				make_tree(g, c_it->root_edge_local_node, c_it->root_edge_ghost_node, forest);
//
//				// the name visited actually means: is_root?
//				// we only want to add a root node once to the list of roots
//				if(!forest.visited(c_it->root_edge_ghost_node))
//				{
//					roots_of_parallel_trees.push_back(c_it->root_edge_ghost_node);
//					forest.set_visited(c_it->root_edge_ghost_node, true);
//				}
//			}
//		}

	}





	void compute_matching_of_local_trees(const Graph &g,
										 const std::vector<NodeIDType> &roots_of_local_trees,
 	   	   	   	    					 std::vector<Edge> &result)
	{
		// compute the tree matchings for the remaining local trees
		for(NodeIDType n = next_local_tree_pos; n<roots_of_local_trees.size(); n++)
		{
			NodeIDType current_node = roots_of_local_trees[n];

			tree_max_weighted_matching(current_node, 0);

			add_matching_edges(current_node, false, result);
		}

		// we've computed the matchings of all local trees -> ensures no multiple
		// computation of matchings of local trees
		next_local_tree_pos = roots_of_local_trees.size();
	}


	/*!
	 * \brief Returns the local ID of the given ghost ID
	 *
	 * We use the global ID of an adjacent local node of the ghost node
	 * as a hint where to search for the local ID of the ghost node.
	 *
	 * Runtime: O(deg(adjacent_global_id)) worst case O(log(num_ghost_vertices))
	 *
	 * @param g
	 * @param global_ghost_id
	 * @param adjacent_global_id
	 * @return
	 */
	NodeIDType get_local_id_of_ghost_vertex(const NodeIDType global_ghost_id,
											const NodeIDType adjacent_global_id)
	{
		NodeIDType adjacent_local_id = g.global_vertex_id_to_local_id_of_local_vertex(adjacent_global_id);


		if(forest.out_degree(adjacent_local_id)>log_num_ghost_vertices)
		{
			calls_to_g_for_gid_to_lid++;
		 	return g.global_vertex_id_to_local_id_of_ghost_vertex(global_ghost_id);
		}


		for(typename Forest::incident_edge_iterator
				e_it=forest.begin_incident_edges(adjacent_local_id),
				end = forest.end_incident_edges(adjacent_local_id);
				e_it!=end;
				e_it++)
		{
			if(g.local_vertex_id_to_global_id(e_it->second)==global_ghost_id)
			{
				return e_it->second;
			}
		}

		// couldn't find the node -> e.g. the given adjacent node id was wrong
		return -1;
	}




	NodeIDType get_pending_msgs_of_root(const NodeIDType root)
	{
		NodeIDType result = 0;

		for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(root),
				 	 	 	 	 	 	 	 	 	end = forest.end_incident_edges(root);
				e_it!=end; e_it++)
		{
			if(g.is_ghost(e_it->second))
			{// we don't check outgoing edges of ghost nodes (if they exist),
			 // they belong to a different tree. the ghost node is the root
			 // of those trees
				result++;
			}
			else
			{
				result += get_pending_msgs_of_root(e_it->second);
			}
		}

		return result;
	}


	void set_pending_msgs_of_roots(const Graph &g,
								   const std::vector<NodeIDType> &roots_of_parallel_trees)
	{
		// assuming that pending_msgs_of_root is zero for every node
		for(typename std::vector<NodeIDType>::const_iterator n_it=roots_of_parallel_trees.begin(),
															 n_it_end=roots_of_parallel_trees.end();
			n_it!=n_it_end; n_it++)
		{
			pending_msgs_of_root[*n_it] = get_pending_msgs_of_root(*n_it);
		}
	}


	/*!
	 * \brief Fills the subtree-table which is used to compute the optimal matching
	 * of the parallel tree and sends messages to the parent-tree (if there's one).
	 *
	 * @param root	The root-node of the subtree who's dynamic table
	 * 				is filled.
	 */
	void fill_subtree_table_for_given_root_and_add_messages_for_parent_tree(const NodeIDType root)
	{
		if(g.is_ghost(root))
		{// not a root tree, thus we have to send some messages to the parent tree
			SubtreeResultMSG msg;

			int proc_of_root = g.get_proc_of_ghost_vertex(root);

			if(proc_requests[proc_of_root].size() == proc_bottom_up_buffers[proc_of_root].size())
			{// both have the same size, thus the last buffer of proc_msg_buffers[proc_of_root]
			 // is from an older communication -> we need a new buffer for the current communication
				proc_bottom_up_buffers[proc_of_root].push_back(std::vector<SubtreeResultMSG>());
			}


			for(typename Forest::incident_edge_iterator
					e_it=forest.begin_incident_edges(root),
					end = forest.end_incident_edges(root);
					e_it!=end; e_it++)
			{
				tree_max_weighted_matching(e_it->second,
										   g.get_edge(e_it->first).weight);

				msg.hint_node    = g.local_vertex_id_to_global_id(root);
				msg.subtree_root = g.local_vertex_id_to_global_id(e_it->second);
				msg.incoming_matched = subtree_table[e_it->second].first.first;
				msg.incoming_not_matched = subtree_table[e_it->second].first.second;
				msg.uses_outgoing_edge = subtree_table[e_it->second].second!=(NodeIDType)-1;


				proc_bottom_up_buffers[proc_of_root].back().push_back(msg);
			}
		}
		else
		{// this is a root-tree, we don't have to send any messages
			// this root doensn't have any incoming edge, thus incoming weight=0
			tree_max_weighted_matching(root, 0);
		}
	}


	void isend_subtree_messages()
	{
		/// \todo use a list of active partners -> don't have check each process
		for(int p=0; p<procs; p++)
		{
			if(proc_requests[p].size() != proc_bottom_up_buffers[p].size())
			{// there are new messages for p
				proc_requests[p].push_back(MPI_Request());

				MPI_Isend(&(proc_bottom_up_buffers[p].back()[0]),
						 proc_bottom_up_buffers[p].back().size(),
						 SubtreeMSGType, p, bottom_up_tag,
						 MPI_COMM_WORLD, &(proc_requests[p].back()) );

				number_of_mpi_messages_sent++;
			}
		}

	}


	void isend_top_down_messages()
	{
		/// \todo use a list of active partners -> don't have check each process
		for(int p=0; p<procs; p++)
		{
			if(proc_requests[p].size() != proc_top_down_buffers[p].size())
			{// there are new messages for p
				proc_requests[p].push_back(MPI_Request());

				MPI_Isend(&(proc_top_down_buffers[p].back()[0]),
						 proc_top_down_buffers[p].back().size(),
						 TopDownMSGType, p, top_down_tag,
						 MPI_COMM_WORLD, &(proc_requests[p].back()) );

				number_of_mpi_messages_sent++;
			}
		}

	}


	void receive_message_and_fill_subtree_table(MPI_Status &status, NodeIDType &remaining_roots)
	{
		// i might receive a msg different than the one corresponding to status
		// because the description by src_proc (and this proc) isn't necessarily
		// unique. But this doesn't matter because we're not loosing the other
		// message

		/// \todo if the src process send 2 messages to this process it's not clear
		/// which message i'll receive, wrong buffer!!! use a counter for the tag??
		/// or make it round based
		/// same problem occures in top down

		int src_proc = status.MPI_SOURCE;
		int msg_size;

		MPI_Get_count(&status, SubtreeMSGType, &msg_size);

		std::vector<SubtreeResultMSG> msg_buffer(msg_size);

		MPI_Status tmp_status;
		int recv_size;

		MPI_Recv(&msg_buffer[0], msg_size, SubtreeMSGType, src_proc, bottom_up_tag, MPI_COMM_WORLD, &tmp_status /*MPI_STATUS_IGNORE*/);
		number_of_mpi_messages_received++;

		MPI_Get_count(&tmp_status, SubtreeMSGType, &recv_size);

		for(int i=0; i<msg_size; i++)
		{
			NodeIDType local_ghost_id = get_local_id_of_ghost_vertex(msg_buffer[i].subtree_root, msg_buffer[i].hint_node);

			subtree_table[local_ghost_id].first.first = msg_buffer[i].incoming_matched;
			subtree_table[local_ghost_id].first.second = msg_buffer[i].incoming_not_matched;

			// if no outgoing edge is used, indicate it with value -1.
			// This helps us to decide wether the node is matched or not,
			// without any extra messages from the subtree.
			subtree_table[local_ghost_id].second = msg_buffer[i].uses_outgoing_edge ? 0 : -1;

			NodeIDType root = forest.get_root_of_tree(local_ghost_id);

			pending_msgs_of_root[root]--;

			if(pending_msgs_of_root[root]==0)
			{
				fill_subtree_table_for_given_root_and_add_messages_for_parent_tree(root);
				remaining_roots--;
			}
		}
	}


	/*!
	 * \brief This functions probes for a MPI-message corresponding to src, tag
	 * and communicator. This function will also compute matchings of local
	 * trees while waiting for such a message.
	 *
	 * This function won't return until there's such a message.
	 *
	 * The template argument rounds specifies the number of matching of local
	 * trees that are computed for each check of such a message.
	 *
	 * @param src		ID of the source process within the communicator
	 * @param tag		The tag of the message
	 * @param comm		The communicator for this communication
	 * @param status	Status object for the probed message
	 * @param result	Vector to store the matching result of the local tree
	 * 					matchings.
	 */
	template <unsigned int rounds = 1>
	void probe_and_compute_local_tree_matching(int src, int tag, MPI_Comm comm,
											   MPI_Status *status,
											   std::vector<Edge> &result)
	{
		int msg_pending = false;

		while(!msg_pending)
		{
			MPI_Iprobe(src, tag, comm, &msg_pending, status);

			if(msg_pending)
			{
				return;
			}

			// hopefully the compiler will do some loop unrolling
			for(unsigned int i=0; i<rounds; i++)
			{
				if(next_local_tree_pos < roots_of_local_trees_ptr->size())
				{
					// no pending messages but use our time to do some computations
					NodeIDType local_root =
							(*roots_of_local_trees_ptr)[next_local_tree_pos];

					tree_max_weighted_matching(local_root, 0);

					add_matching_edges(local_root, false, result);

					next_local_tree_pos++;
				}
			}
		}
	}


	// this should help to send several msgs in a single MPI_MSG
	// needs a better name
	void wait_for_pending_messages_and_fill_subtree_table(NodeIDType &remaining_roots,
														  std::vector<Edge> &result)
	{
		MPI_Status status;

//		MPI_Probe(MPI_ANY_SOURCE, bottom_up_tag, MPI_COMM_WORLD, &status);
		probe_and_compute_local_tree_matching<1>(MPI_ANY_SOURCE,
											  bottom_up_tag,
											  MPI_COMM_WORLD,
											  &status,
											  result);

		receive_message_and_fill_subtree_table(status, remaining_roots);

		int msg_is_waiting;

		MPI_Iprobe(MPI_ANY_SOURCE, bottom_up_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);

		while(msg_is_waiting)
		{
			receive_message_and_fill_subtree_table(status, remaining_roots);

			MPI_Iprobe(MPI_ANY_SOURCE, bottom_up_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);
		}

	}



	void fill_subtree_table_bottom_up_parallel(const Graph &g,
											   const std::vector<NodeIDType> &roots_of_parallel_trees,
											   std::vector<Edge> &result)
	{
		set_pending_msgs_of_roots(g, roots_of_parallel_trees);

		NodeIDType remaining_roots = roots_of_parallel_trees.size();
//		NodeIDType msgs_to_receive = 0;

		// initially send the results for the leave trees
		for(NodeIDType root_pos=0; root_pos<roots_of_parallel_trees.size(); root_pos++)
		{
			const NodeIDType current_root = roots_of_parallel_trees[root_pos];

//			msgs_to_receive += pending_msgs_of_root[current_root];

			if(pending_msgs_of_root[current_root]==0)
			{
				fill_subtree_table_for_given_root_and_add_messages_for_parent_tree(current_root);

				remaining_roots--;
			}

		}

		// send initial messages
		isend_subtree_messages();

		// receive messages and compute tree matchings if all information is
		// available. Also send new messages afterwards
		while(remaining_roots>0)
		{
			wait_for_pending_messages_and_fill_subtree_table(remaining_roots,
															 result);
			isend_subtree_messages();

			/// \todo check for completed messages (to free memory)
		}

		wait_for_message_completion_and_reset_buffers(proc_requests, proc_bottom_up_buffers);
	}

	struct ActiveNodeInfo
	{
		NodeIDType active_node;
		bool parent_matched, used_incoming;

		ActiveNodeInfo()
		{}

		ActiveNodeInfo(NodeIDType active_node, bool parent_matched, bool used_incoming)
		: active_node(active_node), parent_matched(parent_matched), used_incoming(used_incoming)
		{}
	};


	/*!
	 *
	 * @param root 					Must not be a ghost node
	 * @param used_incoming_edge
	 * @param result
	 * @param matched_nodes
	 */
	void add_matched_edges_of_subtree_and_send_top_down_messages(const NodeIDType root,
																 const bool used_incoming_edge,
																 std::vector<Edge> &result,
																 std::vector<NodeIDType> &matched_nodes)
	{
		// root must not be a ghost node!!!
		std::list<ActiveNodeInfo> active_nodes;

		// because root isn't a ghost node, we don't have to know if the parent is matched
		active_nodes.push_back(ActiveNodeInfo(root, used_incoming_edge, used_incoming_edge));


		while(!active_nodes.empty())
		{
			NodeIDType current_node = active_nodes.front().active_node;
			bool used_incoming_edge = active_nodes.front().used_incoming;
			bool parent_matched		= active_nodes.front().parent_matched;
			active_nodes.pop_front();

			if(g.is_ghost(current_node))
			{
				int proc_of_current_node = g.get_proc_of_ghost_vertex(current_node);

				if(proc_requests[proc_of_current_node].size() == proc_top_down_buffers[proc_of_current_node].size())
				{// both have the same size, thus the last buffer of proc_top_down_buffers[proc_of_current_node]
				 // is from an older communication -> we need a new buffer for the current communication
					proc_top_down_buffers[proc_of_current_node].push_back(std::vector<TopDownMSG>());
				}

				if(! used_incoming_edge && subtree_table[current_node].second!=(NodeIDType)-1)
				{// current_node uses an outgoing edge on the subtree
					matched_nodes.push_back(current_node);
				}

				TopDownMSG msg;

				if(used_incoming_edge)
				{
					msg.msg_info = TopDownMSG::use_incoming_edge;
				}
				// we have to provide extra information for the other process.
				// the parent might not be matched
				// \todo actually that's not true - all internal vertices are matched
				else if(parent_matched)
				{
					msg.msg_info = TopDownMSG::parent_matched;
				}
				else
				{
					msg.msg_info = TopDownMSG::parent_not_matched;
				}

				msg.tgt_node = g.local_vertex_id_to_global_id(current_node);

				// add the message to the send buffers
				proc_top_down_buffers[proc_of_current_node].back().push_back(msg);
			}
			else
			{
				if(used_incoming_edge)
				{
					for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(current_node),
																end = forest.end_incident_edges(current_node);
							e_it!=end; e_it++)
					{
						active_nodes.push_back(ActiveNodeInfo(e_it->second, true, false));
					}
				}
				else
				{
					for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(current_node),
																end = forest.end_incident_edges(current_node);
							e_it!=end; e_it++)
					{
						if(subtree_table[current_node].second == e_it->second)
						{
							result.push_back(g.get_edge(e_it->first));

							active_nodes.push_back(ActiveNodeInfo(e_it->second, true, true));
						}
						else
						{
							active_nodes.push_back(ActiveNodeInfo(e_it->second, subtree_table[current_node].second!=(NodeIDType)-1, false));
						}
					}
				}
			}
		}
	}


	void receive_message_and_add_matched_edges(MPI_Status &status,
											   NodeIDType &number_of_incoming_messages,
											   std::vector<Edge> &result,
											   std::vector<NodeIDType> &matched_nodes)
	{
		// i might receive a msg different than the one corresponding to status
		// because the description by src_proc (and this proc) isn't necessarily
		// unique. But this doesn't matter because we're not loosing the other
		// message

		int src_proc = status.MPI_SOURCE;
		int msg_size;

		MPI_Get_count(&status, TopDownMSGType, &msg_size);

		std::vector<TopDownMSG> msg_buffer(msg_size);

		MPI_Recv(&msg_buffer[0], msg_size, TopDownMSGType, src_proc, top_down_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		number_of_mpi_messages_received++;

		// we received msg_size many incoming messages
		number_of_incoming_messages -= msg_size;

		for(int i=0; i<msg_size; i++)
		{
			NodeIDType root_id = g.global_vertex_id_to_local_id_of_local_vertex(msg_buffer[i].tgt_node);
			bool use_incoming_edge = (msg_buffer[i].msg_info == TopDownMSG::use_incoming_edge);

			if(use_incoming_edge)
			{// the incoming edge is matched, but it was added to the result by the
			 // parent partner, thus we don't add the edge to the result, but we still
			 // have to indicate that both end nodes are matched!

				// in this case the root of the tree is also the parent. tgt_nodes are
				// always children of the root of the tree!!
				NodeIDType parent = forest.get_root_of_tree(root_id);
				matched_nodes.push_back(parent);
				matched_nodes.push_back(root_id);
			}
			else if(msg_buffer[i].msg_info == TopDownMSG::parent_matched)
			{// the parent is a matched node
				// in this case the root of the tree is also the parent. tgt_nodes are
				// always children of the root of the tree!!
				NodeIDType parent = forest.get_root_of_tree(root_id);
				matched_nodes.push_back(parent);
			}

			// we don't have to add the incoming edge to the result, that's already done
			// by the parent tree on the other process. we don't want double entries of
			// cross edges

			add_matched_edges_of_subtree_and_send_top_down_messages(root_id, use_incoming_edge,
																	result, matched_nodes);
		}
	}


	void wait_for_pending_top_down_messages_and_add_matched_edges(NodeIDType &number_of_incoming_messages,
																  std::vector<Edge> &result,
																  std::vector<NodeIDType> &matched_nodes)
	{
		MPI_Status status;

//		MPI_Probe(MPI_ANY_SOURCE, top_down_tag, MPI_COMM_WORLD, &status);
		probe_and_compute_local_tree_matching<1>(MPI_ANY_SOURCE,
											  	 top_down_tag,
											  	 MPI_COMM_WORLD,
											  	 &status,
											  	 result);

		receive_message_and_add_matched_edges(status, number_of_incoming_messages,
											  result, matched_nodes);

		int msg_is_waiting;

		MPI_Iprobe(MPI_ANY_SOURCE, top_down_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);

		while(msg_is_waiting)
		{
			receive_message_and_add_matched_edges(status, number_of_incoming_messages,
												  result, matched_nodes);

			MPI_Iprobe(MPI_ANY_SOURCE, top_down_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);
		}

	}




#ifdef MORE_INFORMATION


	void get_subtree_info(const NodeIDType root,
						  NodeIDType &number_of_nodes,
						  NodeIDType &depth,
						  NodeIDType &number_of_components,
						  NodeIDType &max_component_depth)
	{
		number_of_nodes = 1;
		depth = 0;

		for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(root),
													end = forest.end_incident_edges(root);
				e_it!=end; e_it++)
		{
			if(g.is_ghost(e_it->second))
			{
				number_of_nodes += subtree_info[e_it->second].number_of_nodes;
				if(subtree_info[e_it->second].depth>depth)
				{
					depth = subtree_info[e_it->second].depth;
				}

				number_of_components += subtree_info[e_it->second].number_of_components;
				if(subtree_info[e_it->second].component_depth>max_component_depth)
				{
					max_component_depth = subtree_info[e_it->second].component_depth;
				}
			}
			else
			{
				NodeIDType sub_number_of_nodes, sub_depth;
				get_subtree_info(e_it->second, sub_number_of_nodes, sub_depth, number_of_components, max_component_depth);
				number_of_nodes += sub_number_of_nodes;

				if(sub_depth>depth)
				{
					depth = sub_depth;
				}
			}
		}

		depth += (forest.out_degree(root)>0);
	}


	void compute_subtree_information_and_add_messages_for_parent(const NodeIDType root)
	{
		if(g.is_ghost(root))
		{// not a root tree, thus we have to send some messages to the parent tree
			bool not_added = true;

			for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(root),
														end = forest.end_incident_edges(root);
					e_it!=end; e_it++)
			{
				NodeIDType number_of_nodes = 0;
				NodeIDType depth = 0;
				NodeIDType number_of_components = 0;
				NodeIDType component_depth = 0;

				get_subtree_info(e_it->second,
								 number_of_nodes,
								 depth,
								 number_of_components,
								 component_depth);

				number_of_components += not_added;
				not_added = false,

				component_depth += (get_pending_msgs_of_root(e_it->second) > 0);


				ParallelTreeInfo pt_info;
				pt_info.depth = depth;
				pt_info.number_of_nodes = number_of_nodes;
				pt_info.number_of_components = number_of_components;
				pt_info.component_depth = component_depth;

				pt_info.subtree_root = g.local_vertex_id_to_global_id(e_it->second);
				pt_info.hint_node = g.local_vertex_id_to_global_id(root);

				int proc_of_current_node = g.get_proc_of_ghost_vertex(root);

				if(proc_subtree_info_buffers[proc_of_current_node].size() == proc_requests[proc_of_current_node].size())
				{// both have the same size, thus the last buffer of proc_top_down_buffers[proc_of_current_node]
				 // is from an older communication -> we need a new buffer for the current communication
					proc_subtree_info_buffers[proc_of_current_node].push_back(std::vector<ParallelTreeInfo>());
				}

				proc_subtree_info_buffers[proc_of_current_node].back().push_back(pt_info);
			}
		}
		else
		{// this is a root-tree, we don't have to send any messages. store the info about this parallel tree
			NodeIDType number_of_nodes = 0;
			NodeIDType depth = 0;
			NodeIDType number_of_components = 0;
			NodeIDType component_depth = 0;

			get_subtree_info(root,
							 number_of_nodes,
							 depth,
							 number_of_components,
							 component_depth);

			number_of_components++;
			component_depth++;

			outfile_parallel_tree_info << proc_id << " "
									   << number_of_nodes << " "
									   << depth << " "
									   << number_of_components << " "
									   << component_depth << std::endl;
		}
	}




	void receive_message_and_add_parallel_tree_infos(MPI_Status &status, NodeIDType &remaining_roots)
	{
		// i might receive a msg different than the one corresponding to status
		// because the description by src_proc (and this proc) isn't necessarily
		// unique. But this doesn't matter because we're not loosing the other
		// message

		/// \todo if the src process send 2 messages to this process it's not clear
		/// which message i'll receive, wrong buffer!!! use a counter for the tag??
		/// or make it round based
		/// same problem occures in top down

		int src_proc = status.MPI_SOURCE;
		int msg_size;

		MPI_Get_count(&status, ParallelTreeInfoMSGType, &msg_size);

		std::vector<ParallelTreeInfo> msg_buffer(msg_size);


		MPI_Recv(&msg_buffer[0], msg_size, ParallelTreeInfoMSGType, src_proc, parallel_tree_info_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for(int i=0; i<msg_size; i++)
		{
			NodeIDType local_ghost_id = get_local_id_of_ghost_vertex(msg_buffer[i].subtree_root, msg_buffer[i].hint_node);

			subtree_info[local_ghost_id] = msg_buffer[i];

			NodeIDType root = forest.get_root_of_tree(local_ghost_id);

			pending_msgs_of_root[root]--;

			if(pending_msgs_of_root[root]==0)
			{
				compute_subtree_information_and_add_messages_for_parent(root);
				remaining_roots--;
			}
		}
	}


	// this should help to send several msgs in a single MPI_MSG
	// needs a better name
	void wait_for_pending_messages_and_add_parallel_tree_infos(NodeIDType &remaining_roots)
	{
		MPI_Status status;

		MPI_Probe(MPI_ANY_SOURCE, parallel_tree_info_tag, MPI_COMM_WORLD, &status);

		receive_message_and_add_parallel_tree_infos(status, remaining_roots);

		int msg_is_waiting;

		MPI_Iprobe(MPI_ANY_SOURCE, parallel_tree_info_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);

		while(msg_is_waiting)
		{
			receive_message_and_add_parallel_tree_infos(status, remaining_roots);

			MPI_Iprobe(MPI_ANY_SOURCE, parallel_tree_info_tag, MPI_COMM_WORLD, &msg_is_waiting, &status);
		}

	}


	void isend_parallel_tree_info_messages()
	{
		/// \todo use a list of active partners -> don't have check each process
		for(int p=0; p<procs; p++)
		{
			if(proc_requests[p].size() != proc_subtree_info_buffers[p].size())
			{// there are new messages for p
				proc_requests[p].push_back(MPI_Request());

				MPI_Isend(&(proc_subtree_info_buffers[p].back()[0]),
						 proc_subtree_info_buffers[p].back().size(),
						 ParallelTreeInfoMSGType, p, parallel_tree_info_tag,
						 MPI_COMM_WORLD, &(proc_requests[p].back()) );
			}
		}

	}



	void get_parallel_tree_information(const std::vector<NodeIDType> &roots_of_parallel_trees)
	{
		for(typename std::vector<NodeIDType>::const_iterator nit=roots_of_parallel_trees.begin(),
															 nit_end=roots_of_parallel_trees.end();
				nit!=nit_end; nit++)
		{
			pending_msgs_of_root[*nit] = 0;
		}


		set_pending_msgs_of_roots(g, roots_of_parallel_trees);

		NodeIDType remaining_roots = roots_of_parallel_trees.size();
		NodeIDType msgs_to_receive = 0;

		// initially send the results for the leave trees
		for(NodeIDType root_pos=0; root_pos<roots_of_parallel_trees.size(); root_pos++)
		{
			const NodeIDType current_root = roots_of_parallel_trees[root_pos];

			msgs_to_receive += pending_msgs_of_root[current_root];

			if(pending_msgs_of_root[current_root]==0)
			{
				compute_subtree_information_and_add_messages_for_parent(current_root);

				remaining_roots--;
			}

		}

		// send initial messages
		isend_parallel_tree_info_messages();

		// receive messages and compute tree matchings if all information is
		// available. Also send new messages afterwards
		while(remaining_roots>0)
		{
			wait_for_pending_messages_and_add_parallel_tree_infos(remaining_roots);
			isend_parallel_tree_info_messages();

			/// \todo check for completed messages (to free memory)
		}

		wait_for_message_completion_and_reset_buffers(proc_requests, proc_subtree_info_buffers);


		for(typename std::vector<NodeIDType>::const_iterator nit=roots_of_parallel_trees.begin(),
															 nit_end=roots_of_parallel_trees.end();
				nit!=nit_end; nit++)
		{
			pending_msgs_of_root[*nit] = 0;
		}
	}

#endif





	void get_matching_of_parallel_forest(const std::vector<NodeIDType> &roots_of_parallel_trees,
										 std::vector<Edge> &result,
										 std::vector<NodeIDType> &matched_nodes)
	{
		NodeIDType number_of_incoming_messages = 0;

		// send initial top down messages and compute number os messages to be received
		for(NodeIDType pos=0; pos<roots_of_parallel_trees.size(); pos++)
		{
			NodeIDType root = roots_of_parallel_trees[pos];

			if(g.is_local(root))
			{// that's a root tree we don't receive any messages for this tree
			 // but we might have to send a few messages
				// no incoming edge is used - because there's none
				add_matched_edges_of_subtree_and_send_top_down_messages(root, false, result, matched_nodes);
			}
			else
			{// have to wait for incoming messages
				number_of_incoming_messages += forest.out_degree(root);
			}
		}

		// send initial messages
		isend_top_down_messages();


		// receive messages and compute tree matchings if all information is
		// available. Also send new messages afterwards
		while(number_of_incoming_messages>0)
		{
			wait_for_pending_top_down_messages_and_add_matched_edges(number_of_incoming_messages,
																	 result,
																	 matched_nodes);
			isend_top_down_messages();

			/// \todo check for completed messages (to free memory)
		}

		wait_for_message_completion_and_reset_buffers(proc_requests, proc_top_down_buffers);
	}


	/*!
	 *
	 *
	 * @param g
	 * @param roots_of_parallel_trees
	 * @param result
	 */
	void compute_matching_of_parallel_trees(const std::vector<NodeIDType> &roots_of_parallel_trees,
 	   	   	   	    						std::vector<Edge> &result,
 	   	   	   	    						std::vector<NodeIDType> &matched_nodes)
	{
#ifdef MORE_INFORMATION
		double start, end;

		start = MPI_Wtime();
#endif

		fill_subtree_table_bottom_up_parallel(g, roots_of_parallel_trees, result);

#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run fill_subtree_table_bottom_up_parallel=" << (end-start) << " sec" << std::endl;

		start = MPI_Wtime();
#endif

		get_matching_of_parallel_forest(roots_of_parallel_trees, result, matched_nodes);

#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run get_matching_of_parallel_forest=" << (end-start) << " sec" << std::endl;
#endif

	}


	void get_parallel_tree_matching(const std::vector<NodeIDType> &roots_of_local_trees,
	  	  	  	   	   	   	   	    const std::vector<NodeIDType> &roots_of_parallel_trees,
	  	  	  	   	   	   	   	    std::vector<Edge> &result,
	  	  	  	   	   	   	   	    std::vector<NodeIDType> &matched_nodes)
	{
		// compute the matching of parallel trees
		compute_matching_of_parallel_trees(roots_of_parallel_trees, result, matched_nodes);


#ifdef MORE_INFORMATION
		double start, end;

		start = MPI_Wtime();
#endif

		// compute the matching of local_trees
		compute_matching_of_local_trees(g, roots_of_local_trees, result);

#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run compute_matching_of_local-trees=" << (end-start) << " sec" << std::endl;
#endif
	}



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
	inline void
	tree_max_weighted_matching(const NodeIDType &tgt_node,
							   const WeightType &weight_to_node)
	{
		std::list<NodeIDType> active_nodes;
		std::vector<std::pair<NodeIDType, WeightType> > bfs_traversal;
		active_nodes.push_back(tgt_node);
		bfs_traversal.push_back(std::make_pair(tgt_node, weight_to_node));

		// create BFS-traversal
		while(!active_nodes.empty())
		{
			NodeIDType current_node   = active_nodes.front();
			active_nodes.pop_front();

			for(typename Forest::incident_edge_iterator
					e_it=forest.begin_incident_edges(current_node),
					end = forest.end_incident_edges(current_node);
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

			if(forest.out_degree(current_node)==0)
			{
				// values for ghost_nodes are set during communication phase
				// so we don't have to set them here
				if(! g.is_ghost(current_node))
				{
					subtree_table[current_node].first = std::make_pair(current_weight, 0);
					subtree_table[current_node].second = -1; // there are no outgoing edges to be used!!
				}
			}
			else
			{
				// case: edge to this subtree is matched
				WeightType matching_weight1 = current_weight;

				for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(current_node),
															end = forest.end_incident_edges(current_node);
						e_it!=end; e_it++)
				{
					// already computed the local result from each following node -> traversing the tree bottom up
					matching_weight1 += subtree_table[e_it->second].first.second;
				}

				// case: edge to this subtree is not matched
				WeightType matching_weight2 = 0;

				WeightType old_diff = 0;
				subtree_table[current_node].second = -1; // set the outgoing partner initially to -1,
														 // in case we're not using any outgoing partner.

				for(typename Forest::incident_edge_iterator
						e_it=forest.begin_incident_edges(current_node),
						end = forest.end_incident_edges(current_node);
						e_it!=end; e_it++)
				{
					const std::pair<WeightType, WeightType> current_pair = subtree_table[e_it->second].first;
					const WeightType current_diff = current_pair.first-current_pair.second;

					if(current_diff>old_diff)
					{
						matching_weight2 += current_pair.first - old_diff;
						subtree_table[current_node].second = e_it->second; // set the target-node, who's edge has been used by the matching
						old_diff = current_diff;
					}
					else
					{
						matching_weight2 += current_pair.second;
					}
				}

				subtree_table[current_node].first = std::make_pair(matching_weight1, matching_weight2);
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
	inline void add_matching_edges(const NodeIDType &current_node,
								   const bool &used_incoming_edge,
								   std::vector<Edge> &result)
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
				for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(current_node),
															end = forest.end_incident_edges(current_node);
						e_it!=end; e_it++)
				{
					active_nodes.push_back(std::make_pair(e_it->second, false));
				}
			}
			else
			{
				for(typename Forest::incident_edge_iterator e_it=forest.begin_incident_edges(current_node),
															end = forest.end_incident_edges(current_node);
						e_it!=end; e_it++)
				{
					if(subtree_table[current_node].second == e_it->second)
					{
						result.push_back(g.get_edge(e_it->first));

						active_nodes.push_back(std::make_pair(e_it->second, true));
					}
					else
					{
						active_nodes.push_back(std::make_pair(e_it->second, false));
					}
				}
			}

		}
	}




	void reset(const std::vector<EdgeIDType> &local_edges,
			   const std::vector<EdgeIDType> &cross_edges,
			   const std::vector<NodeIDType> &parallel_tree_roots)
	{
		for(typename std::vector<EdgeIDType>::const_iterator eit=local_edges.begin(),
															 eit_end=local_edges.end();
			eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			uf.reset(e.n1);
			uf.reset(e.n2);

			connection_component_ptr_of_vertex[e.n1] = NULL;
			connection_component_ptr_of_vertex[e.n2] = NULL;

			forest.clear(e.n1);
			forest.clear(e.n2);
		}

		for(typename std::vector<EdgeIDType>::const_iterator eit=cross_edges.begin(),
															 eit_end=cross_edges.end();
				eit!=eit_end; eit++)
		{
			const Edge &e = g.get_edge(*eit);

			uf.reset(e.n1); // we have to reset all local nodes
			connection_component_ptr_of_vertex[e.n1] = NULL;

			forest.clear(e.n1);
			forest.clear(e.n2);
		}

		for(typename std::vector<NodeIDType>::const_iterator nit=parallel_tree_roots.begin(),
															 nit_end=parallel_tree_roots.end();
				nit!=nit_end; nit++)
		{
			pending_msgs_of_root[*nit] = 0;
		}
	}

public:


	ParallelForestMatching(const Graph &g)
	:	g(g),
		uf(g.num_local_vertices()),
	 	connection_component_ptr_of_vertex(g.num_local_vertices(), NULL),
		num_local_vertices(g.num_local_vertices()),
		num_ghost_vertices(g.num_local_ghost_vertices()),
		log_num_ghost_vertices(std::log(num_ghost_vertices)/std::log(2.)),
		forest(num_local_vertices+num_ghost_vertices),
		pending_msgs_of_root(num_local_vertices+num_ghost_vertices),
		subtree_table(num_local_vertices+num_ghost_vertices),
		number_of_mpi_messages_sent(0),
		number_of_mpi_messages_received(0),
		calls_to_g_for_gid_to_lid(0)
	{
		MPI_Type_contiguous(sizeof(SubtreeResultMSG), MPI_BYTE, &SubtreeMSGType);
		MPI_Type_commit(&SubtreeMSGType);

		MPI_Type_contiguous(sizeof(TopDownMSG), MPI_BYTE, &TopDownMSGType);
		MPI_Type_commit(&TopDownMSGType);

		MPI_Type_contiguous(sizeof(ComponentMSG), MPI_BYTE, &ComponentMSGType);
		MPI_Type_commit(&ComponentMSGType);

		MPI_Type_contiguous(sizeof(BorderComponentMSG), MPI_BYTE, &BorderComponentMSGType);
		MPI_Type_commit(&BorderComponentMSGType);

		// set number of processes
		MPI_Comm_size(MPI_COMM_WORLD, &procs);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

		proc_requests.resize(procs);
		proc_bottom_up_buffers.resize(procs);
		proc_top_down_buffers.resize(procs);

		proc_is_active_src.resize(procs, false);

		std::stringstream outfilename;
		outfilename << "out_remaining_nodes_" << proc_id;

#ifdef MORE_INFORMATION
		MPI_Type_contiguous(sizeof(ParallelTreeInfo), MPI_BYTE, &ParallelTreeInfoMSGType);
		MPI_Type_commit(&ParallelTreeInfoMSGType);

		subtree_info.resize(num_local_vertices+num_ghost_vertices);
		proc_subtree_info_buffers.resize(procs);

		outfilename.str("");
		outfilename << "parallel_tree_infos_" << proc_id;
		outfile_parallel_tree_info.open(outfilename.str().c_str());


		outfilename.str("");
		outfilename << "std_out_" << proc_id;
		outfile_std_out.open(outfilename.str().c_str());
#endif
	}

	~ParallelForestMatching()
	{
#ifdef MORE_INFORMATION
		outfile_parallel_tree_info.close();
		outfile_std_out.close();
#endif
	}

	/*!
	 *
	 * @param g
	 * @param local_edges
	 * @param cross_edges
	 * @param result
	 * @param matched_nodes a list that is used for matched nodes that aren't incident
	 * 		  				to any matched edges that locally have been added to the result
	 */
	void get_matching(const Graph &g,
					  const std::vector<EdgeIDType> &local_edges,
				  	  const std::vector<EdgeIDType> &cross_edges,
				  	  std::vector<Edge> &result,
				  	  std::vector<NodeIDType> &matched_nodes)
	{
	//	MPI_Type_contiguous(sizeof(ComponentMSG<>), MPI_BYTE, &ComponentMSGType);
#ifdef MORE_INFORMATION
		double start, end;
#endif

		std::vector<NodeIDType> roots_of_local_trees;
//		ConnectionComponents connection_components;
		BorderComponents border_components;
		get_border_components(local_edges, cross_edges, border_components, roots_of_local_trees);
//		get_structure_of_parallel_trees(local_edges, cross_edges, connection_components, roots_of_local_trees);

		std::vector<NodeIDType> roots_of_parallel_trees;
		decide_on_global_roots(g, border_components);
#ifdef MORE_INFORMATION
		start = MPI_Wtime();
#endif

		build_forest(g, local_edges, cross_edges,
					roots_of_local_trees, roots_of_parallel_trees,
					border_components);

//		for(int p=0; p<procs; p++)
//		{
//			if(proc_id==p)
//			{
//				std::cout << "proc " << p << ":" << std::endl;
//				for(unsigned int i=0; i<roots_of_parallel_trees.size(); i++)
//				{
//					print_local_tree(g, roots_of_parallel_trees[i], forest);
//				}
//			}
//			MPI_Barrier(MPI_COMM_WORLD);
//			sleep(1);
//		}


#ifdef MORE_INFORMATION
		end = MPI_Wtime();

		outfile_std_out << "proc " << proc_id << ": time to run build_forest=" << (end-start) << " sec" << std::endl;

		get_parallel_tree_information(roots_of_parallel_trees);
#endif

		// set the information about the local trees
		roots_of_local_trees_ptr = &roots_of_local_trees;
		next_local_tree_pos = 0;

		get_parallel_tree_matching(roots_of_local_trees,
									roots_of_parallel_trees,
									result,
									matched_nodes);

		reset(local_edges, cross_edges, roots_of_parallel_trees);
	}


}; // end of class ParallelTreeMatching




}// end of namespace dipl


template <typename T>
std::ostream& operator<< (std::ostream &out, const typename dipl::ParallelForestMatching<T>::ConnectionComponents &c)
{
	for(typename dipl::ParallelForestMatching<T>::ConnectionComponents::const_iterator c_it=c.begin(); c_it!=c.end(); c_it++)
	{
		out << *c_it << std::endl;
	}

	return out;
}


template <typename T>
std::ostream& operator<< (std::ostream &out, const typename dipl::ParallelForestMatching<T>::ConnectionComponent &c)
{
	unsigned int size = c.size();
	if(size==0) return out;

	typename dipl::ParallelForestMatching<T>::Partner p;

	unsigned int pos=0;

	out << c.root_candidate << "(";
	if(c.root_edge_ghost_node == -1 && c.root_edge_local_node==-1)
	{
		out << "NULL): { ";
	}
	else
	{
		out << c.root_edge_ghost_node << ", " << c.root_edge_local_node << "): { ";
	}
	for(typename dipl::ParallelForestMatching<T>::ConnectionComponent::const_iterator c_it=c.begin(); pos<size-1; c_it++, pos++)
	{
		p = *c_it;
		out << p << ", ";
	}

	p = c.back();
	out << p << " }";

	return out;
}




#endif
