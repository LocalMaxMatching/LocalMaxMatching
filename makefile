

ifneq ($(CXX),clang++)
ifneq ($(CXX),icpc)
CXX = g++
endif
endif

ifeq ($(BOOST_INC),)
BOOST_INC=.
endif

ifeq ($(BOOST_LIB),)
BOOST_LIB=.
endif

MPI_CXX = mpicxx


INC_PATHS=-I$(PWD) -I$(BOOST_INC)
LIB_PATHS=-L$(BOOST_LIB)

#OPTIMIZATION = -O3
#OPTIMIZATION = -O2 -fno-inline
OPTIMIZATION = -g

COMPILE_OPTS = $(OPTIMIZATION) -std=c++0x -Wall $(LIB_PATHS)
#COMPILE_OPTS = -g -Wall -fno-inline -D_GLIBCXX_DEBUG $(LIB_PATHS)
MPI_COMPILE_OPTS = $(COMPILE_OPTS) -D_MPI
LIBS= -Wl,-Bstatic -lboost_program_options -Wl,-Bdynamic
#LIBS= -Wl,-static -lboost_program_options 
PROGRAM_PATH = programs
GRAPH_PATH = graphs
MATCHING_PATH = matching_algorithms
COMMON_PATH = common
BIN_PATH = bin

COMMON_FILES = $(COMMON_PATH)/edge.h $(COMMON_PATH)/hash_functions.h\
				$(COMMON_PATH)/linear_buffered_reader.h \
				$(COMMON_PATH)/metis_dimacs_reader.h\
				$(COMMON_PATH)/parallel_metis_dimacs_reader.h\
				$(COMMON_PATH)/parallel_metis_dimacs_reader_evenly_distributed.h\
				$(COMMON_PATH)/rating_functions.h \
				$(COMMON_PATH)/grid_graph.h\
				$(COMMON_PATH)/grid_graph_2d.h \
				$(COMMON_PATH)/grid_graph_general.h\
				$(COMMON_PATH)/union_find.h $(COMMON_PATH)/k_n_graph.h\
				$(COMMON_PATH)/get_graph.h \
				$(COMMON_PATH)/read_command_line_arguments.h \
				$(COMMON_PATH)/read_matrix_market.h \
				$(COMMON_PATH)/input_info.h \
				$(COMMON_PATH)/mpi_system_support.h \
				makefile


$(BIN_PATH)/test_tree_matching: $(PROGRAM_PATH)/test_tree_matching.cpp $(GRAPH_PATH)/edge_graph.h $(GRAPH_PATH)/tree.h $(MATCHING_PATH)/tree_matching.h 
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


$(BIN_PATH)/mem_usg_of_hypercube_matching: $(PROGRAM_PATH)/mem_usg_of_hypercube_matching.cpp $(GRAPH_PATH)/parallel_edge_graph.h $(COMMON_FILES)
	$(MPI_CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


$(BIN_PATH)/parallel_local_tree_matching: $(PROGRAM_PATH)/parallel_local_tree_matching.cpp $(PROGRAM_PATH)/parallel_matching.h $(GRAPH_PATH)/parallel_edge_graph.h $(GRAPH_PATH)/parallel_forest.h $(MATCHING_PATH)/parallel_local_tree_maximum_matching.h $(MATCHING_PATH)/parallel_forest_matching.h  $(COMMON_FILES)
	$(MPI_CXX) $(MPI_COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)

$(BIN_PATH)/parallel_local_tree_matching_more_info: $(PROGRAM_PATH)/parallel_local_tree_matching.cpp $(PROGRAM_PATH)/parallel_matching.h $(GRAPH_PATH)/parallel_edge_graph.h $(GRAPH_PATH)/parallel_forest.h $(MATCHING_PATH)/parallel_local_tree_maximum_matching.h $(MATCHING_PATH)/parallel_forest_matching.h  $(COMMON_FILES)
	$(MPI_CXX) $(MPI_COMPILE_OPTS) -DMORE_INFORMATION  $(INC_PATHS) $< -o $@ $(LIBS)


$(BIN_PATH)/parallel_local_max_matching: $(PROGRAM_PATH)/parallel_local_max_matching.cpp $(PROGRAM_PATH)/parallel_matching.h $(GRAPH_PATH)/parallel_edge_graph.h $(MATCHING_PATH)/parallel_local_maximum_matching.h $(COMMON_FILES)
	$(MPI_CXX) $(MPI_COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)

$(BIN_PATH)/parallel_local_max_matching_more_info: $(PROGRAM_PATH)/parallel_local_max_matching.cpp $(PROGRAM_PATH)/parallel_matching.h $(GRAPH_PATH)/parallel_edge_graph.h $(MATCHING_PATH)/parallel_local_maximum_matching.h $(COMMON_FILES)
	$(MPI_CXX) $(MPI_COMPILE_OPTS) -DMORE_INFORMATION  $(INC_PATHS) $< -o $@ $(LIBS)


$(BIN_PATH)/local_max_matching: $(PROGRAM_PATH)/local_max_matching.cpp $(PROGRAM_PATH)/matching.h $(GRAPH_PATH)/edge_graph.h $(MATCHING_PATH)/local_maximum_matching.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)

$(BIN_PATH)/local_max_matching_more_info: $(PROGRAM_PATH)/local_max_matching.cpp $(PROGRAM_PATH)/matching.h $(GRAPH_PATH)/edge_graph.h $(MATCHING_PATH)/local_maximum_matching.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) -DMORE_INFO_LOCAL_MAX $(INC_PATHS) $< -o $@ $(LIBS)
	

$(BIN_PATH)/local_tree_matching: $(PROGRAM_PATH)/local_tree_matching.cpp $(PROGRAM_PATH)/matching.h $(GRAPH_PATH)/edge_graph.h $(MATCHING_PATH)/local_tree_maximum_matching.h $(MATCHING_PATH)/tree_matching.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)

$(BIN_PATH)/local_tree_matching_more_info: $(PROGRAM_PATH)/local_tree_matching.cpp $(PROGRAM_PATH)/matching.h $(GRAPH_PATH)/edge_graph.h $(MATCHING_PATH)/local_tree_maximum_matching.h $(MATCHING_PATH)/tree_matching.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) -DMORE_INFO_LOCAL_TREE $(INC_PATHS) $< -o $@ $(LIBS)
	

$(BIN_PATH)/karp_sipser_matching: $(PROGRAM_PATH)/karp_sipser_matching.cpp $(PROGRAM_PATH)/matching.h  $(GRAPH_PATH)/adjacency_array.h $(GRAPH_PATH)/adjacency_array_graph_with_active_edges_vector.h $(MATCHING_PATH)/karp_sipser.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)


$(BIN_PATH)/mixed-karp_sipser-local_max: $(PROGRAM_PATH)/mixed-karp_sipser-local_max.cpp $(PROGRAM_PATH)/matching.h  $(GRAPH_PATH)/adjacency_array.h $(GRAPH_PATH)/adjacency_array_graph_with_active_edges_vector.h $(MATCHING_PATH)/mixed-karp_sipser_and_local_max.h $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@ $(LIBS)



$(BIN_PATH)/add_random_weights_to_metis_graph: $(PROGRAM_PATH)/add_random_weights_to_metis_graph.cpp $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


$(BIN_PATH)/create_node_index_for_metis_graph: $(PROGRAM_PATH)/create_node_index_for_metis_graph.cpp $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


$(BIN_PATH)/arrange_nodes_along_a_curve: $(PROGRAM_PATH)/arrange_nodes_along_a_curve.cpp $(COMMON_FILES) $(GRAPH_PATH)/adjacency_list.h
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


$(BIN_PATH)/matrix_market_to_bipartite_metis: $(PROGRAM_PATH)/matrix_market_to_bipartite_metis.cpp $(COMMON_FILES) $(GRAPH_PATH)/adjacency_list.h
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@

$(BIN_PATH)/compute_edge_distribution: $(PROGRAM_PATH)/compute_edge_distribution.cpp $(COMMON_FILES)
	$(CXX) $(COMPILE_OPTS) $(INC_PATHS) $< -o $@


bin/optimal_lemon_matching: optimal_matchings/optimal_lemon_matching.cpp
	$(CXX) $(OPTIMIZATION) $(INC_PATHS) -I$(HOME)/Development/lemon/include $< -o $@ -L$(HOME)/Development/lemon/lib -lemon $(LIBS)


all: test test_inactive_edges

.PHONY:
clean:
	@rm -rf bin/*
