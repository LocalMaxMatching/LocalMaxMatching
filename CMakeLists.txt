cmake_minimum_required(VERSION 3.12)
project(LocalMaxMatching LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/deploy)

# --- Dependencies ---
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(MPI QUIET)

# All targets include the project root for header resolution
include_directories(${CMAKE_SOURCE_DIR})

# --- Sequential matching programs (Boost required) ---
set(SEQ_BOOST_TARGETS
    local_max_matching
    local_tree_matching
    karp_sipser_matching
)

foreach(target ${SEQ_BOOST_TARGETS})
    add_executable(${target} programs/${target}.cpp)
    target_link_libraries(${target} PRIVATE Boost::program_options)
endforeach()

add_executable(mixed-karp_sipser-local_max programs/mixed-karp_sipser-local_max.cpp)
target_link_libraries(mixed-karp_sipser-local_max PRIVATE Boost::program_options)

# --- Sequential utility programs (no Boost) ---
add_executable(add_random_weights_to_metis_graph programs/add_random_weights_to_metis_graph.cpp)
add_executable(create_node_index_for_metis_graph programs/create_node_index_for_metis_graph.cpp)
add_executable(arrange_nodes_along_a_curve programs/arrange_nodes_along_a_curve.cpp)
add_executable(matrix_market_to_bipartite_metis programs/matrix_market_to_bipartite_metis.cpp)
add_executable(compute_edge_distribution programs/compute_edge_distribution.cpp)

# --- MPI programs ---
if(MPI_FOUND)
    message(STATUS "MPI found -- building parallel targets")

    add_executable(parallel_local_max_matching programs/parallel_local_max_matching.cpp)
    target_compile_definitions(parallel_local_max_matching PRIVATE _MPI)
    target_link_libraries(parallel_local_max_matching PRIVATE MPI::MPI_CXX Boost::program_options)

    add_executable(parallel_local_max_matching_more_info programs/parallel_local_max_matching.cpp)
    target_compile_definitions(parallel_local_max_matching_more_info PRIVATE _MPI MORE_INFORMATION)
    target_link_libraries(parallel_local_max_matching_more_info PRIVATE MPI::MPI_CXX Boost::program_options)
else()
    message(STATUS "MPI not found -- skipping parallel targets")
endif()

# --- Optional: LEMON-based optimal matching ---
find_library(LEMON_LIB lemon)
find_path(LEMON_INC lemon/list_graph.h)

if(LEMON_LIB AND LEMON_INC)
    message(STATUS "LEMON found -- building optimal_lemon_matching")
    add_executable(optimal_lemon_matching optimal_matchings/optimal_lemon_matching.cpp)
    target_include_directories(optimal_lemon_matching PRIVATE ${LEMON_INC})
    target_link_libraries(optimal_lemon_matching PRIVATE ${LEMON_LIB} Boost::program_options)
else()
    message(STATUS "LEMON not found -- skipping optimal_lemon_matching")
endif()
