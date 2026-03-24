
#include <common/read_matrix_market.h>

template <typename NodeIDType, typename EdgeIDType, typename WeightType>
void write_graph(const std::string &src_graph,
				 const std::string &filename, const EdgeIDType num_edges,
				 const WeightType weight_adjustment,
				 const dipl::AdjacencyList<NodeIDType, WeightType> &graph)
{
	unsigned int one_per_cent = graph.num_vertices()/100;
	one_per_cent += (one_per_cent==0);
	unsigned int percent = 0;

	std::ofstream outfile(filename.c_str());
	if(outfile.good())
	{
		// write some information about the original graph

		outfile << "% This file was created from the graph: " << src_graph << std::endl;
		outfile << "% This graph represents the bipartite version of this source matrix" << std::endl;
		outfile << "% Row and Column Nodes are interleaved." << std::endl;
		outfile << "% The edge weights are adjusted if negative edge weights exists," << std::endl;
		outfile << "% Such that the smallest edge weight is 0." << std::endl;

		// write header
		outfile << graph.num_vertices() << " " << num_edges << " 10" << std::endl;

		// write nodes
		for(NodeIDType i=0; i<graph.num_vertices(); i++)
		{
			if(i%one_per_cent == 0)
			{
				std::cout << percent++ << "%  ";
				std::cout.flush();
			}
			NodeIDType current_node = i;

			for(typename dipl::AdjacencyList<NodeIDType, WeightType>::incident_edge_iterator e_it=graph.begin_incident_edges(current_node);
																					e_it!=graph.end_incident_edges(current_node);
																					e_it++)
			{
				// if we have negative edge weights adjust them, such that they start with weight 0
				outfile << " " << (e_it->first+1) << " " << (e_it->second-weight_adjustment);
			}
			outfile << std::endl;
		}

		std::cout << std::endl;
	}
}


int main(int argc, char **argv)
{
	std::string infile_graph, outfile_graph;

	if(argc!=3)
	{
		std::cout << "infile_graph outfile_graph" << std::endl;
		return 1;
	}

	infile_graph = argv[1];
	outfile_graph = argv[2];


	dipl::AdjacencyList<unsigned int, double> g;

	unsigned int num_edges;
	double weight_adjustment;

	std::cout << "Reading graph " << infile_graph << std::endl;

	dipl::get_matrix_market_graph(infile_graph, g, num_edges, weight_adjustment);

	std::cout << "Read graph " << infile_graph << std::endl;

	write_graph(infile_graph, outfile_graph, num_edges, weight_adjustment, g);

	return 0;
}


