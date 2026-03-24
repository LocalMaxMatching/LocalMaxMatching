#ifndef INPUT_INFO_H_
#define INPUT_INFO_H_

#include <sstream>

namespace dipl
{

template <typename NodeIDType>
std::string input_info(const bool read_graph, const std::string &path,
						const bool create_kn, const NodeIDType num_vertices,
						const bool create_grid, const unsigned int dim, const unsigned int dim_length,
						const std::string &rating)
{
	std::string result = "";

	if(read_graph)
	{
		result = path + "::" + rating;
	}
	else if(create_kn)
	{
		std::stringstream stream_num_vertices;
		stream_num_vertices << num_vertices;
		result = "Kn_" + stream_num_vertices.str() + "::" + rating;
 	}
	else if(create_grid)
	{
		std::stringstream stream_dim, stream_dim_length;
		stream_dim << dim;
		stream_dim_length << dim_length;
		result = "grid_" + stream_dim.str() + "_" + stream_dim_length.str() + "::" + rating;
	}

	return result;
}

}

#endif
