#ifndef RATING_FUNCTIONS_H_
#define RATING_FUNCTIONS_H_

#include <common/math_funcs.h>

namespace dipl
{


template <typename WeightType>
inline WeightType const_1_weight(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight)
{
	return 1;
}


template <typename WeightType>
inline WeightType identity_edge_weight(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight)
{
	return edge_weight;
}


template <typename WeightType>
inline WeightType expansionstar2(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight)
{
	return (edge_weight*edge_weight)/(n1_weight*n2_weight);
}

template <typename WeightType>
inline WeightType rand_weight(const WeightType edge_weight, const WeightType n1_weight, const WeightType n2_weight)
{
	return random_weight();
}


} // end namespace dipl

#endif
