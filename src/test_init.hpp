#ifndef TEST_POLYLINE_HPP__
#define TEST_POLYLINE_HPP__


#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include <cstddef>
#include <cstdlib>


#include "sdf.hpp"


template <typename RBF_T>
void gen_1(UM::PolyLine& pl, SDF<RBF_T>& sdf) {
	pl.points.create_points(3);
	pl.create_edges(3);

	pl.points[0] = UM::vec2(-1., 1.).xy0();
	pl.points[1] = UM::vec2(0., 0.).xy0();
	pl.points[2] = UM::vec2(1., 1.).xy0();

	pl.vert(0, 0) = 0;
	pl.vert(0, 1) = 1;
	pl.vert(1, 0) = 1;
	pl.vert(1, 1) = 2;
	pl.vert(2, 0) = 2;
	pl.vert(2, 1) = 0;

	pl.connect();

	int t = 3;

	sdf.p.resize(t);
	sdf.alpha.resize(t);
	sdf.beta.resize(t);
	sdf.sigma.resize(t);

	

	--t;
	sdf.p		[t] = {-1, 0};
	sdf.alpha	[t] = 1;
	sdf.beta	[t] = {1, 0};
	sdf.sigma	[t] = 1;

	--t;
	sdf.p		[t] = {1, 0};
	sdf.alpha	[t] = 1;
	sdf.beta	[t] = {-1, 0};
	sdf.sigma	[t] = 1;
	
	--t;
	sdf.p		[t] = {0, 1};
	sdf.alpha	[t] = 1;
	sdf.beta	[t] = {0, -1};
	sdf.sigma	[t] = 1;

	/*
	auto rd = [&]() { return ((double)rand()) / RAND_MAX; };

	for (int i = 0; i < t; ++i) {
		sdf.p		[i] = {2 - rd()*4, 2 - rd()*4};
		sdf.alpha	[i] = rd()*5;
		sdf.beta	[i] = {1-rd()*2, 1-rd()*2};
		sdf.sigma	[i] = rd()*2;
	}
	*/
}

#endif
