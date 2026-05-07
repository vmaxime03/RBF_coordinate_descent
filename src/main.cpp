#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <tuple>
#include <ultimaille/all.h>
#include <utility>
#include <vector>

#include "geo_polyline.hpp"
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"
#include "output.hpp"

using namespace UM;


typedef std::pair<vec2, vec2> PointNormal;
typedef std::vector<PointNormal> Samples;

// points along the polyline
void  compute_samples_normals(PolyLine& pl, Samples& samples, size_t n = 10) {
	for (const auto& e : pl.iter_edges()) {
		vec2 d = e.to().pos().xy() - e.from().pos().xy();
		
		vec2 normal = vec2(-d.y, d.x).normalized();

		for (size_t i = 0; i < n; ++i) {
			vec2 p = e.from().pos().xy() + ((double)i / (double)n) * d;
			samples.push_back({p, normal});

		}
	}
}

double error(SDF& sdf, PointNormal& p) {
	double d = sdf.distance(p.first); // distance 
	vec2 g = sdf.gradient(p.first); // gradiant
	
	vec2 diff = (g - p.second);

	// || grad f || = 1
	double eikonal = g.norm() - 1.;


	// pi on polyline 
	// distance -> 0 
	// gradient -> ni
	return d * d + diff.norm2();
}



// total error score
double error_total(SDF& sdf, Samples& samples) {
	double t = 0.;
	for (size_t i = 0; i < samples.size(); ++i) {
		t += error(sdf, samples[i]);	
	}
	return t;
}

bool try_step(
		double& var,  // coordiante to optimize
		double& step, // step for the coord
		double lb, double ub, // bound for the coord
		SDF& sdf, // function
		size_t i, // current point 
		Samples& samples, // samples
		double& curr_err // current error score
		) {

	const double old = var;
	
	// try +step
    var = std::clamp(old + step, lb, ub);
    double errp = error_total(sdf, samples);

    if (errp < curr_err) { 
        curr_err = errp;
        return true;
    }

    // try -step
    var = std::clamp(old - step, lb, ub);
    double errm = error_total(sdf, samples);

    if (errm < curr_err) { 
        curr_err = errm;
        return true;
    }

    var = old;
    return false;
}



int main(int argc, char** argv) {

	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;


	PolyLine pl;
	auto rbf = std::make_unique<Gaussian>();
	SDF sdf(std::move(rbf));

	//PolyLineGenerator::regular_polygon(pl, sdf, 8, 6);
	
	PolyLineGenerator::random_polygon(pl, sdf, 8, 5);
	
	//PolyLineGenerator::read_from_file(pl, sdf, input_dir + "duck.geogram", 50);


	// 	   lower bound		upper bound 	step
	double lb_point = -5.,  ub_point = 5.,  	step_point = 0.01;
	double lb_alpha =  -5.,  ub_alpha = 10., 	step_alpha = 0.01;
	double lb_beta  = -5.,  ub_beta  = 5.,  	step_beta  = 0.01;
	double lb_sigma = 0.5,	ub_sigma = 10.,  	step_sigma = 0.01;




	size_t max_it = (argc > 1) ? std::stoul(argv[1]) : 1000;

	Samples samples;
	compute_samples_normals(pl, samples, 10);

	bool c = true;
	for (size_t it = 0; it < max_it; ++it) {

			
		if (it % std::max<size_t>(1, max_it / 20) == 0 || !c) {
			double total_err = error_total(sdf, samples);
			std::cout << it << ": err : " << total_err << "\t"  << sdf.to_string() << std::endl;

			// std::cout << it << "\t/" << max_it << "\r";	std::cout.flush();
			std::ofstream(output_dir + std::to_string(it) + "sdf.json") << sdf.to_json().dump(2);
		}

		if (!c) break;
		c = false;

		for (size_t i = 0; i < sdf.p.size(); ++i) {


			double ce = error_total(sdf, samples);


			c |= try_step(sdf.alpha[i],  step_alpha, lb_alpha, ub_alpha, sdf, i, samples, ce);
			c |= try_step(sdf.beta[i].x, step_beta,  lb_beta,  ub_beta,  sdf, i, samples, ce);
			c |=try_step(sdf.beta[i].y, step_beta,  lb_beta,  ub_beta,  sdf, i, samples, ce);
			c |= try_step(sdf.p[i].x,    step_point, lb_point, ub_point, sdf, i, samples, ce);
			c |= try_step(sdf.p[i].y,    step_point, lb_point, ub_point, sdf, i, samples, ce);
			c |= try_step(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, sdf, i, samples, ce);

			
		}

	}

	export_polyline(pl, output_dir + "polyline.csv");
	std::ofstream(output_dir + "sdf.json") << sdf.to_json().dump(2);
	debug_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");

	return 0;
}
