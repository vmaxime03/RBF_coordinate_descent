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
#include "ultimaille/sparse/linexpr.h"

using namespace UM;
using namespace UM::Linear;


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

double error_on_polyline(SDF& sdf, PointNormal& p) {
	double d = sdf.distance(p.first); // distance 
	vec2 g = sdf.gradient(p.first); // gradiant
	
	vec2 diff = (g - p.second);

	// pi on polyline 
	// distance -> 0 
	// gradient -> ni
	return d * d + diff.norm2();
}


// TODO 
double error_off_polyline(SDF& sdf, size_t npoint = 8) {
	const double Rp = 1.5;
	const double anglep = 2 * std::numbers::pi / npoint;

	double t = 0.;
	
	for (auto i = 0; i < npoint; ++i) {
		auto p = UM::vec2(Rp * cos(i * anglep), Rp * sin(i * anglep));
	
		double e = (1 / sdf.distance(p));
		t += e;
	}
	return t;
}



// total error score
double error_total(SDF& sdf, Samples& samples) {
	double t = 0.;
	for (size_t i = 0; i < samples.size(); ++i) {
		t += error_on_polyline(sdf, samples[i]);	
	}

	// TODO
	// t += error_off_polyline(sdf);
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

	PolyLineGenerator::regular_polygon(pl, sdf, 5, 5);
	
	//PolyLineGenerator::random_polygon(pl, sdf, 20, 10);
	
	//PolyLineGenerator::read_from_file(pl, sdf, input_dir + "duck.geogram", 20);


	// 	   lower bound		upper bound 	step
	// double lb_alpha =  -5.,  ub_alpha = 10., 	step_alpha = 0.01;
	// double lb_beta  = -5.,  ub_beta  = 5.,  	step_beta  = 0.01;
	
	double lb_point = -5.,  ub_point = 5.,  	step_point = 0.01;
	double lb_sigma = 0.5,	ub_sigma = 10.,  	step_sigma = 0.01;


	size_t max_it = 100;
	size_t max_it_descent = 20;
	size_t snapshot = 20;


	Samples samples;
	compute_samples_normals(pl, samples, 10);

	bool c = true;
	for (size_t it = 0; it < max_it; ++it) {

		// OUTPUT	
		if (it % (max_it / snapshot) == 0 || !c) {
			double total_err = error_total(sdf, samples);
			std::cout << it << ": err : " << total_err << "\t"  << sdf.to_string() << std::endl;
			std::ofstream(output_dir + std::to_string(it) + "sdf.json") << sdf.to_json().dump(2);
		}

		// convergence check
		if (!c) break;
		c = false;


		// find alphas & betas
		LeastSquares ls(3*sdf.p.size());


		for (const auto& s : samples) {

			LinExpr value_res;
			LinExpr grad_res[2];

			for (size_t i = 0; i < sdf.p.size(); ++i) { 

				vec2 p = s.first - sdf.p[i]; // x - p_i 
				double l = p.norm();


				double phi 		= sdf.rbf->f  (l, sdf.sigma[i]);
				double dphi 	= sdf.rbf->df (l, sdf.sigma[i]);
				double ddphi	= sdf.rbf->ddf(l, sdf.sigma[i]);

				// f(p_i) == 0
				LinExpr term = X(i*3) * phi +	X(i*3+1) * (dphi * p.x/p.norm()) + X(i*3+2) * (dphi * p.y/p.norm());
				
				value_res += term;

				// grad f(p_i) == n_i
				auto dotprod = [&](size_t i) { return (X(i*3 + 1) * p[0] + X(i*3 + 2) * p[1]);};

				for (size_t k = 0; k < 2; ++k) {

					LinExpr alpha_term = X(i*3) * dphi * ( p[k] / l );

					LinExpr beta1_term = ddphi * ( p[k] / (l * l) ) * dotprod(i);

					LinExpr beta2_term = (dphi / l) * X(i*3 + 1 + k) ;

					LinExpr beta3_term = ((dphi * p[k])/(l*l*l)) * dotprod(i);

					grad_res[k] += ((alpha_term + beta1_term + beta2_term - beta3_term));
				}


			}

			vec2 n = s.second; // n_i

			ls.add_to_energy(value_res);
			ls.add_to_energy(grad_res[0] - n[0]);
			ls.add_to_energy(grad_res[1] - n[1]);

		}



		ls.solve();


		for (size_t i = 0; i < sdf.p.size(); ++i) { 

			sdf.alpha[i] = ls.value(i*3);
			sdf.beta[i].x = ls.value(i*3+1);
			sdf.beta[i].y = ls.value(i*3+2);

		}


		for (size_t t = 0; t < max_it_descent; ++t) {		
			c = false;
			for (size_t i = 0; i < sdf.p.size(); ++i) {

				double ce = error_total(sdf, samples);

				//c = true;

				//c |= try_step(sdf.alpha[i],  step_alpha, lb_alpha, ub_alpha, sdf, i, samples, ce);
				//c |= try_step(sdf.beta[i].x, step_beta,  lb_beta,  ub_beta,  sdf, i, samples, ce);
				//c |= try_step(sdf.beta[i].y, step_beta,  lb_beta,  ub_beta,  sdf, i, samples, ce);

				c |= try_step(sdf.p[i].x,    step_point, lb_point, ub_point, sdf, i, samples, ce);
				c |= try_step(sdf.p[i].y,    step_point, lb_point, ub_point, sdf, i, samples, ce);
				c |= try_step(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, sdf, i, samples, ce);


			}
			if (!c) break;
		}
	}


	// OUTPUT
	export_polyline(pl, output_dir + "polyline.csv");
	std::ofstream(output_dir + "sdf.json") << sdf.to_json().dump(2);
	debug_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");

	return 0;
}
