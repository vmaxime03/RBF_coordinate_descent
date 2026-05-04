#include "ultimaille/polyline.h"
#include <fstream>
#include <iostream>
#include <ostream>
#include <ultimaille/all.h>
#include <vector>


#include "geo_polyline.hpp"
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"

using namespace UM;

template <typename RBF_T>
void debug_sdf(SDF<RBF_T>& sdf, double minx, double miny, double maxx, double maxy, const std::string& fname, int step = 100) {
	std::ofstream out(fname);
	int w = step;
	int h = step;

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {

            const double x = minx + (double(i) / double(w - 1)) * (maxx - minx);
            const double y = miny + (double(j) / double(h - 1)) * (maxy - miny);

			double dist = sdf.distance({x, y});

			out << x << "," << y << "," << dist;
			if (j < h -1) out << ",";
		}
		out << "\n";
	}
	out.close();

}




template <typename RBF_T>
double err(PolyLine& pl, SDF<RBF_T>& sdf, double ntest = 10) {
	double t = 0.;

	for (const auto& e : pl.iter_edges()) {
		vec2 d = e.to().pos().xy() - e.from().pos().xy();

		for (size_t i = 0; i < ntest; ++i) {
			vec2 p = e.from().pos().xy() + ((double)i/(double)ntest) * d;

			double diff = signed_distance_point_polyline_2d(pl, p) - sdf.distance(p);
			t += diff * diff;
		}
	}
	return t;
}



int main(int argc, char** argv) {

	
	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::create_directories(output_dir);


	PolyLine pl;
	SDF<Gaussian> sdf;
	gen_1(pl, sdf);


	double lb_point = -5.;
	double ub_point = 5.;
	double step_point = 0.001;

	double lb_alpha = 0.; 
	double ub_alpha = 10.;
	double step_alpha = 0.001;
	
	double lb_beta = -5.;
	double ub_beta = 5.;
	double step_beta = 0.001;

	double lb_sigma = 0.01;
	double ub_sigma = 5.;
	double step_sigma = 0.001;

	auto try_step = [&](double& var, double step, double lb, double ub) {
		double curr_err = err(pl, sdf);
		double old = var;

		// cas + step meilleur
		var = std::clamp(old + step, lb, ub);
		if (err(pl, sdf) < curr_err) return;

		// cas - step meilleur
		var = std::clamp(old - step, lb, ub);
		if (err(pl, sdf) < curr_err) return;

		// aucun des deux
		var = old;
	};
	
	size_t max_it = 2500; // TODO use argv
	size_t it;
	for (it = 0; it < max_it; ++it) {
		
		double old_err = err(pl, sdf);

		if (it % (max_it/100) == 0) {
			std::cout << it << ": err : " << old_err << "	"
				<< sdf.to_string()
				<< std::endl;

			std::ofstream(output_dir + std::to_string(it) + "sdf.json") << sdf.to_json().dump(2);
		}

		for (size_t i = 0; i < sdf.p.size(); ++i) {
			try_step(sdf.alpha[i], step_alpha, lb_alpha, ub_alpha);
			try_step(sdf.beta[i].x, step_beta, lb_beta, ub_beta);
			try_step(sdf.beta[i].y, step_beta, lb_beta, ub_beta);
			try_step(sdf.p[i].x, step_point, lb_point, ub_point);
			try_step(sdf.p[i].y, step_point, lb_point, ub_point);
			try_step(sdf.sigma[i], step_sigma, lb_sigma, ub_sigma);

		}


		

	}

	std::ofstream(output_dir + "sdf.json") << sdf.to_json().dump(2);
	debug_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");

	
	




	return 0;
}
