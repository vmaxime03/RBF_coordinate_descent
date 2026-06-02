#ifndef ALGO_IMP3__
#define ALGO_IMP3__ 

#include "algo_base.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/knn.h"
#include <algorithm>
#include <cstddef>
#include <numbers>
#include <ostream>
#include "output.hpp"


namespace sdffitting {


struct DecimateFitter : Fitter {
	using Fitter::Fitter;
	
	double alpha_remove_threshold = 1e-5;
	double beta_remove_threshold = 1e-5;

	void decimate(double err) {

		/* TODO optimise : only update neighbors
		std::vector<vec2> samples_points;
		samples_points.reserve(samples.size());
		std::transform(samples.begin(), samples.end(), samples_points.begin(), [](const auto& a) { return a.point; });

		auto knn = UM::KNN<2>(samples_points);
		*/

		for (size_t i = sdf.fonctions.size(); i-- > 0;) {

			// test if candidate for removal
			if (sdf.fonctions[i].alpha < alpha_remove_threshold 
				|| sdf.fonctions[i].beta.norm() < beta_remove_threshold 
				) {
				
				sdf.active[i] = false;

				resolve_ls();

				double nerr = error_total();

				if (err < nerr) {
					sdf.active[i] = true;
				} else {
					sdf.delete_point(i);
				}
			}
		}
	}

	void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {

		run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool terminated = false;

				auto s = sdf.fonctions.size();

				resolve_ls();
				double err = error_total();
				decimate(err);

				terminated = s = sdf.fonctions.size();

				return terminated;

				});
	}


};

}

#endif
