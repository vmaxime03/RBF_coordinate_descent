#ifndef SDF_HPP__
#define SDF_HPP__

#include "rbf.hpp"
#include "ultimaille/algebra/vec.h"
#include <cstddef>
#include <sstream>
#include <vector>

#include "nlohmann/json.hpp"
using json = nlohmann::json;


struct SDF {
	std::vector<UM::vec2> p;
	std::vector<double> alpha;
	std::vector<UM::vec2> beta;
	std::vector<double> sigma;

	std::unique_ptr<RBF> rbf;

	SDF(std::unique_ptr<RBF> _rbf) : rbf(std::move(_rbf)) {};

	inline double contrib_dist_i(size_t i, UM::vec2 pos) {
		double n = (pos - p[i]).norm();
		if (n == 0) return 0.;
		return alpha[i] * rbf->f(n, sigma[i]) 
				+ rbf->df(n, sigma[i]) *  (beta[i] * ((pos - p[i])/n));

	}

	double distance(UM::vec2 pos) {
		double t = 0.;

		for (size_t i = 0; i < p.size(); ++i) {
			t += contrib_dist_i(i, pos);
				}
		return t;
	}

	// https://mobile.rodolphe-vaillant.fr/images/pdfs/hrbf.pdf#section.3
	
	inline UM::vec2 contrib_grad_i(size_t i, UM::vec2 pos) {
		auto diff = pos - p[i];
		double l = diff.norm();

		if (l <= 0.0000001f) return {0., 0.};

		auto diff_normalized = diff.normalized();

		double dphi = rbf->df(l, sigma[i]);
		double ddphi = rbf->ddf(l, sigma[i]);

		double alpha_dphi = alpha[i] * dphi;
		double beta_dot_diff_l = (beta[i] * diff)/l;
		double squared_l = diff.norm2();

		return alpha_dphi * diff_normalized
			+	beta_dot_diff_l * (ddphi * diff_normalized - diff * dphi / squared_l)
			+	beta[i] * dphi / l;


	}
	
	UM::vec2 gradient(UM::vec2 pos) {
		UM::vec2 grad(0., 0.);

		for (size_t i = 0; i < p.size(); ++i) {
			grad += contrib_grad_i(i, pos);	
		
		}
		return grad;
	}
	




	std::string to_string() {
		std::stringstream ss;
		ss << "SDF(\n";
		for (size_t i = 0; i < p.size(); ++i) {
			ss << "\tpi=(" << p[i].x << ", " << p[i].y << "), ai="	<< alpha[i] << ", bi=(" << beta[i].x << ", " << beta[i].y << "), s=" << sigma[i] << ";\n";
		}

		return ss.str();
	}

	json to_json() {
		json j;
		j["rbf"] = rbf.name();
		j["param"] = sigma[0]; // TODO ? 

		j["points"] = json::array();

		for (size_t i = 0; i < p.size(); ++i) {
			j["points"].push_back({
					{"px", p[i].x},        
					{"py", p[i].y},
					{"alpha", alpha[i]},
					{"bx", beta[i].x},
					{"by", beta[i].y},
					{"s", sigma[i]}
					});
		}
		return j;
	}
};

#endif
