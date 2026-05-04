#ifndef SDF_HPP__
#define SDF_HPP__

#include "ultimaille/algebra/vec.h"
#include <cstddef>
#include <sstream>
#include <vector>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

using namespace UM; 

template <typename RBF_t>
struct SDF {
	std::vector<vec2> p;
	std::vector<double> alpha;
	std::vector<vec2> beta;
	std::vector<double> sigma;

	RBF_t rbf = RBF_t();

	double distance(vec2 pos) {
		double t = 0.;

		for (size_t i = 0; i < p.size(); ++i) {
			double n = (pos - p[i]).norm();
			t += alpha[i] * rbf.f(n, sigma[i]) 
				 + rbf.df(n, sigma[i]) *  (beta[i] * ((pos - p[i])/n));
		}
		return t;
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
