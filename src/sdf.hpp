#ifndef SDF_HPP__
#define SDF_HPP__

#include "rbf.hpp"
#include "ultimaille/algebra/vec.h"
#include <concepts>
#include <cstddef>
#include <sstream>
#include <vector>

#include "sdf_base.hpp"

struct FonctionCirculaire {
	UM::vec2 point;
	double alpha;
	UM::vec2 beta;
	double sigma;
};

struct SDF : SDF_Base<FonctionCirculaire> {
	using SDF_Base::SDF_Base;

	inline double fi(size_t i, UM::vec2 pos) const override {
		double n = (pos - fonctions[i].point).norm();
		if (n == 0) return 0.;
		return fonctions[i].alpha * rbf->f(n, fonctions[i].sigma) 
				+ rbf->df(n, fonctions[i].sigma) *  (fonctions[i].beta * ((pos - fonctions[i].point)/n));

	}

	// https://mobile.rodolphe-vaillant.fr/images/pdfs/hrbf.pdf#section.3
	
	inline UM::vec2 gi(size_t i, UM::vec2 pos) const override {
		auto diff = pos - fonctions[i].point;
		double l = diff.norm();

		if (l <= 0.0000001f) return {0., 0.};

		auto diff_normalized = diff.normalized();

		double dphi = rbf->df(l, fonctions[i].sigma);
		double ddphi = rbf->ddf(l, fonctions[i].sigma);

		double alpha_dphi = fonctions[i].alpha * dphi;
		double beta_dot_diff_l = (fonctions[i].beta * diff)/l;
		double squared_l = diff.norm2();

		return alpha_dphi * diff_normalized
			+	beta_dot_diff_l * (ddphi * diff_normalized - diff * dphi / squared_l)
			+	fonctions[i].beta * dphi / l;


	}
	


	std::string to_string() const override {
		std::stringstream ss;
		ss << "SDF:\n";
		for (size_t i = 0; i < fonctions.size(); ++i) {
			ss << "\tpi=(" << fonctions[i].point.x << ",\t" << fonctions[i].point.y << "),\tai="	<< fonctions[i].alpha << ",\tbi=(" << fonctions[i].beta.x << ",\t" << fonctions[i].beta.y << "),\ts=" << fonctions[i].sigma << ";\n";
		}

		return ss.str();
	}

	nlohmann::json to_json() const override {
		nlohmann::json j;
		j["rbf"] = rbf->name();
		j["points"] = nlohmann::json::array();

		for (size_t i = 0; i < fonctions.size(); ++i) {
			j["points"].push_back({
					{"px", fonctions[i].point.x},        
					{"py", fonctions[i].point.y},
					{"alpha", fonctions[i].alpha},
					{"bx", fonctions[i].beta.x},
					{"by", fonctions[i].beta.y},
					{"s", fonctions[i].sigma}
					});
		}
		return j;
	}

};

#endif
