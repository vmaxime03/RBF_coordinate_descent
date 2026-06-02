#ifndef SDF_ELLIPTIC_HPP__
#define SDF_ELLIPTIC_HPP__


#include "ultimaille/algebra/mat.h"
#include "ultimaille/algebra/vec.h"
#include "rbf.hpp"
#include <memory>
#include <sstream>
#include <vector>
#include "sdf_base.hpp"

struct FonctionElliptic {
	UM::vec2 point;
	double alpha;
	UM::vec2 beta;
	UM::vec2 ellipse_minor;
	UM::vec2 ellipse_major;
};

struct SDF_Elliptic : SDF_Base<FonctionElliptic> {
	using SDF_Base::SDF_Base;

	inline double fi(size_t i, UM::vec2 pos) const override {
		auto d = pos - fonctions[i].point;

		auto A = fonctions[i].ellipse_major;
		auto B = fonctions[i].ellipse_minor;

		auto a2 = A.norm2();
		auto b2 = B.norm2();

		if (a2 <= 1e-14 || b2 <= 1e-14) return 0.;

		UM::mat2x2 M = {{{A.x / a2, A.y / a2}, {B.x / b2, B.y / b2 }}};

		auto y = M * d;

		auto Mt = M.transpose();

		double e = y.norm();

		auto grad_e = (Mt * y) / y.norm();

		return fonctions[i].alpha * rbf->f(e, 1.) + rbf->df(e, 1.) * (fonctions[i].beta * grad_e);
	}

	inline UM::vec2 gi(size_t i, UM::vec2 pos) const override {
		auto d = pos - fonctions[i].point;

		auto A = fonctions[i].ellipse_major;
		auto B = fonctions[i].ellipse_minor;

		auto a2 = A.norm2();
		auto b2 = B.norm2();

		if (a2 <= 1e-14 || b2 <= 1e-14) return {0., 0.};

		UM::mat2x2 M = {{{A.x / a2, A.y / a2}, {B.x / b2, B.y / b2 }}};

		auto y = M * d;

		auto Mt = M.transpose();

		double e = y.norm();

		auto grad_e = (Mt * y) / e;

		UM::mat<2, 1> y_as_mat = {{{y.x}, {y.y}}};

		auto yyt = y_as_mat * y_as_mat.transpose();

		auto hess_e = (1/e) * Mt * M - 1/(e*e*e) * Mt * yyt * M;

		return fonctions[i].alpha * rbf->df(e, 1.) * grad_e
			+ rbf->df(e, 1.) * hess_e * fonctions[i].beta
			+ rbf->ddf(e, 1.) * grad_e * (fonctions[i].beta * grad_e);

	}



	std::string to_string() const override {
		std::stringstream ss;
		ss << "SDF:\n";
		for (size_t i = 0; i < fonctions.size(); ++i) {
			ss 	<<	"\tpi=(" << fonctions[i].point.x << ",\t" << fonctions[i].point.y 
				<< "),\tai="	<< fonctions[i].alpha 
				<< ",\tbi=(" << fonctions[i].beta.x << ",\t" << fonctions[i].beta.y 
				<< "),\tA=" << fonctions[i].ellipse_major.x << ",\t" << fonctions[i].ellipse_major.y
				<< "),\tB=" << fonctions[i].ellipse_minor.x << ",\t" << fonctions[i].ellipse_minor.y
				<< ";\n";
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
					{"s", {
						{"Ax", fonctions[i].ellipse_major.x},
						{"Ay", fonctions[i].ellipse_major.y},
						{"Bx", fonctions[i].ellipse_minor.x},
						{"By", fonctions[i].ellipse_minor.y}
						}}
					});
		}
		return j;
	}
	

};






#endif // !SDF_NEW_HPP__
