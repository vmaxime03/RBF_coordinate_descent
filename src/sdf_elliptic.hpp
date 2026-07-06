#ifndef SDF_ELLIPTIC_HPP__
#define SDF_ELLIPTIC_HPP__


#include "debug_macros.hpp"
#include "samples.hpp"
#include "ultimaille/all.h"
#include "ultimaille/algebra/vec.h"
#include "rbf.hpp"
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>
#include "sdf_base.hpp"
#include "ultimaille/helpers/knn.h"
#include "ultimaille/sparse/least_squares.h"
#include "ultimaille/sparse/linexpr.h"

struct FunctionElliptic : FunctionBase {
	UM::vec2 ellipse_major;
	UM::vec2 ellipse_minor;
};

struct SDF_Elliptic : SDF_Base<FunctionElliptic> {
	using SDF_Base::SDF_Base;



	bool optimized = false;
	std::vector<std::tuple<UM::mat2x2, UM::mat2x2>> cache;
	double max_ellipse_radius2 = 0.;

	std::vector<UM::vec2> knn_points;
	std::unique_ptr<UM::KNN<2>> knn;

	int neighborhood_size = 10;


	void optimize() {
		// CACHE
		cache.resize(fonctions.size());
#pragma omp parallel for reduction(max:max_ellipse_radius2)
		for (size_t i = 0; i < fonctions.size(); ++i) {
			auto A = fonctions[i].ellipse_major;
			auto B = fonctions[i].ellipse_minor;

			auto a2 = A.norm2();
			auto b2 = B.norm2();

			max_ellipse_radius2 = std::max(a2, max_ellipse_radius2);
			max_ellipse_radius2 = std::max(b2, max_ellipse_radius2);

			UM::mat2x2 M = {{{A.x / a2, A.y / a2}, {B.x / b2, B.y / b2 }}};
			auto Mt = M.transpose();

			cache[i] = {M, Mt};
		}
		DEBUG("max radius : " << std::sqrt(max_ellipse_radius2));

		init_knn();

		optimized = true;
	}

	void init_knn() {

		knn_points.resize(fonctions.size());

		std::transform(fonctions.begin(), fonctions.end(), knn_points.begin(), [](const auto& a) {
				return a.point;
				});

		knn = std::make_unique<UM::KNN<2>>(knn_points);

	}






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

		if (e <= 1e-14 || e >= 1)  return 0.;

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


		if (e <= 1e-14 || e >= 1)  return {0., 0.};

		auto grad_e = (Mt * y) / e;

		UM::mat<2, 1> y_as_mat = {{{y.x}, {y.y}}};

		auto yyt = y_as_mat * y_as_mat.transpose();

		auto hess_e = (1/e) * Mt * M - 1/(e*e*e) * Mt * yyt * M;

		return fonctions[i].alpha * rbf->df(e, 1.) * grad_e
			+ rbf->df(e, 1.) * hess_e * fonctions[i].beta
			+ rbf->ddf(e, 1.) * grad_e * (fonctions[i].beta * grad_e);

	}





	std::tuple<double, UM::vec2> evali(const size_t i, const UM::vec2 pos) const override {
		auto d = pos - fonctions[i].point;

		auto& [M, Mt] = cache[i];

		auto y = M * d;

		double e = y.norm();

		if (e >= 1)  return {0., {0., 0.}};
		if (e <= 1e-14) { // TODO this calculus assume WendlandC2 rbf, fix in refac
			return {
				fonctions[i].alpha,
				-20.0 * (Mt * M * fonctions[i].beta)
			};
		}

		auto grad_e = (Mt * y) / e;

		UM::mat<2, 1> y_as_mat = {{{y.x}, {y.y}}};

		auto yyt = y_as_mat * y_as_mat.transpose();

		auto hess_e = (1/e) * Mt * M - 1/(e*e*e) * Mt * yyt * M;

		return { 
			fonctions[i].alpha * rbf->f(e, 1.) + rbf->df(e, 1.) * (fonctions[i].beta * grad_e),
			fonctions[i].alpha * rbf->df(e, 1.) * grad_e
				+ rbf->df(e, 1.) * hess_e * fonctions[i].beta
				+ rbf->ddf(e, 1.) * grad_e * (fonctions[i].beta * grad_e)
		};
	}



	std::tuple<double, UM::vec2> eval(const UM::vec2 pos) override {
		if (!optimized) optimize();

		int expand_counter = 0;


		int current_n_size = neighborhood_size;

		auto neighbors = knn->query(pos, current_n_size);

		while (neighbors.size() < fonctions.size()) {

			int last  = neighbors.back();
			if ((pos - fonctions[last].point).norm2() < max_ellipse_radius2) {
				current_n_size *= 2;
				neighbors = knn->query(pos, current_n_size);
				expand_counter++;
			} else {
				break;
			}
		}

		double d = 0.;
		double gx = 0.;
		double gy = 0.;

#pragma omp parallel for reduction(+:d,gx,gy) schedule(guided)
		// for (size_t i = 0; i < fonctions.size(); ++i) {
		for (size_t ni = 0; ni < neighbors.size(); ++ni) {
			size_t i = neighbors[ni];
			if (!active[i]) continue;
			auto [di, gi] =  evali(i, pos);
			d += di;
			gx += gi.x;
			gy += gi.y;

		}
		return {d, {gx, gy}};
	}



















	std::string to_string() const override {
		std::stringstream ss;
		ss << "SDF:\n";
		for (size_t i = 0; i < fonctions.size(); ++i) {
			ss 	<<	"\tp" << i << "=(" << fonctions[i].point.x << ",\t" << fonctions[i].point.y 
				<< "),\ta" << i << "="	<< fonctions[i].alpha 
				<< ",\tb" << i << "=(" << fonctions[i].beta.x << ",\t" << fonctions[i].beta.y 
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

struct UDF_Elliptic : SDF_Elliptic {
	using SDF_Elliptic::SDF_Elliptic;


	double distance(UM::vec2 pos) override {
		double f = SDF_Elliptic::distance(pos);
		return f * f;
	}

	UM::vec2 gradient(UM::vec2 pos) override {
		auto [f, g] = SDF_Elliptic::eval(pos);
		return 2.0 * f * g; 
	}


	std::tuple<double, UM::vec2> eval(const UM::vec2 pos) override {
		auto [d, g] = SDF_Elliptic::eval(pos);
		return {d * d, 2. * d * g};
	}

};








#endif // !SDF_NEW_HPP__
