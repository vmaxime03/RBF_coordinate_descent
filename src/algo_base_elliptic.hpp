#ifndef ALGO_BASE_ELLIPSE_HPP__
#define ALGO_BASE_ELLIPSE_HPP__

#include "debug_macros.hpp"
#include "sdf.hpp"
#include "samples.hpp"
#include "ultimaille/algebra/mat.h"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include "ultimaille/helpers/knn.h"
#include "ultimaille/sparse/least_squares.h"
#include "ultimaille/sparse/linexpr.h"
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <csignal>
#include <cstddef>
#include <fstream>
#include <vector>
#include "sdf_elliptic.hpp"

namespace sdffitting_elliptic {

	using namespace UM;
	using namespace UM::Linear;

		struct Fitter_Elliptic {

			SDF_Elliptic&        sdf;
			samples::Samples samples;

			bool   fix_alpha_zero  = false;
			bool   fix_beta_zero   = false;


			// error coefs 
			double lambda_distance = 1.;
			double lambda_gradient = 1.;

		explicit Fitter_Elliptic(SDF_Elliptic& _sdf, samples::Samples _samples)
        : sdf(_sdf), samples(std::move(_samples)) {}



		virtual ~Fitter_Elliptic() = default;

		// abstract
		virtual void fit(size_t max_it, size_t snapshot, const std::string& output_dir) = 0;


		void resolve_ls() {


			TIME_DEBUG_INIT();
			std::vector<std::tuple<UM::mat2x2, UM::mat2x2>> cache(sdf.fonctions.size());

#pragma omp parallel for
			for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
				auto A = sdf.fonctions[i].ellipse_major;
				auto B = sdf.fonctions[i].ellipse_minor;

				auto a2 = A.norm2();
				auto b2 = B.norm2();

				UM::mat2x2 M = {{{A.x / a2, A.y / a2}, {B.x / b2, B.y / b2 }}};
				auto Mt = M.transpose();

				cache[i] = {M, Mt};
			}

			TIME_DEBUG("cache init");

			LeastSquares ls(3 * sdf.fonctions.size());

			if (fix_alpha_zero) {
				for (size_t i = 0; i < sdf.fonctions.size(); ++i)
					ls.fix(i * 3, 0.);
			}
			if (fix_beta_zero) {
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					ls.fix(i * 3 + 1, 0.);
					ls.fix(i * 3 + 2, 0.);
				}
			}


			DEBUG("LS init");
			TIME_DEBUG_INIT();

#pragma omp parallel for schedule(dynamic)
			for (const auto& s : samples) {
				LinExpr value_res;
				LinExpr grad_res[2];


				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					if (!sdf.active[i]) continue;

					// VARIABLES 
					auto alpha = X(i*3);
					auto betax = X(i*3 + 1);

					auto betay = X(i*3 + 2);

					// VALUE TERM 
					auto d = s.point - sdf.fonctions[i].point;

					auto [M, Mt] = cache[i];

					auto y = M * d;

					double e = y.norm();

					if (e >= 1) continue; // phi will be 0 anyway
										  // case when e->0 (center is on sample)
					if (e <= 1e-14) { // TODO this calculus assume WendlandC2 rbf, fix in refac
						value_res += alpha;
						auto MtM = Mt * M;
						for (size_t k = 0; k < 2; ++k) {
							auto MtMk = MtM.rows[k];
							grad_res[k] += -20 * (MtMk.x * betax + MtMk.y * betay);
						}
						continue;
					} 

						

					auto grad_e = (Mt * y) / y.norm();

					auto phi = sdf.rbf->f(e, 1.);
					auto dphi = sdf.rbf->df(e, 1.);
					auto ddphi = sdf.rbf->ddf(e, 1.);

					LinExpr term = alpha * phi + dphi * (betax * grad_e.x + betay * grad_e.y);

					value_res += term;


					// GRAD TERM

					UM::mat<2, 1> y_as_mat = {{{y.x}, {y.y}}};

					auto yyt = y_as_mat * y_as_mat.transpose();

					auto hess_e = (1/e) * Mt * M - 1/(e*e*e) * Mt * yyt * M;


					for (size_t k = 0; k < 2; ++k) {
						LinExpr alpha_term = alpha * dphi * grad_e[k];
						LinExpr beta1_term = dphi * (hess_e[k][0] * betax + hess_e[k][1] * betay);
						LinExpr beta2_term = ddphi * ( grad_e[k] * (betax * grad_e.x + betay * grad_e.y) );
						grad_res[k] += (alpha_term + beta1_term + beta2_term);
					}
				}


				vec2 n = s.normal;

				LinExpr equality_x = (grad_res[0] - n[0]);
				LinExpr equality_y = (grad_res[1] - n[1]);

#pragma omp critical
				{
					ls.add_to_energy(lambda_distance * value_res);
					ls.add_to_energy(equality_x);
					ls.add_to_energy(equality_y);
				}
			}

			TIME_DEBUG("LS init terminated");

			ls.solve();


			TIME_DEBUG("LS solve");


			for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
				sdf.fonctions[i].alpha   = ls.value(i*3);
				sdf.fonctions[i].beta.x  = ls.value(i*3+1);
				sdf.fonctions[i].beta.y  = ls.value(i*3+2);
			}
		}




		double error_on_sample(samples::PointNormal& p) const {
			double d   = sdf.distance(p.point);
			vec2   g   = sdf.gradient(p.point);
			vec2   diff = g - p.normal;
			return lambda_distance * d*d + lambda_gradient * diff.norm2();
		}

		double error_total(std::vector<size_t>* sorted_idx = nullptr) {
			double total = 0.;
			std::vector<std::pair<size_t, double>> errors;
			if (sorted_idx) errors.reserve(samples.size());

			for (size_t i = 0; i < samples.size(); ++i) {
				double err = error_on_sample(samples[i]);
				if (sorted_idx) errors.emplace_back(i, err);
				total += err;
			}

			if (sorted_idx) {
				std::sort(errors.begin(), errors.end(),
					[](const auto& a, const auto& b){ return a.second > b.second; });
				sorted_idx->resize(errors.size());
				std::transform(errors.begin(), errors.end(), sorted_idx->begin(),
						[](const auto& p){ return p.first; });
			}
			return total;
		}

	};
//
// 	struct TestEllipse : Fitter_Elliptic {
// 		using Fitter_Elliptic::Fitter_Elliptic;
//
// 		void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {
//
// 		}
//
// };


} // namespace sdffitting
  


#endif
