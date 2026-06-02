#ifndef ALGO_BASE_HPP__
#define ALGO_BASE_HPP__

#include "debug_macros.hpp"
#include "sdf.hpp"
#include "samples.hpp"
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

namespace sdffitting {

	using namespace UM;
	using namespace UM::Linear;

	struct Fitter {

		SDF&    sdf;
		samples::Samples samples;

		double lb_sigma  = 0.,   ub_sigma  = DBL_MAX,  	step_sigma = 0.01;
		double lb_point  = -5.,  ub_point  = 5.,   		step_point = 0.01;

		double ADD_POINT_ERR_THRESHOLD = 5.;
		double default_sigma_add       = 1.;

		bool   fix_alpha_zero  = false;
		bool   fix_beta_zero   = false;
		double MIN_IMPROVEMENT = 1e-15;

		// adaptative step
		double EXPAND   = 1.2;
		double SHRINK   = 0.5;
		double MIN_STEP = 0.001;


		// error coefs 
		double lambda_distance = 1.;
		double lambda_gradient = 1.;


		explicit Fitter(SDF& _sdf, samples::Samples _samples)
			: sdf(_sdf), samples(std::move(_samples)) {}

		virtual ~Fitter() = default;

		// abstract
		virtual void fit(size_t max_it, size_t snapshot, const std::string& output_dir) = 0;


		void resolve_ls() {

			if (sdf.fonctions.empty()) return;

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

			for (const auto& s : samples) {
				LinExpr value_res;
				LinExpr grad_res[2];

				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {

					if (!sdf.active[i]) continue;

					vec2   p     = s.point - sdf.fonctions[i].point;
					double l     = p.norm();
					if (l == 0.) continue;

					double phi   = sdf.rbf->f  (l, sdf.fonctions[i].sigma);
					double dphi  = sdf.rbf->df (l, sdf.fonctions[i].sigma);
					double ddphi = sdf.rbf->ddf(l, sdf.fonctions[i].sigma);

					LinExpr term =
						X(i*3)   * phi +
						X(i*3+1) * (dphi * p.x / p.norm()) +
						X(i*3+2) * (dphi * p.y / p.norm());
					value_res += term;

					auto dotprod = [&](size_t j) {
						return (X(j*3+1) * p[0] + X(j*3+2) * p[1]);
					};

					for (size_t k = 0; k < 2; ++k) {
						LinExpr alpha_term = X(i*3)       * dphi * (p[k] / l);
						LinExpr beta1_term = ddphi        * (p[k] / (l*l)) * dotprod(i);
						LinExpr beta2_term = (dphi / l)   * X(i*3+1+k);
						LinExpr beta3_term = ((dphi*p[k]) / (l*l*l)) * dotprod(i);
						grad_res[k] += (alpha_term + beta1_term + beta2_term - beta3_term);
					}
				}

				vec2 n = s.normal;
				ls.add_to_energy( std::sqrt(lambda_distance) * value_res);
				ls.add_to_energy( std::sqrt(lambda_gradient) * (grad_res[0] - n[0]));
				ls.add_to_energy( std::sqrt(lambda_gradient) * (grad_res[1] - n[1]));
			}

			ls.solve();

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

		bool coordinate_descent(double& var, double& step, double lb, double ub, double& curr_err) {
			const double old       = var;
			const auto   old_f = sdf.fonctions;

			var = std::clamp(old + step, lb, ub);
			resolve_ls();
			double errp = error_total();

			if (curr_err - errp > MIN_IMPROVEMENT) { curr_err = errp; return true; }

			var = std::clamp(old - step, lb, ub);
			resolve_ls();
			double errm = error_total();

			if (curr_err - errm > MIN_IMPROVEMENT) { curr_err = errm; return true; }

			var       = old;
			sdf.fonctions = old_f;
			return false;
		}

		bool coordinate_descent_adaptive(double& var, double& step, double lb, double ub, double& ce) {
			const double old       = var;
			const auto   old_f = sdf.fonctions;

			var = std::clamp(old + step, lb, ub);
			resolve_ls();
			double errp = error_total();
			const auto oldp = sdf.fonctions;

			var = std::clamp(old - step, lb, ub);
			resolve_ls();
			double errm = error_total();

			double best_err = std::min(errp, errm);

			if (ce - best_err > MIN_IMPROVEMENT) {
				if (errp <= errm) {
					var       = std::clamp(old + step, lb, ub);
					sdf.fonctions = oldp;
				}
				ce   = best_err;
				step = std::min(step * EXPAND, (ub - lb));
				return true;
			}

			var       = old;
			sdf.fonctions = old_f;
			step      = std::max(step * SHRINK, MIN_STEP);
			return step > MIN_STEP;
		}

		void run_loop(size_t max_it, size_t snapshot, const std::string& output_dir,
		              std::function<bool(size_t)> loop_step) {

			static std::atomic<bool> interupted{false};
			auto prev_handler = std::signal(SIGINT, [](int) { 
					if (interupted.load()) exit(1);
					interupted.store(true); 
				});


			std::ofstream(output_dir + "init.json") << sdf.to_json().dump(2);
			std::cout << std::setprecision(15)
			          << "init: err : " << error_total() << "\t" << sdf.to_string() << "\n";

			for (size_t it = 0; it < max_it; ++it) {
				bool terminated = loop_step(it);

				if (it % (max_it / std::min(snapshot, max_it)) == 0 || terminated || interupted.load()) {
					double total_err = error_total();
					std::cout << std::setprecision(15)
					          << it << ": err : " << total_err << "\t" << sdf.to_string() << std::endl;
					std::ofstream(output_dir + std::to_string(it) + "sdf.json")
					          << sdf.to_json().dump(2);
				}

				if (terminated) {
					std::cout << "Terminated in " << it << " iterations" << std::endl;
					break;
				}
				if (interupted.load()) {
					std::cout << "Interrupted at iteration " << it << std::endl;
					break;
				}
			}
			std::signal(SIGINT, prev_handler);
		}


		bool is_too_close(const UM::vec2& pos) const {
			for (size_t i = 0; i < sdf.fonctions.size(); ++i)
				if ((pos - sdf.fonctions[i].point).norm() < 1e-8) return true;
			return false;
		}

		bool try_add_point(const std::vector<size_t>& idx) {
			if (idx.empty()) return false;
			size_t k = 0;
			while (k < idx.size() && is_too_close(samples[idx[k]].point)) ++k;
			if (k == idx.size()) return false;

			sdf.add_func({samples[idx[k]].point, 1., samples[idx[k]].normal, default_sigma_add});
			return true;
		}
	};

	struct TestCircular : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {

		}
	};

} // namespace sdffitting

#endif
