#ifndef ALGO_HPP__
#define ALGO_HPP__

#include "debug_macros.hpp"
#include "types.hpp"
#include "sdf.hpp"
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
		Samples samples;

		double lb_sigma  = 0.5,  ub_sigma  = 10.,  step_sigma = 0.01;
		double lb_point  = -5.,  ub_point  = 5.,   step_point = 0.01;

		double ADD_POINT_ERR_THRESHOLD = 5.;
		double default_sigma_add       = 1.;

		bool   fix_alpha_zero  = false;
		bool   fix_beta_zero   = false;
		double MIN_IMPROVEMENT = 1e-15;

		// adaptative step
		double EXPAND   = 1.2;
		double SHRINK   = 0.5;
		double MIN_STEP = 0.001;

		explicit Fitter(SDF& _sdf, Samples _samples)
			: sdf(_sdf), samples(std::move(_samples)) {}

		virtual ~Fitter() = default;

		// abstract
		virtual void fit(size_t max_it, size_t snapshot, const std::string& output_dir) = 0;


		void resolve_ls() {
			LeastSquares ls(3 * sdf.p.size());

			if (fix_alpha_zero) {
				for (size_t i = 0; i < sdf.p.size(); ++i)
					ls.fix(i * 3, 0.);
			}
			if (fix_beta_zero) {
				for (size_t i = 0; i < sdf.p.size(); ++i) {
					ls.fix(i * 3 + 1, 0.);
					ls.fix(i * 3 + 2, 0.);
				}
			}

			for (const auto& s : samples) {
				LinExpr value_res;
				LinExpr grad_res[2];

				for (size_t i = 0; i < sdf.p.size(); ++i) {
					vec2   p     = s.first - sdf.p[i];
					double l     = p.norm();
					if (l == 0.) continue;

					double phi   = sdf.rbf->f  (l, sdf.sigma[i]);
					double dphi  = sdf.rbf->df (l, sdf.sigma[i]);
					double ddphi = sdf.rbf->ddf(l, sdf.sigma[i]);

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

				vec2 n = s.second;
				ls.add_to_energy(value_res);
				ls.add_to_energy(grad_res[0] - n[0]);
				ls.add_to_energy(grad_res[1] - n[1]);
			}

			ls.solve();

			for (size_t i = 0; i < sdf.p.size(); ++i) {
				sdf.alpha[i]   = ls.value(i*3);
				sdf.beta[i].x  = ls.value(i*3+1);
				sdf.beta[i].y  = ls.value(i*3+2);
			}
		}

		double error_on_sample(PointNormal& p) const {
			double d   = sdf.distance(p.first);
			vec2   g   = sdf.gradient(p.first);
			vec2   diff = g - p.second;
			return d*d + diff.norm2();
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
			const auto   old_alpha = sdf.alpha;
			const auto   old_beta  = sdf.beta;

			var = std::clamp(old + step, lb, ub);
			resolve_ls();
			double errp = error_total();

			if (curr_err - errp > MIN_IMPROVEMENT) { curr_err = errp; return true; }

			var = std::clamp(old - step, lb, ub);
			resolve_ls();
			double errm = error_total();

			if (curr_err - errm > MIN_IMPROVEMENT) { curr_err = errm; return true; }

			var       = old;
			sdf.alpha = old_alpha;
			sdf.beta  = old_beta;
			return false;
		}

		bool coordinate_descent_adaptive(double& var, double& step, double lb, double ub, double& ce) {
			const double old       = var;
			const auto   old_alpha = sdf.alpha;
			const auto   old_betas = sdf.beta;

			var = std::clamp(old + step, lb, ub);
			resolve_ls();
			double errp = error_total();
			const auto alphap = sdf.alpha;
			const auto betap  = sdf.beta;

			var = std::clamp(old - step, lb, ub);
			resolve_ls();
			double errm = error_total();

			double best_err = std::min(errp, errm);

			if (ce - best_err > MIN_IMPROVEMENT) {
				if (errp <= errm) {
					var       = std::clamp(old + step, lb, ub);
					sdf.alpha = alphap;
					sdf.beta  = betap;
				}
				ce   = best_err;
				step = std::min(step * EXPAND, (ub - lb));
				return true;
			}

			var       = old;
			sdf.alpha = old_alpha;
			sdf.beta  = old_betas;
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

				if (it % (max_it / snapshot) == 0 || terminated || interupted.load()) {
					double total_err = error_total();
					std::cout << std::setprecision(15)
					          << it << ": err : " << total_err << "\t" << sdf.to_string() << "\n";
					std::ofstream(output_dir + std::to_string(it) + "sdf.json")
					          << sdf.to_json().dump(2);
				}

				if (terminated) {
					std::cout << "Terminated in " << it << " iterations\n";
					break;
				}
				if (interupted.load()) {
					std::cout << "Interrupted at iteration " << it << "\n";
					break;
				}
			}
			std::signal(SIGINT, prev_handler);
		}


		bool is_too_close(const UM::vec2& pos) const {
			for (size_t i = 0; i < sdf.p.size(); ++i)
				if ((pos - sdf.p[i]).norm() < 1e-8) return true;
			return false;
		}

		bool try_add_point(const std::vector<size_t>& idx) {
			if (idx.empty()) return false;
			size_t k = 0;
			while (k < idx.size() && is_too_close(samples[idx[k]].first)) ++k;
			if (k == idx.size()) return false;

			sdf.add_func(samples[idx[k]].first, 1., samples[idx[k]].second, default_sigma_add);
			return true;
		}
	};

	// =========================================================================
	// Fitters
	// =========================================================================

	struct FixedStepFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i)
					terminated &= !coordinate_descent(sdf.sigma[i], step_sigma, lb_sigma, ub_sigma, ce);
				return terminated;
			});
		}
	};

	struct FixedStepMovingFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent(sdf.sigma[i], step_sigma, lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent(sdf.p[i].x,   step_point, lb_point, ub_point, ce);
					terminated &= !coordinate_descent(sdf.p[i].y,   step_point, lb_point, ub_point, ce);
				}
				return terminated;
			});
		}
	};

	struct FixedStepAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			if (sdf.p.empty())
				sdf.add_func(samples[0].first, 1., samples[0].second, default_sigma_add);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i)
					terminated &= !coordinate_descent(sdf.sigma[i], step_sigma, lb_sigma, ub_sigma, ce);

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func(samples[idx[0]].first, 1., samples[idx[0]].second, default_sigma_add);
						terminated = false;
					}
				}
				return terminated;
			});
		}
	};

	struct FixedStepMovingAddPointFitter : Fitter {
		using Fitter::Fitter;

		double last_add_err = DBL_MAX;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			if (sdf.p.empty())
				sdf.add_func(samples[0].first, 1., samples[0].second, default_sigma_add);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent(sdf.sigma[i], step_sigma, lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent(sdf.p[i].x,   step_point, lb_point, ub_point, ce);
					terminated &= !coordinate_descent(sdf.p[i].y,   step_point, lb_point, ub_point, ce);
				}

				if (terminated) {
					std::vector<size_t> idx;
					double err = error_total(&idx);
					if (err > ADD_POINT_ERR_THRESHOLD && last_add_err - err > ADD_POINT_ERR_THRESHOLD) {
						last_add_err = err;
						sdf.add_func(samples[idx[0]].first, 1., samples[idx[0]].second, default_sigma_add);
						terminated = false;
					}
				}
				return terminated;
			});
		}
	};

	struct AdaptiveStepFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.p.size(), step_sigma);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i)
					terminated &= !coordinate_descent_adaptive(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
				return terminated;
			});
		}
	};

	struct AdaptiveStepAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.p.size(), step_sigma);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i)
					terminated &= !coordinate_descent_adaptive(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func(samples[idx[0]].first, 1., samples[idx[0]].second, default_sigma_add);
						steps_sigma.push_back(step_sigma);
						terminated = false;
					}
				}
				return terminated;
			});
		}
	};

	struct AdaptiveStepMovingFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.p.size(), step_sigma);
			std::vector<double> steps_px   (sdf.p.size(), step_point);
			std::vector<double> steps_py   (sdf.p.size(), step_point);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent_adaptive(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptive(sdf.p[i].x,   steps_px[i],   lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptive(sdf.p[i].y,   steps_py[i],   lb_point, ub_point, ce);
				}
				return terminated;
			});
		}
	};

	struct AdaptiveStepMovingAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.p.size(), step_sigma);
			std::vector<double> steps_px   (sdf.p.size(), step_point);
			std::vector<double> steps_py   (sdf.p.size(), step_point);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent_adaptive(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptive(sdf.p[i].x,   steps_px[i],   lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptive(sdf.p[i].y,   steps_py[i],   lb_point, ub_point, ce);
				}

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func(samples[idx[0]].first, 1., samples[idx[0]].second, default_sigma_add);
						steps_sigma.push_back(step_sigma);
						steps_px.push_back(step_point);
						steps_py.push_back(step_point);
						terminated = false;
					}
				}
				return terminated;
			});
		}
	};

	struct AlphaBetaOnlyAddFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {

			int K = 2;

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
					resolve_ls();

					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD && it < (max_it-1)) {

						auto knn = UM::KNN<2>(sdf.p);

						auto p = samples[idx[0]].first;
						auto n = samples[idx[0]].second;

						double new_sigma = 0.;

						for (const int& i : knn.query(p, K)) {
							new_sigma += (p - sdf.p[i]).norm();

							double modif_sigma = 0.;

							for (const auto& j : knn.query(sdf.p[i], K)) {
								modif_sigma += (sdf.p[i] - sdf.p[j]).norm();	
							}
							
							sdf.sigma[i] = modif_sigma / K;
						}


						sdf.add_func(p, 0., n, new_sigma/ K);
						
						return false;
					}

					return true; // below threshold → done
					});



		}
	};

} // namespace sdffitting

#endif
