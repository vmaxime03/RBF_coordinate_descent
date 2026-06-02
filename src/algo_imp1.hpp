#ifndef ALGO_IMP1__
#define ALGO_IMP1__

#include "algo_base.hpp"

namespace sdffitting {
	// =========================================================================
	// Fitters
	// =========================================================================

	struct FixedStepFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i)
					terminated &= !coordinate_descent(sdf.fonctions[i].sigma, step_sigma, lb_sigma, ub_sigma, ce);
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
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					terminated &= !coordinate_descent(sdf.fonctions[i].sigma, step_sigma, lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent(sdf.fonctions[i].point.x,   step_point, lb_point, ub_point, ce);
					terminated &= !coordinate_descent(sdf.fonctions[i].point.y,   step_point, lb_point, ub_point, ce);
				}
				return terminated;
			});
		}
	};

	struct FixedStepAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			if (sdf.fonctions.empty())
				sdf.add_func({samples[0].point, 1., samples[0].normal, default_sigma_add});

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i)
					terminated &= !coordinate_descent(sdf.fonctions[i].sigma, step_sigma, lb_sigma, ub_sigma, ce);

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func({samples[idx[0]].point, 1., samples[idx[0]].normal, default_sigma_add});
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
			if (sdf.fonctions.empty())
				sdf.add_func({samples[0].point, 1., samples[0].normal, default_sigma_add});

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					terminated &= !coordinate_descent(sdf.fonctions[i].sigma, step_sigma, lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent(sdf.fonctions[i].point.x,   step_point, lb_point, ub_point, ce);
					terminated &= !coordinate_descent(sdf.fonctions[i].point.y,   step_point, lb_point, ub_point, ce);
				}

				if (terminated) {
					std::vector<size_t> idx;
					double err = error_total(&idx);
					if (err > ADD_POINT_ERR_THRESHOLD && last_add_err - err > ADD_POINT_ERR_THRESHOLD) {
						last_add_err = err;
						sdf.add_func({samples[idx[0]].point, 1., samples[idx[0]].normal, default_sigma_add});
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
			std::vector<double> steps_sigma(sdf.fonctions.size(), step_sigma);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i)
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].sigma, steps_sigma[i], lb_sigma, ub_sigma, ce);
				return terminated;
			});
		}
	};

	struct AdaptiveStepAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.fonctions.size(), step_sigma);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i)
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].sigma, steps_sigma[i], lb_sigma, ub_sigma, ce);

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func({samples[idx[0]].point, 1., samples[idx[0]].normal, default_sigma_add});
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
			std::vector<double> steps_sigma(sdf.fonctions.size(), step_sigma);
			std::vector<double> steps_px   (sdf.fonctions.size(), step_point);
			std::vector<double> steps_py   (sdf.fonctions.size(), step_point);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].sigma, steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].point.x,   steps_px[i],   lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].point.y,   steps_py[i],   lb_point, ub_point, ce);
				}
				return terminated;
			});
		}
	};

	struct AdaptiveStepMovingAddPointFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {
			std::vector<double> steps_sigma(sdf.fonctions.size(), step_sigma);
			std::vector<double> steps_px   (sdf.fonctions.size(), step_point);
			std::vector<double> steps_py   (sdf.fonctions.size(), step_point);

			run_loop(max_it, snapshot, output_dir, [&](size_t it) {
				bool   terminated = true;
				double ce         = error_total();
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].sigma, steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].point.x,   steps_px[i],   lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptive(sdf.fonctions[i].point.y,   steps_py[i],   lb_point, ub_point, ce);
				}

				if (terminated) {
					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {
						sdf.add_func({samples[idx[0]].point, 1., samples[idx[0]].normal, default_sigma_add});
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
}

#endif
