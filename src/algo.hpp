#ifndef ALGO_HPP__
#define ALGO_HPP__


#include "output.hpp"
#include "sdf.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include "ultimaille/sparse/least_squares.h"
#include "ultimaille/sparse/linexpr.h"
#include <cfloat>
#include <fstream>
#include <vector>


namespace sdffitting { // TODO

	using namespace UM;
	using namespace UM::Linear;


	typedef std::pair<vec2, vec2> PointNormal;
	typedef std::vector<PointNormal> Samples;

	Samples compute_samples_normals(PolyLine& pl, size_t n = 10) {
		Samples samples;
		for (const auto& e : pl.iter_edges()) {
			vec2 d = e.to().pos().xy() - e.from().pos().xy();

			vec2 normal = vec2(-d.y, d.x).normalized();

			for (size_t i = 0; i < n; ++i) {
				vec2 p = e.from().pos().xy() + ((double)i / (double)n) * d;
				samples.push_back({p, normal});
			}
		}
		return samples;
	}

	struct Fitter {

		PolyLine& pl;
		SDF& sdf;
		Samples samples;


		explicit Fitter(	
				PolyLine& _pl,
				SDF& _sdf,
				Samples _samples
				) 
			: pl(_pl), sdf(_sdf), samples(_samples) {}



		double lb_sigma = 0.5,	ub_sigma = 10.,  	step_sigma = 0.01;
		double lb_point = -5.,	ub_point = 5.,  	step_point = 0.01;

		double ADD_POINT_ERR_THRESHOLD = 5.;	


		void resolve_ls() {
			// LEAST SQUARE
			// find alphas & betas
			LeastSquares ls(3*sdf.p.size());

			for (const auto& s : samples) {

				LinExpr value_res;
				LinExpr grad_res[2];

				for (size_t i = 0; i < sdf.p.size(); ++i) { 

					vec2 p = s.first - sdf.p[i]; // x - p_i 
					double l = p.norm();

					// TODO mimic compact support ?  
					if (l == 0.) continue;

					double phi 		= sdf.rbf->f  (l, sdf.sigma[i]);
					double dphi 	= sdf.rbf->df (l, sdf.sigma[i]);
					double ddphi	= sdf.rbf->ddf(l, sdf.sigma[i]);

					// f(p_i) == 0
					LinExpr term = X(i*3) * phi +	X(i*3+1) * (dphi * p.x/p.norm()) + X(i*3+2) * (dphi * p.y/p.norm());

					value_res += term;

					// grad f(p_i) == n_i
					auto dotprod = [&](size_t i) { return (X(i*3 + 1) * p[0] + X(i*3 + 2) * p[1]);};

					for (size_t k = 0; k < 2; ++k) {

						LinExpr alpha_term = X(i*3) * dphi * ( p[k] / l );
						LinExpr beta1_term = ddphi * ( p[k] / (l * l) ) * dotprod(i);
						LinExpr beta2_term = (dphi / l) * X(i*3 + 1 + k) ;
						LinExpr beta3_term = ((dphi * p[k])/(l*l*l)) * dotprod(i);

						grad_res[k] += ((alpha_term + beta1_term + beta2_term - beta3_term));
					}
				}

				vec2 n = s.second; // n_i

				ls.add_to_energy(value_res);
				ls.add_to_energy(grad_res[0] - n[0]);
				ls.add_to_energy(grad_res[1] - n[1]);
			}

			ls.solve();

			for (size_t i = 0; i < sdf.p.size(); ++i) { 

				double a = ls.value(i*3);
				double bx = ls.value(i*3+1);
				double by = ls.value(i*3+2);

				sdf.alpha[i]  = a;
				sdf.beta[i].x = bx;
				sdf.beta[i].y = by;
			}
		};

		double error_on_polyline(PointNormal& p) {
			double d = sdf.distance(p.first); // distance 
			vec2 g = sdf.gradient(p.first); // gradiant

			vec2 diff = (g - p.second);

			return d * d + diff.norm2();
		}


		// total error score
		double error_total(size_t* pemax = nullptr) {
			double t = 0.;

			double emax = -DBL_MAX;

			for (size_t i = 0; i < samples.size(); ++i) {
				double err = error_on_polyline(samples[i]);	

				if (pemax && err > emax) {
					emax = err;
					*pemax = (int)i;
				}

				t += err;
			}

			return t;
		}

		bool coordinate_descent(
				double& var,  // coordiante to optimize
				double& step, // step for the coord
				double lb, double ub, // bound for the coord
				double& curr_err // current error score
				) {

			const double old = var;

			const auto old_alpha = sdf.alpha;
			const auto old_beta  = sdf.beta;


			// try +step
			var = std::clamp(old + step, lb, ub);
			resolve_ls();

			double errp = error_total();

			if (errp < curr_err) { 
				curr_err = errp;
				return true;
			}

			// try -step
			var = std::clamp(old - step, lb, ub);
			resolve_ls();

			double errm = error_total();

			if (errm < curr_err) { 
				curr_err = errm;
				return true;
			}

			
			var = old;
			sdf.alpha = old_alpha;
			sdf.beta = old_beta;

			return false;
		}

		void run_loop(size_t max_it, size_t snapshot, const std::string& output_dir, std::function<bool()> loop_step) {
			bool terminated = false;
			for (size_t it = 0; it < max_it; ++it) {

				// SNAPSHOT OUTPUT	
				if (it % (max_it / snapshot) == 0 || terminated) {
					double total_err = error_total();
					std::cout << it << ": err : " << total_err << "\t"  << sdf.to_string() << std::endl;
					std::ofstream(output_dir + std::to_string(it) + "sdf.json") << sdf.to_json().dump(2);
				}

				if (terminated) {
					std::cout << "Terminated in " << it << " iterations" << std::endl;
					break;
				}

				terminated = loop_step();
			}

			// OUTPUT
			export_polyline(pl, output_dir + "polyline.csv");
			export_sdf(sdf, output_dir + "sdf_params.csv");
			sample_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");
		}


		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) {
			run_loop(max_it, snapshot, output_dir, 
					[&]() {	
						bool terminated = true;
						double ce = error_total();
						for (size_t i = 0; i < sdf.p.size(); ++i) {

							terminated &= !coordinate_descent(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, ce);

						}
						return terminated;
					});

		}


		void fit_moving_points(size_t max_it, size_t snapshot, const std::string& output_dir) {
			run_loop(max_it, snapshot, output_dir, 
					[&]() {	
						bool terminated = true;
						double ce = error_total();

						for (size_t i = 0; i < sdf.p.size(); ++i) {

							terminated &= !coordinate_descent(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, ce);

							terminated &= !coordinate_descent(sdf.p[i].x,  step_point, lb_point, ub_point, ce);
							terminated &= !coordinate_descent(sdf.p[i].y,  step_point, lb_point, ub_point, ce);

						}
						return terminated;
					});

		}


		// TODO
		void fit_add_point(size_t max_it, size_t snapshot, const std::string& output_dir) {

			sdf.add_func(samples[0].first, 1., samples[0].second, 1.);

			run_loop(max_it, snapshot, output_dir, 
					[&]() {	
					bool terminated = true;
					double ce = error_total();


					for (size_t i = 0; i < sdf.p.size(); ++i) {

						terminated &= !coordinate_descent(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, ce);

					}

					if (terminated) {
						size_t idx;
						double err = error_total(&idx);

						if (err > ADD_POINT_ERR_THRESHOLD) {
							sdf.add_func(samples[idx].first, 1, samples[idx].second, 1.);
							terminated = false;
						}
					}
					return terminated;
					});


		}



		void fit_moving_add_points(size_t max_it, size_t snapshot, const std::string& output_dir) {

			double last_add_err = DBL_MAX;
			sdf.add_func(samples[0].first, 1., samples[0].second, 1.);


			run_loop(max_it, snapshot, output_dir, 
					[&]() {	
						bool terminated = true;
						double ce = error_total();
						for (size_t i = 0; i < sdf.p.size(); ++i) {


							terminated &= !coordinate_descent(sdf.sigma[i],  step_sigma, lb_sigma, ub_sigma, ce);

							terminated &= !coordinate_descent(sdf.p[i].x,  step_point, lb_point, ub_point, ce);
							terminated &= !coordinate_descent(sdf.p[i].y,  step_point, lb_point, ub_point, ce);

						}

						if (terminated) {
							size_t idx;
							double err = error_total(&idx);

							if (err > ADD_POINT_ERR_THRESHOLD && last_add_err - err > ADD_POINT_ERR_THRESHOLD) {
								last_add_err = err;
								sdf.add_func(samples[idx].first, 1, samples[idx].second, 1.);
								terminated = false;
							}
						}
						return terminated;
					});

		}


		// ADAPTATIVE STEP FITTER

		double EXPAND = 1.2;
		double SHRINK = 0.5;
		double MIN_STEP = 0.001;



		bool coordinate_descent_adaptative(
				double& var, 
				double& step, 
				double lb, double ub, 
				double& ce
				) {
				
			const double old = var;
			const auto old_alpha = sdf.alpha;
			const auto old_betas = sdf.beta;

			// try +
			var = std::clamp(old + step, lb, ub);
			resolve_ls();
			double errp = error_total();

			const auto alphap = sdf.alpha;
			const auto betap = sdf.beta;


			// try -
			var = std::clamp(old - step, lb, ub);
			resolve_ls();
			double errm = error_total();

			double best_err = std::min(errp, errm);


			// improvment case 
			if (best_err < ce) {

				if (errp <= errm) {
					var = std::clamp(old + step, lb, ub);
					sdf.alpha = alphap;
					sdf.beta = betap;
				}

				ce = best_err;
				step = std::min(step * EXPAND, (ub - lb));
				return true;
			}

			// no improvment 
			var = old;
			sdf.alpha = old_alpha;
			sdf.beta = old_betas;

			step = std::max(step * SHRINK, MIN_STEP);
			return step > MIN_STEP;
		}

		void fit_adaptative_step(size_t max_it, size_t snapshot, const std::string& output_dir) {
			std::vector<double> steps_sigma;
			steps_sigma.assign(sdf.p.size(), step_sigma);

			run_loop(max_it, snapshot, output_dir, [&]() {
					bool terminated = true;
					double ce = error_total();
					for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent_adaptative(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
					}
					return terminated;
					});
		}

		void fit_adaptative_step_moving_points(size_t max_it, size_t snapshot, const std::string& output_dir) {
			std::vector<double> steps_sigma;
			std::vector<double> steps_px;
			std::vector<double> steps_py;
			steps_sigma.assign(sdf.p.size(), step_sigma);
			steps_px.assign(sdf.p.size(), step_point);
			steps_py.assign(sdf.p.size(), step_point);

			run_loop(max_it, snapshot, output_dir, [&]() {
					bool terminated = true;
					double ce = error_total();
					for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent_adaptative(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptative(sdf.p[i].x, steps_px[i], lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptative(sdf.p[i].y, steps_py[i], lb_point, ub_point, ce);
					}
					return terminated;
					});
		}

		void fit_adaptative_step_moving_add_points(size_t max_it, size_t snapshot, const std::string& output_dir) {

			double last_add_err = DBL_MAX;
			sdf.add_func(samples[0].first, 1., samples[0].second, 1.);

			std::vector<double> steps_sigma;
			std::vector<double> steps_px;
			std::vector<double> steps_py;
			steps_sigma.assign(sdf.p.size(), step_sigma);
			steps_px.assign(sdf.p.size(), step_point);
			steps_py.assign(sdf.p.size(), step_point);

			run_loop(max_it, snapshot, output_dir, 
					[&]() {	
					bool terminated = true;
					double ce = error_total();
					for (size_t i = 0; i < sdf.p.size(); ++i) {
					terminated &= !coordinate_descent_adaptative(sdf.sigma[i], steps_sigma[i], lb_sigma, ub_sigma, ce);
					terminated &= !coordinate_descent_adaptative(sdf.p[i].x, steps_px[i], lb_point, ub_point, ce);
					terminated &= !coordinate_descent_adaptative(sdf.p[i].y, steps_py[i], lb_point, ub_point, ce);
					}

					if (terminated) {
					size_t idx;
					double err = error_total(&idx);

					if (err > ADD_POINT_ERR_THRESHOLD && last_add_err - err > ADD_POINT_ERR_THRESHOLD) {
					last_add_err = err;
					sdf.add_func(samples[idx].first, 1, samples[idx].second, 1.);

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
