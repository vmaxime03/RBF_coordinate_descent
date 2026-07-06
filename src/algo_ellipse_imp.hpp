#ifndef ALGO_ELLIPSE_IMP__
#define ALGO_ELLIPSE_IMP__

#include <algorithm>
#include <cfloat>
#include <cstddef>
#include <numbers>
#include <sstream>
#include <vector>
#include "algo_base_elliptic.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/knn.h"

#include "debug_macros.hpp"
#include "ultimaille/polyline.h"
#include "ultimaille/primitive_geometry.h"
#include "ultimaille/all.h"


namespace sdffitting_elliptic {


		struct TODONAMEFitter : Fitter_Elliptic {
			using Fitter_Elliptic::Fitter_Elliptic;

			void adapt_minor(double min_radius_coef = 1.0) {

#pragma omp parallel for schedule(guided)
			for (size_t i = 0; i < sdf.fonctions.size(); ++i) {

				auto& f = sdf.fonctions[i];
				auto final_new_minor = f.ellipse_major; // start with cicle and then shrink


				// for (const auto& j : neigh) {
				for (size_t j = 0; j < sdf.fonctions.size(); ++j) {
					if (j == i) continue;

					// find the maximal minor radius that avoid interaction with function that have opposites normals
					// assertion : minor radius have same direction as normals 

					auto d = sdf.fonctions[j].point - f.point;

					auto n1 = f.ellipse_minor;
					auto n2 = sdf.fonctions[j].ellipse_minor;

					auto n1_norm = n1.normalized();
					auto n1_max = n1_norm * f.ellipse_major.norm();
					auto n2_norm = n2.normalized();

					double cos = n1_norm * n2_norm;

					auto hd = d / min_radius_coef; // minimal minor radius

					auto m = std::abs(hd * n1_norm) * n1_norm;

					auto r = n1_max - m;

					double scale = (cos + 1.0) / 2.0;

					auto new_minor = m + r * scale;
					if (new_minor.norm() < final_new_minor.norm()) final_new_minor = new_minor;

				}

				f.ellipse_minor = final_new_minor;
			}
			}



			struct triplet { vec2 centroid, major, minor; };

			inline triplet as_triplet(const FunctionElliptic& f) {
				return {f.point, f.ellipse_major, f.ellipse_minor};
			};
			bool are_ellipse_colliding(size_t id1, size_t id2) {
				// we assume that every ellipse are per segment with major being longer than segment
				auto e1 = sdf.fonctions[id1];
				auto e2 = sdf.fonctions[id2];

				// check if major radius "segment" cross
				Segment2 s1 = {e1.point + e1.ellipse_major, e1.point - e1.ellipse_major};
				Segment2 s2 = {e2.point + e2.ellipse_major, e2.point - e2.ellipse_major};

				auto _cross = [](vec2 a, vec2 b, vec2 c) {
					return UM::cross( (b - a).xy0(), (c-a).xy0()).z;
				};
				auto o1 = _cross(s1.a, s1.b, s2.a);
				auto o2 = _cross(s1.a, s1.b, s2.b);
				auto o3 = _cross(s2.a, s2.b, s1.a);
				auto o4 = _cross(s2.a, s2.b, s1.b);
				
				// check segme:q
				// nt overlap
				auto _overlap = [](double a, double b, double c, double d) {
					return std::max(std::min(a, b), std::min(c, d)) <= std::min(std::max(a, b), std::max(c, d));
				};

				if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) {
					return _overlap(s1.a.x, s1.b.x, s2.a.x, s2.b.x) && _overlap(s1.a.y, s1.b.y, s2.a.y, s2.b.y);
				}

				return (o1 * o2 <= 0) && (o3 * o4 <= 0);
			};





			void decimate(double angle_threshold, bool pass2) {

				std::vector<std::vector<size_t>> clusters;
				std::vector<vec2> clusters_angle;

			

				auto is_in_angular_threshold = [&](const vec2& n1, const vec2& n2) -> bool {
					return n1.normalized() * n2.normalized() >= std::cos(angle_threshold);
				};

				auto compute_cluster_angle = [&](size_t cid) -> void {
					clusters_angle.resize(clusters.size());
					vec2 t = {0., 0.};
					for (size_t idx : clusters[cid])
						t += sdf.fonctions[idx].ellipse_minor;
					clusters_angle[cid] = t / clusters[cid].size();
				};


				// // VERBOSE  
				// static std::atomic<int> decimate2_call_counter{0};
				// int call_id = decimate2_call_counter++;
				// std::string csv_path = output_dir + "decimate2_clusters_" + std::to_string(call_id) + ".csv";
				// std::ofstream csv(csv_path);
				// csv << "step,action,element_id,cluster_id,px,py\n";
				// int step = 0;
				//
				// // Flush the current cluster state for every ellipse to the CSV
				// auto dump_state = [&](const std::string& action, size_t changed_element) {
				// 	for (size_t cl = 0; cl < clusters.size(); ++cl) {
				// 		for (size_t idx : clusters[cl]) {
				// 			const auto& pt = sdf.fonctions[idx].point;
				// 			csv << step << ","
				// 				<< action << ","
				// 				<< changed_element << ","
				// 				<< cl << ","
				// 				<< pt.x << ","
				// 				<< pt.y << "\n";
				// 		}
				// 	}
				// 	++step;
				// };



				// Pass 1	
				clusters.push_back({0});
				compute_cluster_angle(0);

				for (size_t i = 1; i < sdf.fonctions.size(); ++i) {
					auto& n = sdf.fonctions[i];
					size_t cid = clusters.size() - 1;

					bool compatible = is_in_angular_threshold(n.ellipse_minor, clusters_angle[cid])
						&& are_ellipse_colliding(i, clusters[cid].back());

					if (compatible) {
						clusters[cid].push_back(i);
					} else {
						clusters.push_back({i});
						cid = clusters.size() - 1;
					}
					compute_cluster_angle(cid);

					// dump_state("pass1", i);
				}

//
// 				if (pass2) { // todo PARALELIZE
// 				// pass2 : merge cluster that should be merged
// 				for (size_t i = 0; i < clusters.size(); ++i) {
// 					for (size_t j = 0; j < clusters.size(); ++j) {
// 						DEBUG(i << " / " << clusters.size() << ", " << j);
// 						if (i == j) continue;
// 						auto last_id  = clusters[i].back();
// 						auto first_id = clusters[j].front();
// 						if (are_ellipse_colliding(last_id, first_id)
// 								&& is_in_angular_threshold(sdf.fonctions[last_id].ellipse_minor,sdf.fonctions[first_id].ellipse_minor)
// 								&& is_in_angular_threshold(clusters_angle[i], clusters_angle[j])
// 						   ) {
// 							clusters[i].insert(clusters[i].end(),
// 									clusters[j].begin(),
// 									clusters[j].end());
// 							clusters.erase(clusters.begin() + j);
// 							clusters_angle.erase(clusters_angle.begin() + j);
//
// 							compute_cluster_angle(i > 0 ? i : 0);  // recompute merged cluster angle
// // dump_state("pass2_merge", last_id);     //  dump after every merge
// 							i = 0; j = 0; // restart
// 						}
// 					}
// 				}
// 				}


				std::vector<triplet> new_points(clusters.size());


				for (size_t cluster_id = 0; cluster_id < clusters.size(); ++cluster_id) {
					vec2 centroid = {0, 0}, minor = {0, 0};
					for (size_t j : clusters[cluster_id]) {
						centroid += sdf.fonctions[j].point;
						minor    += sdf.fonctions[j].ellipse_minor;
					}
					centroid /= clusters[cluster_id].size();
					minor    /= clusters[cluster_id].size();

					double major_radius = -DBL_MAX;
					for (size_t j : clusters[cluster_id]) {
						auto& fj = sdf.fonctions[j];
						major_radius = std::max(major_radius, (centroid - (fj.point + fj.ellipse_major)).norm());
						major_radius = std::max(major_radius, (centroid - (fj.point - fj.ellipse_major)).norm());
					}
					auto major = vec2(-minor.y, minor.x).normalized() * major_radius;
					new_points[cluster_id] = {centroid, major, minor};


				}

				DEBUG(sdf.fonctions.size());
				DEBUG(clusters.size());

				sdf.fonctions.resize(clusters.size());
				sdf.active.resize(clusters.size());
				for (size_t i = 0; i < clusters.size(); ++i) {
					auto& c = new_points[i];
					sdf.fonctions[i] = {{c.centroid, 0., {0., 0.}}, c.major, c.minor};
					sdf.active[i] = true;
				}


				// csv.close();
				// std::cout << "[decimate2] cluster evolution written to " << csv_path << "\n";
			}


			void adjust_minor(double radius_coef = 1.0) {
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					sdf.fonctions[i].ellipse_minor *= radius_coef;
				}

			}

			
			void adjust_major(double radius_coef = 1.0) {
				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					sdf.fonctions[i].ellipse_major *= radius_coef;
				}
			}


			// ELLIPSE FITTING ============================================================

			void fit_ellipses_radius(double W, int K = 2, int min_K = 2) {

				int N = sdf.fonctions.size();

				std::vector<double> X(N * 2);

				std::vector<vec2> A_dirs(N);
				std::vector<vec2> B_dirs(N);

#pragma omp parallel for
				for (int i = 0; i < N; ++i) {
					X[i*2 + 0] = sdf.fonctions[i].ellipse_major.norm();
					X[i*2 + 1] = sdf.fonctions[i].ellipse_minor.norm();

					A_dirs[i] = sdf.fonctions[i].ellipse_major.normalized();
					B_dirs[i] = sdf.fonctions[i].ellipse_minor.normalized();
				}

				// precompute neighborhoods

				if (!sdf.optimized) sdf.init_knn();
				std::vector<std::vector<int>> neighborhoods(N);

				for (int i = 0; i < N; ++i) {
					auto& fi = sdf.fonctions[i];

					// TODO maybe query more than K and filter out opposed ellipse
					auto neighborhood = sdf.knn->query(fi.point, 1 + K);

					std::vector<int> filtered;
					filtered.reserve(neighborhood.size());

					int right = 0, left = 0;
					for (int n : neighborhood) {
						if (n == i) {
							continue;
						}
						auto& fn = sdf.fonctions[n];

						if ((fn.point - fi.point) * fi.ellipse_major > 0) { // right side neighbor
							++right;
							if (right > K/2 && filtered.size() > static_cast<size_t>(min_K)) {
								continue;
							}
						} else { // left side neighbor
							++left;
							if (left > K/2 && filtered.size() > static_cast<size_t>(min_K)) {
								continue;
							}
						}

						filtered.push_back(n);

						if (filtered.size() > static_cast<size_t>(K)) break;

					}

					neighborhoods[i] = filtered;
				}


				constexpr double lambda_W = 1000.;
				constexpr double lambda_Area = 1.;
				constexpr double lambda_K = 100.;

				double lambda_K_weigth = (lambda_K / (2*K));


				// energy : minor value = W + ellipse area + contains Ci-1 & Ci+1

				const auto func = [&](const std::vector<double>&x, double& f, std::vector<double>& g) {
					f = 0.0; 
					std::fill(g.begin(), g.end(), 0.0);

					// ellise area 
#pragma omp parallel for reduction(+:f)
					for (int i = 0; i < N; ++i) {
						auto& fi = sdf.fonctions[i];

						auto a = x[2*i];
						auto b = x[2*i+1];

						auto& g1 = g[2*i];
						auto& g2 = g[2*i+1];

						// 
						f += lambda_W * ((b - W)*(b - W));

						g1 += 0;
						g2 += lambda_W * 2 * (b - W);

						//
						f += lambda_Area * (a*a*b*b);

						g1 += lambda_Area * 2 * a * b * b;
						g2 += lambda_Area * 2 * b * a * a;


						for (int j : neighborhoods[i]) {
							if (i == j) continue;

							auto& fj = sdf.fonctions[j];

							double cosinus_normalized = ((B_dirs[i] * B_dirs[j]) + 1.0) / 2.0;
							double dir_weigth = cosinus_normalized * cosinus_normalized;
							double current_lambda_K = lambda_K_weigth * dir_weigth;

							for (auto dir : {-1.0, 0.0}) {
								// auto point_to_contain = fj.point;
								auto point_to_contain = fj.point + dir * fj.ellipse_minor.normalized() * (W * 0.5) * cosinus_normalized;

								auto diff = point_to_contain - fi.point;

								double x_prime = diff * A_dirs[i];
								double y_prime = diff * B_dirs[i];

								double ell_dist = ((x_prime * x_prime) / (a*a) + (y_prime * y_prime) / (b*b) - 1.0);
								// f += max(0, ell_dist)^2
								if (ell_dist > 0.) {
									f += current_lambda_K * (ell_dist * ell_dist);

									g1 += current_lambda_K * 2.0 * ell_dist * (-2.0 * x_prime * x_prime / (a * a * a));
									g2 += current_lambda_K * 2.0 * ell_dist * (-2.0 * y_prime * y_prime / (b * b * b));

								}
							}

						}
					}



				};

				double E_prev, E;
				std::vector<double> trash(X.size());
				func(X, E_prev, trash);


				STLBFGS::Optimizer opt(func);
				opt.run(X);

				func(X, E, trash);
    			DEBUG("E: " <<  E_prev << " --> " << E);

				for (int i = 0; i < N; ++i) {
					sdf.fonctions[i].ellipse_major = X[i*2] * A_dirs[i];
					sdf.fonctions[i].ellipse_minor = X[i*2+1] * B_dirs[i];
				}



			}


			void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {

			}

		};


	}




#endif 
