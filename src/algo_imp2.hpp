#ifndef ALGO_IMP2__
#define ALGO_IMP2__ 

#include "algo_base.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/knn.h"
#include <algorithm>
#include <cstddef>
#include <numbers>
#include <ostream>
#include "output.hpp"


namespace sdffitting {
	struct AlphaBetaOnlyAddFitter : Fitter {
		using Fitter::Fitter;

		void fit(size_t max_it, size_t snapshot, const std::string& output_dir) override {

			resolve_ls();
			int K = 2;

			run_loop(max_it-1, snapshot, output_dir, [&](size_t it) {

					std::vector<size_t> idx;
					if (error_total(&idx) > ADD_POINT_ERR_THRESHOLD) {


					std::vector<vec2> pts(sdf.fonctions.size());
					std::transform(sdf.fonctions.begin(), sdf.fonctions.end(), pts.begin(), [](const auto& f) {
							return f.point;
							});

					auto knn = UM::KNN<2>(pts);

					auto p = samples[idx[0]].point;
					auto n = samples[idx[0]].normal;

					double new_sigma = 0.;

					for (const int& i : knn.query(p, K)) {
					new_sigma += (p - sdf.fonctions[i].point).norm();

					double modif_sigma = 0.;

					for (const auto& j : knn.query(sdf.fonctions[i].point, K)) {
						modif_sigma += (sdf.fonctions[i].point - sdf.fonctions[j].point).norm();	
					}

					sdf.fonctions[i].sigma = modif_sigma / K;
					}


					sdf.add_func({p, 0., n, new_sigma/ K});

					resolve_ls();

					return false;
					}

					return true; // below threshold → done
			});



		}
	};


	struct ClusteringFitter : Fitter {
		using Fitter::Fitter;

		size_t K = 2;
		double margin = 1.;

		// update each point sigma so they have only K point in their support
		void adapt_sigma() {

			std::vector<vec2> pts(sdf.fonctions.size());
			std::transform(sdf.fonctions.begin(), sdf.fonctions.end(), pts.begin(), [](const auto& f) {
					return f.point;
					});

			auto knn = UM::KNN<2>(pts);

			for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
				auto neigh = knn.query(sdf.fonctions[i].point, std::min(K, sdf.fonctions.size()));

				double max_dist = 0.;

				for (const auto& j : neigh) {
					max_dist = std::max(max_dist, (sdf.fonctions[i].point - sdf.fonctions[j].point).norm());
				}

				sdf.fonctions[i].sigma = max_dist > 0. ? max_dist * margin : default_sigma_add; // prevent 0 sigma when only 1 point
			}
		}

		vec2 nearest_point_on_polyline(const vec2& pos) { // TODO
			std::vector<vec2> points(samples.size());
			std::transform(samples.begin(), samples.end(), points.begin(), [](const auto& s) { return s.point; });

			auto knn = UM::KNN<2>(points);

			return points[knn.query(pos, 1)[0]];
		}


		// CLUSTERING ERROR SAMPLES 

		double angular_threshold = std::numbers::pi/30; // maximum difference in normal angles between samples in a cluster
		double error_threshold = 0.2; // minimal error 
		double min_distance_points = 0.1;

		size_t cluster_min_size = 2;

		std::vector<vec2> where_to_add_point() {

			// Samples error computation
			std::vector<double> error_samples(samples.size());

			for (size_t i = 0; i < samples.size(); ++i) {
				error_samples[i] = error_on_sample(samples[i]);
			}

			// effective threshold for ignored error
			double effective_threshold = *std::max_element(error_samples.begin(), error_samples.end()) * error_threshold;


			std::vector<std::vector<size_t>> clusters;
			std::vector<size_t> current;

			auto push_cluster = [&]() {
				if (!current.empty()) {
					clusters.push_back(current);
					current.clear();
				}
			};

			// clustering
			size_t first_loop_cluster = 0;
			size_t i = 0;
			while (i < samples.size()) {
				double err = error_samples[i];

				// if error is not big enough
				if (err < effective_threshold) {
					push_cluster();
					++i;
					continue;
				}

				// if angle difference is too big
				if (!current.empty()) {
					size_t prev = current.back();
					double angle = (samples[i].normal * samples[prev].normal) / (samples[i].normal.norm() * samples[prev].normal.norm());
					if (angle < std::cos(angular_threshold)) {
						push_cluster();
					}
				}

				current.push_back(i);

				// check if change loop
				if (samples[i].next != (i+1)) {
					// check if last cluster should be merged with first cluster 
					if (clusters.size() > first_loop_cluster && !current.empty()) {
						size_t last = current.back();
						size_t first = clusters[first_loop_cluster].front();

						double angle = (samples[first].normal * samples[last].normal) / (samples[first].normal.norm() * samples[last].normal.norm());

						if (angle >= std::cos(angular_threshold)) {
							// merge clusters
							clusters[first_loop_cluster].insert(clusters[first_loop_cluster].begin(), current.begin(), current.end());
							current.clear();
						}
					}

					push_cluster();

					first_loop_cluster = clusters.size();
				}
				++i;
			}

			// handler last cluster merge
			// check if last cluster should be merged with first cluster 
			if (clusters.size() > first_loop_cluster && !current.empty()) {
				size_t last = current.back();
				size_t first = clusters[first_loop_cluster].front();

				double angle = (samples[first].normal * samples[last].normal) / (samples[first].normal.norm() * samples[last].normal.norm());

				if (angle >= std::cos(angular_threshold)) {
					// merge clusters
					clusters[first_loop_cluster].insert(clusters[first_loop_cluster].begin(), current.begin(), current.end());
					current.clear();
				}
			}
			push_cluster();


			///////////////////////////////////////////////////////////////////////////


			// compute centroid of each clusters; 
			// centroid weighted by error
			std::vector<vec2> result;

			for (auto& cluster : clusters) {
				if (cluster.size() < cluster_min_size) continue;

				vec2 centroid = {0., 0.};
				//vec2 centroid_err = {0., 0.};
				//double total_err = 0.;

				for (size_t si : cluster) {
					double err = error_samples[si];

					centroid += samples[si].point;

					//centroid_err += samples[si].point * err;
					//total_err += err;
				}

				// weight by error
				//centroid_err /= total_err;
				centroid /= cluster.size();

				// find the nearest point of centroid on polyline 

				bool ok = true;
				for (auto& f : sdf.fonctions) {
					if ((centroid - f.point).norm() < min_distance_points) {
						ok = false;
						break;
					}
				}
				if (!ok) continue;

				result.push_back(nearest_point_on_polyline(centroid));

			}

			return result;
		}  

		void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {

			std::ofstream(output_dir + "0sdf.json") << sdf.to_json().dump(2);

			bool converged = false;

			for (size_t it = 0; it < max_it; ++it) {

				adapt_sigma();
				resolve_ls();

				std::cout << "iteration " << it << " - err : " << error_total() << " - " << sdf.fonctions.size() << " points" << std::endl; 

				std::ofstream(output_dir + std::to_string(it + 1) + "sdf.json") << sdf.to_json().dump(2);

				auto pts = where_to_add_point();
				if (pts.empty()) {
					converged = true;
					break;
				}
				for (auto& p : pts) {
					sdf.add_func({p, 0., {0., 0.}, 0.});

				}

			}

			if (!converged) {
				adapt_sigma();
				resolve_ls();
				std::cout << "iteration " << max_it << " - err : " << error_total() << " - " << sdf.fonctions.size() << " points" << std::endl; 
			}
		}
	};


}

#endif
