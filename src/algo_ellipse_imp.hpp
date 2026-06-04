#ifndef ALGO_ELLIPSE_IMP__
#define ALGO_ELLIPSE_IMP__

#include <algorithm>
#include <cfloat>
#include <numbers>
#include <sstream>
#include "algo_base_elliptic.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/knn.h"

#include "debug_macros.hpp"



namespace sdffitting_elliptic {


struct TODONAMEFitter : Fitter_Elliptic {
	using Fitter_Elliptic::Fitter_Elliptic;


	void adapt_minor() {

		std::vector<vec2> pts(sdf.fonctions.size());
		std::transform(sdf.fonctions.begin(), sdf.fonctions.end(), pts.begin(), [](const auto& f) {
				return f.point;
				});

		// auto knn = UM::KNN<2>(pts);
		
		for (size_t i = 0; i < sdf.fonctions.size(); ++i) {

			auto& f = sdf.fonctions[i];

			// TODO find better way to find function we dont want to interact with
			// auto neigh = knn.query(f.point, 2);

			// for (const auto& j : neigh) {
			for (size_t j = 0; j < sdf.fonctions.size(); ++j) {
				if (j == i) continue;
		
				// find the maximal minor radius that avoid interaction with function that have opposites normals
				// assertion : minor radius have same direction as normals 

				auto d = sdf.fonctions[j].point - f.point;

				auto n1 = f.ellipse_minor;
				auto n2 = sdf.fonctions[j].ellipse_minor;

				auto n1_norm = n1.normalized();
				auto n1_max = n1_norm * (f.ellipse_major.norm());
				auto n2_norm = n2.normalized();

				double cos = n1_norm * n2_norm;

				auto hd = d/2.0;

				auto m = std::abs(hd * n1_norm) * n1_norm;
				
				auto r = n1_max - m;

				double scale = (cos + 1.0) / 2.0;

				auto new_minor = m + r * scale;
				if (new_minor.norm() < f.ellipse_minor.norm()) f.ellipse_minor = new_minor;

			}
		}
		}



		void decimate(double angle_threshold) {

			std::vector<bool> visited(sdf.fonctions.size(), false);

			std::vector<std::vector<size_t>> clusters;
			struct triplet { vec2 centroid, major, minor; };
			std::vector<triplet> clusters_centroid_major_minor;


			std::function<void(size_t)> expand_cluster = [&](size_t fi) -> void {
				auto& f = sdf.fonctions[fi];

				// TODO replace with finding compatible cluster
				size_t cid = clusters.size();
				clusters.push_back({fi});
				clusters_centroid_major_minor.push_back({f.point, f.ellipse_major, f.ellipse_minor});

				for (size_t i = 0; i < sdf.fonctions.size(); ++i) {
					if (visited[i]) continue;
					visited[i] = true;

					auto& n = sdf.fonctions[i];


					// compare n with cluster centroid
					
					// check angle compatibility
					double angle = n.ellipse_minor.normalized() * clusters_centroid_major_minor[cid].minor.normalized();
					if (angle < std::cos(angle_threshold)) {
						expand_cluster(i);
						return;
					}

					// check distance 
					
					auto center_dist = (f.point - n.point);

					auto d_major = std::abs(center_dist * f.ellipse_major.normalized());
					auto d_minor = std::abs(center_dist * f.ellipse_minor.normalized());

					// check if both bounding rectangle collide
					if (!(d_major <= (f.ellipse_major.norm() + n.ellipse_major.norm()) && d_minor <= (f.ellipse_minor.norm() + n.ellipse_minor.norm()))) {
						expand_cluster(i);
						return;
					}


					// add fi to cluster
					clusters[cid].push_back(i);

					// update cluster centroid, major and minor
					vec2 centroid = {0, 0};
					vec2 minor = {0, 0};
					for (size_t j : clusters[cid]) {
						auto& fj = sdf.fonctions[j];

						centroid += fj.point;
						minor += fj.ellipse_minor;
					}
					centroid /= clusters[cid].size();
					minor /= clusters[cid].size();
				
					double major_radius = -DBL_MAX;
					for (size_t j : clusters[cid]) {
						auto& fj = sdf.fonctions[j];
						
						auto p1 = fj.point + fj.ellipse_major;
						auto p2 = fj.point - fj.ellipse_major;

						major_radius = std::max(major_radius, (centroid - p1).norm());
						major_radius = std::max(major_radius, (centroid - p2).norm());
					}
					auto major = vec2(-minor.y, minor.x).normalized() * major_radius;

					clusters_centroid_major_minor[cid] = {centroid, major, minor};

				}
				
			};


			visited[0] = true;
			expand_cluster(0);


			DEBUG(sdf.fonctions.size());
			DEBUG(clusters.size());

			sdf.fonctions.resize(clusters.size());
			sdf.active.resize(clusters.size());

			for (size_t i = 0; i < clusters.size(); ++i) {
				auto& c = clusters_centroid_major_minor[i];
				sdf.fonctions[i] = {{c.centroid, 0., {0., 0.}}, c.major, c.minor};
				sdf.active[i] = true;
			}

		}



	void fit(size_t max_it, size_t snapshot, const std::string &output_dir) override {

	}

};


}




#endif 
