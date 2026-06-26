#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include "ultimaille/primitive_geometry.h"
#include "ultimaille/sparse/least_squares.h"
#include "ultimaille/sparse/linexpr.h"
#include <numbers>


void test_rbf_ellipse(UM::PolyLine& pl) {

	struct Cluster {
		std::vector<UM::vec2> points;

		UM::vec2 ellipse_dir; // normal to the cluster => opposite of direction with circle center
		UM::vec2 ellipse_center; // center point


		double ellipse_A; // big radius = cluster space between start and end point
		double ellipse_B; // small radius


		double ellipse_R; // circular radius if circle fit, -1 if line fit

	};

	std::vector<Cluster> clusters;
	
	


	double LINE_THRESHOLD = 1e-6; // MAX LINE ERROR
	double CIRCLE_MAX_ANGLE = std::numbers::pi/30; // MAX CIRCLE CLUSTER ANGLE DIFFERENCE
	double CIRCLE_THRESHOLD = 1e-3; // MAX CIRCLE ERROR
	double CIRCLE_MIN_RADIUS = 1.; // MIN CIRCLE RADIUS TODO use


	// INIT
	auto it = pl.iter_edges().begin();
	auto& last_edge = it.h;
	++it;



	std::vector<UM::PolyLine::Edge> current_cluster;



	while (it != pl.iter_edges().end()) {
		auto& edge = it.h;

		UM::Segment3 last_edge_segment = last_edge;
		UM::Segment3 edge_segment = edge;
		
		auto cross = UM::cross(last_edge_segment.vector(), edge_segment.vector()).z;

		if (cross < LINE_THRESHOLD) { // COLINEAR => LINE CLUSTER
			if (current_cluster.empty()) {
				current_cluster.push_back(std::move(edge));
				// TODO append to cluster
			} else {
				// TODO new cluster
			}
		}



		if (last_edge_segment.vector().xy().normalized() * edge_segment.vector().xy().normalized() < std::cos(CIRCLE_MAX_ANGLE)) {
			// TODO new cluster, angle toosharp
		}


		// FIT CIRCLE 
		// TODO possibly skip if cluster is size 1 or 2 (<= 3 points)

		// we fit a circumcircle to the clusters points
		
		UM::LeastSquares ls(3);

		auto add_point_to_ls = [&](const UM::vec3& p) {
			UM::LinExpr expr = p.x*p.x + p.y*p.y - UM::Linear::X(0) * p.x - UM::Linear::X(1) * p.y - UM::Linear::X(2);
			ls.add_to_energy(expr);
		};

		bool first = true;
		for (const auto& cluster_edge : current_cluster) {
			if (first)  {
				add_point_to_ls(cluster_edge.from().pos());
				first = false;
			}
			add_point_to_ls(cluster_edge.to().pos());
		}

		add_point_to_ls(edge.to().pos());

		ls.solve();


		double cx = ls.value(0)/2;
		double cy = ls.value(1)/2;
		double R = std::sqrt(cx*cx + cy*cy + ls.value(2));


		if (R < CIRCLE_MIN_RADIUS) {
			// TODO new cluster
		}


		auto error_dist_to_circle = [&](const UM::vec3& p) {
			double dx = p.x - cx;
			double dy = p.y - cy;
			return std::sqrt(dx*dx + dy*dy) - R;
		};

		bool is_new_circle_good = true;

		first = true;
		for (const auto& cluster_edge : current_cluster) {
			if (first)  {
				if (error_dist_to_circle(cluster_edge.from().pos()) > CIRCLE_THRESHOLD) {
					is_new_circle_good = false;
					break;
				}
			}
			if (error_dist_to_circle(cluster_edge.to().pos()) > CIRCLE_THRESHOLD) {
				is_new_circle_good = false;
				break;
			}
		}

		if (is_new_circle_good) {
			// TODO add cluster
		} else {
			// TODO new cluster
		}

		++ it;
	}



}
