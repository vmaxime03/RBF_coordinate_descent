#ifndef TEST_POLYLINE_HPP__
#define TEST_POLYLINE_HPP__


#include "ultimaille/algebra/vec.h"
#include "ultimaille/io/by_extension.h"
#include "ultimaille/polyline.h"
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <numbers>
#include <random>


#include "sdf_base.hpp"
#include "ultimaille/primitive_geometry.h"


#define _TEMPLATE_SDF_ \
	template<typename S, typename F>\
	requires 	std::derived_from<S, SDF_Base<typename S::function_type>>\
	&& 		std::same_as<std::invoke_result_t<F, UM::vec2, UM::vec2, double>, typename S::function_type>\



namespace SDFPointInit {

	_TEMPLATE_SDF_
		void circle(S& sdf, int npoint, const F& func, double R = 1.0, double offset = 0.) { 
			const double anglep = 2 * std::numbers::pi / npoint;

			for (int i = 0; i < npoint; ++i) {
				UM::vec2 p = {R * cos(i * anglep + offset), R * sin(i * anglep + offset)};
				sdf.add_func(func(p, -p.normalized(), 2*std::numbers::pi*R / npoint));
			}
		}

	_TEMPLATE_SDF_
	void shape(S& sdf, UM::PolyLine& pl, const F& func) {
		int npoint = pl.nedges();

		int i = 0;

		for (const auto& e : pl.iter_edges()) {
			auto v = e.to().pos() - e.from().pos();
			sdf.add_func(func((e.from().pos() + (v / 2)).xy(), UM::vec2(-v.y, v.x).normalized(), v.norm()));
			++i;
		}

	}


	// n point per edge with sigma = edge / (n + 1) * sigma_coef
	_TEMPLATE_SDF_
	void shape_multiple(S& sdf, UM::PolyLine& pl, size_t n, const F& func) {
		for (const auto& e : pl.iter_edges()) {
			UM::vec2 d = e.to().pos().xy() - e.from().pos().xy();
			auto step = d / double(n+1);
			UM::vec2 normal = UM::vec2(-d.y, d.x).normalized();
			for (size_t i = 1; i <= n; ++i) {
				UM::vec2 p = e.from().pos().xy() + i * step;
				sdf.add_func(func(p, normal, step.norm()));
			}
		}
	}


	// Equally spaced point along the polyline 

	_TEMPLATE_SDF_
	void equaly_spaced_shape(S& sdf, UM::PolyLine& pl, int n, const F& func) {

		double total = 0.0;
		for (const auto& e : pl.iter_edges()) {
			UM::Segment3 s = e;
			total += s.length();
		}

		if (n <= 0 || total <= 0.0) return;

		double step = total / double(n);
		double traveled = step/2; // offset
		double next = step + step/2;

		for (const auto& e : pl.iter_edges()) {
			UM::Segment3 s = e;
			auto v = s.xy().vector();
			double len = s.length();

			while (next <= traveled + len) {
				double t = (next - traveled) / len;
				sdf.add_func(func(
						s.a.xy() + t * v,
						UM::vec2(-v.y, v.x).normalized(),
						step
						));
				next += step;
			}

			traveled += len;
		}

	}
}

namespace PolyLineGenerator {

void regular_polygon(UM::PolyLine& pl, int nedge) {
	assert(nedge > 2);

	const double R = 1.;
	const double angle = 2 * std::numbers::pi / nedge;

	pl.points.create_points(nedge);
	pl.create_edges(nedge);



	for (auto i = 0; i < nedge; ++i) {
		UM::vec2 p = {R * cos(i * angle), R * sin(i * angle)}; // clockwise for normal pointing outward
		pl.points[i] = p.xy0();	
		pl.vert(i, 0) = i;
		pl.vert(i, 1) = (i+1)%nedge;
	}
	pl.connect();
}


void read_from_file(UM::PolyLine& pl, const std::string& fname) {
	UM::Triangles m;
	UM::read_by_extension(fname, m);
	m.connect();
	std::vector<bool> visited(m.nfacets()*3, false);

	double minx = DBL_MAX;
	double miny = DBL_MAX;
	double maxx = -DBL_MAX;
	double maxy = -DBL_MAX;


	for (const auto& che : m.iter_halfedges()) {

		if (!che.on_boundary() || visited[che]) continue;

		UM::Surface::Halfedge he = che;
		std::vector<UM::Surface::Vertex> curr_pl;

		do {
			visited[he] = true;

			curr_pl.push_back(he.from());
			UM::Surface::Halfedge nhe = he.next();

			while (!nhe.on_boundary()) {
				nhe = nhe.opposite().next();
			}
			he = nhe;
		} while (he != che);

		int start = pl.points.size();

		pl.points.create_points(curr_pl.size());
		pl.create_edges(curr_pl.size());

		for (uint i = 0; i < curr_pl.size(); ++i) {
			auto p = m.points[curr_pl[i]];

			minx = std::min(minx, p.x);
			miny = std::min(miny, p.y);
			maxx = std::max(maxx, p.x);
			maxy = std::max(maxy, p.y);

			pl.points[start + i] = p;
			pl.vert(start + i, 0) = start + i;
			pl.vert(start + i, 1) = start + ((i+1) % curr_pl.size());
		}
	}

	double range = std::max(maxx - minx, maxy - miny);
	double b = 2.0;
	double cx = (maxx + minx) / 2.0;
	double cy = (maxy + miny) / 2.0;

	for (auto& p : pl.points) {
		p.x = ((p.x - cx) / range) * b;
		p.y = ((p.y - cy) / range) * b;
	}

	pl.connect();

}


void isoceles_triangle(UM::PolyLine& pl, double angle) {
	assert (angle < std::numbers::pi && angle > 0);

	pl.points.create_points(3);

	double hsin = std::sin(angle/2);
	double hcos = std::cos(angle/2);
	


	pl.points[0] = UM::vec3{-hsin, -hcos/2, 0.}.normalized();
	pl.points[1] = UM::vec3{hsin, -hcos/2, 0.}.normalized();
	pl.points[2] = UM::vec3{0., hcos/2, 0.}.normalized();

	pl.create_edges(3);

	pl.vert(0, 0) = 0;
	pl.vert(0, 1) = 1;
	pl.vert(1, 0) = 1;
	pl.vert(1, 1) = 2;
	pl.vert(2, 0) = 2;
	pl.vert(2, 1) = 0;

	pl.connect();

}


// claude
void random_polygon(UM::PolyLine& pl, int nedge, unsigned int seed = 42) {
    assert(nedge > 2);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Generate random points, then sort by angle around centroid
    // — sorting by angle guarantees a simple polygon
    std::vector<UM::vec2> pts(nedge);
    for (auto& p : pts)
        p = {dist(rng), dist(rng)};

    // Compute centroid
    UM::vec2 C = {0, 0};
    for (auto& p : pts) C = C + p;
    C = C * (1.0 / nedge);

    // Sort by angle around centroid → simple convex-ish polygon
    std::sort(pts.begin(), pts.end(), [&](const UM::vec2& a, const UM::vec2& b) {
        return std::atan2(a.y - C.y, a.x - C.x) < std::atan2(b.y - C.y, b.x - C.x);
    });

    pl.points.create_points(nedge);
    pl.create_edges(nedge);

    for (int i = 0; i < nedge; ++i) {
        pl.points[i] = pts[i].xy0();
        pl.vert(i, 0) = i;
        pl.vert(i, 1) = (i + 1) % nedge;
    }


    pl.connect();
}


}


#endif
