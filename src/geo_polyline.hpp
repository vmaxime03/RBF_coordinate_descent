#ifndef GEO_POLYLINE_HPP__
#define GEO_POLYLINE_HPP__

#include "ultimaille/polyline.h"
#include <cfloat>

namespace UM {

	// distance min point - polyline non signee
	double unsigned_distance_point_polyline_2d(PolyLine& p, const vec2& pos) {

		double mind = DBL_MAX;
		for (const auto& edge : p.iter_edges()) {
			auto p1 = edge.from().pos();
			auto p2 = edge.to().pos();

			mind = std::min(mind, Segment2(p1.xy(), p2.xy()).distance(pos));
		}
		return mind;

	}
	int crossing_number_polyline_2d(PolyLine& p, vec2 pos) {
		int cn = 0;    
		for (auto e : p.iter_edges()) {   
			vec2 a = e.from().pos().xy();
			vec2 b = e.to().pos().xy();

			cn += (((a.y <= pos.y) && (b.y > pos.y)) || ((a.y > pos.y) && (b.y <=  pos.y))) // upward/downward crossing
				&& (pos.x < a.x + ((pos.y - a.y)/ (b.y - a.y)) * (b.x - a.x)); // crossing a droite de pos

		}
		return (cn & 1); // 0 out 1 in
	}
	int winding_number_polyline_2d(PolyLine& p, vec2 pos) { 	
		int wn = 0;
		for (auto e : p.iter_edges()) {
			vec2 a = e.from().pos().xy();
			vec2 b = e.to().pos().xy();

			double orientation = cross((b - a).xy0(), (pos - a).xy0()).z;
			wn += 1 * ((a.y <= pos.y) && (b.y > pos.y) && (orientation > 0)); // upward crossing
			wn -= 1 * ((a.y > pos.y) && (b.y <= pos.y) && (orientation < 0)); // downward crossing
		}
		return wn;
	}

	double signed_distance_point_polyline_2d(PolyLine& p, vec2 pos) {
		double mind = DBL_MAX;
		int wn = 0;
		for (auto edge : p.iter_edges()) {
			auto a = edge.from().pos().xy();
			auto b = edge.to().pos().xy();

			double d = Segment2(a, b).distance(pos);

			mind = std::min(d, mind);

			double orientation = cross((b - a).xy0(), (pos - a).xy0()).z;
			wn += 1 * ((a.y <= pos.y) && (b.y > pos.y) && (orientation > 0));
			wn -= 1 * ((a.y > pos.y) && (b.y <= pos.y) && (orientation < 0));

		}
		return ((wn == 0)?1:-1) * mind;
	}

	// aire signée d'une polyline fermee https://cs108.epfl.ch/archive/15/p02_geometry.html
	double closed_polyline_signed_area(PolyLine& pl) {
		double r = 0.;
		for (const auto& e : pl.iter_edges()) {
			r += (e.from().pos().x * e.to().pos().y - e.to().pos().x * e.from().pos().y);
		}
		return 0.5 * r;
	}

}

#endif
