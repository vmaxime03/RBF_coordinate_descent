#ifndef OUTPUT_HPP__
#define OUTPUT_HPP__ 

#include "sdf.hpp"
#include "ultimaille/polyline.h"

#include <fstream>
void debug_sdf(SDF& sdf, double minx, double miny, double maxx, double maxy, const std::string& fname, int step = 100) {
	std::ofstream out(fname);
	int w = step;
	int h = step;
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			const double x = minx + (double(i) / double(w - 1)) * (maxx - minx);
			const double y = miny + (double(j) / double(h - 1)) * (maxy - miny);
			double dist = sdf.distance({x, y});
			UM::vec2 grad = sdf.gradient({x, y});
			out << x << "," << y << "," << dist << "," << grad.x << "," << grad.y ;

			if (i < w - 1) out << ",";
		}
		out << "\n";
	}
	out.close();
}

void export_polyline(UM::PolyLine& pl, const std::string& fname) {
	std::ofstream out(fname);
	for (const auto& e : pl.iter_edges()) {
		auto f = e.from().pos(), t = e.to().pos();
		out << f.x << "," << f.y << "," << t.x << "," << t.y << "\n";
	}
	out.close();
}

#endif 
