#ifndef SAMPLES_HPP__
#define SAMPLES_HPP__

#include "types.hpp"
#include "ultimaille/polyline.h"

namespace samples {

	Samples compute_edges_samples_normals(UM::PolyLine& pl, size_t n = 10, bool flip_normals = false) {
		Samples samples;
		for (const auto& e : pl.iter_edges()) {
			UM::vec2 d = e.to().pos().xy() - e.from().pos().xy();
			auto step = d / double(n+1);
			UM::vec2 normal = UM::vec2(-d.y, d.x).normalized();
			for (size_t i = 1; i <= n; ++i) {
				UM::vec2 p = e.from().pos().xy() + i * step;
				samples.push_back({p, (flip_normals ? -1 : 1) * normal});
			}
		}
		return samples;
	}

	Samples compute_equally_spaced_samples_normals(UM::PolyLine& pl, int n, bool flip_normals = false) {
		Samples samples;
		double total = 0.0;
		for (const auto& e : pl.iter_edges()) {
			UM::Segment3 s = e;
			total += s.length();
		}
		double step = total / double(n);
		double traveled = step / 2;
		double next = step + step / 2;

		for (const auto& e : pl.iter_edges()) {
			UM::Segment3 s = e;
			auto v = s.xy().vector();
			double len = s.length();
			while (next <= traveled + len) {
				double t = (next - traveled) / len;
				samples.push_back({
						s.a.xy() + t * v,
						(flip_normals ? -1 : 1) * UM::vec2(-v.y, v.x).normalized()
						});
				next += step;
			}
			traveled += len;
		}
		return samples;
	}

} // namespace samples

#endif
