#ifndef SAMPLES_HPP__
#define SAMPLES_HPP__

#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"

namespace samples {

	struct PointNormal {
		UM::vec2 point; 
		UM::vec2 normal;   
		size_t   loop_id; 
		size_t   next;     
		size_t   prev;     

		PointNormal(UM::vec2 p, UM::vec2 n, size_t lid)
			: point(p), normal(n), loop_id(lid), next(-1), prev(-1) {}
	};

	typedef std::vector<PointNormal> Samples;


	Samples compute_edges_samples_normals(UM::PolyLine& pl, size_t n = 10, bool flip_normals = false) {
		Samples samples;

		size_t loop_count = 0;
		auto prev_e = pl.iter_edges().begin().h;

		for (const auto& e : pl.iter_edges()) {

			// if previous edge closed a loop, start a new one
			if (prev_e.to() != e.from()) {
				++loop_count;
			}

			auto d = UM::Segment3(e).vector().xy();
			auto step = d / double(n+1);
			auto normal = UM::vec2(-d.y, d.x).normalized();

			for (size_t i = 1; i <= n; ++i) {
				auto p = e.from().pos().xy() + i * step;

				samples.push_back({
						p,
						(flip_normals ? -1 : 1) * normal,
						loop_count
						});
			}

			prev_e = e;
		}

		return samples;

	
	}


	Samples compute_edges_samples_normals_include_end(UM::PolyLine& pl, size_t n = 10, bool flip_normals = false) {
		Samples samples;

		size_t loop_count = 0;
		auto prev_e = pl.iter_edges().begin().h;

		for (const auto& e : pl.iter_edges()) {

			// if previous edge closed a loop, start a new one
			if (prev_e.to() != e.from()) {
				++loop_count;
			}

			auto d = UM::Segment3(e).vector().xy();
			auto step = d / double(n-1);
			auto normal = UM::vec2(-d.y, d.x).normalized();

			for (size_t i = 0; i < n; ++i) {
				auto p = e.from().pos().xy() + i * step;

				samples.push_back({
						p,
						(flip_normals ? -1 : 1) * normal,
						loop_count
						});

				
			}

			prev_e = e;
		}

		return samples;

	
	}


	Samples compute_equally_spaced_samples_normals(UM::PolyLine& pl, int n, bool flip_normals = false) {
		Samples samples;

		size_t loop_count = 0;
		auto prev_e = pl.iter_edges().begin().h;



		double total = 0.0;
		for (const auto& e : pl.iter_edges()) {
			UM::Segment3 s = e;
			total += s.length();
		}
		double step = total / double(n);
		double traveled = step / 2;
		double next = step + step / 2;

		
		for (const auto& e : pl.iter_edges()) {


			// if previous edge closed a loop, start a new one
			if (prev_e.to() != e.from()) {
				++loop_count;
			}


			UM::Segment3 s = e;
			auto v = s.xy().vector();
			double len = s.length();
			while (next <= traveled + len) {
				double t = (next - traveled) / len;
				samples.push_back({
						s.a.xy() + t * v,
						(flip_normals ? -1 : 1) * UM::vec2(-v.y, v.x).normalized(),
						loop_count
						});
				next += step;
			}
			traveled += len;


			prev_e = e;
		
		}

		return samples;
	}

} // namespace samples

#endif
