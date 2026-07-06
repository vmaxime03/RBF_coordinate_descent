#ifndef SDF_SMOOTHMIN_HPP
#define SDF_SMOOTHMIN_HPP


#include "sdf_base.hpp"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/polyline.h"
#include "ultimaille/primitive_geometry.h"
#include <cfloat>
#include <functional>


struct SMIN {
	virtual double smin(double a, double b, double k) const = 0;
	virtual UM::vec2 dsmin(double a, double b, double k) const = 0;
};


struct EXPSMIN : SMIN {
double smin(double a, double b, double k) const override {
	return -k * std::log2( std::exp2(-a/k) + std::exp2(-b/k) );
}
UM::vec2 dsmin(double a, double b, double k) const override {
	double va = std::exp2(-a/k);
	double vb = std::exp2(-b/k);
	return  {va / (va+vb), vb / (va+vb)};
}
};


struct RFUNC_MIN : SMIN {
	double smin(double a, double b, double k) const override {
		assert(k > -1 && k <= 1);
		if (a >= DBL_MAX - 1.0) return b;
		if (b >= DBL_MAX - 1.0) return a;

		return a + b - std::sqrt(a*a + b*b);

	}

	UM::vec2 dsmin(double a, double b, double k) const override {
		assert(k > -1 && k <= 1);
		if (a >= DBL_MAX - 1.0) return {0.0, 1.0};
		if (b >= DBL_MAX - 1.0) return {1.0, 0.0};

		double sqrt = std::sqrt(a*a + b*b);
		return {1 - a/sqrt, 1 - b/sqrt};
	}
};

struct SMOOTHMIN : SMIN {
	typedef std::function<double(double)> Kernel;

	Kernel g;
	Kernel dg;
	double g0;

	SMOOTHMIN(Kernel kernel, Kernel kernel_derivative) : g(kernel), dg(kernel_derivative), g0(kernel(0.0)) {}

	double smin(double a, double b, double k) const override {
		double l = k/g0;
		double z = (b-a)/l;
		return b - l * g(z);
	}

	UM::vec2 dsmin(double a, double b, double k) const override {
		double l = k/g0;
		double z = (b-a)/l;
		return {dg(z), 1-dg(z)};
	}

};


struct SDF_Smoothmin  {

    UM::PolyLine& polyline;       
	SMIN& smin;                     

	explicit SDF_Smoothmin(UM::PolyLine& pl, SMIN& _smin) : polyline(pl), smin(_smin) {}

	virtual ~SDF_Smoothmin() = default;

	double k = 1e-2;
	
	// TODO use KNN struct to speed up
	double distance(UM::vec2 pos) const {

		double distance = DBL_MAX;

		for (const auto& e : polyline.iter_edges()) {
			UM::Segment2 s = {e.from().pos().xy(), e.to().pos().xy()};

			distance = smin.smin(s.distance(pos), distance, k);
		}
		return distance;
	}

	UM::vec2 gradient(UM::vec2 pos) const {
		double distance = DBL_MAX;
		UM::vec2 gradient = {0, 0};



		for (const auto& e : polyline.iter_edges()) {
			UM::Segment2 s = {e.from().pos().xy(), e.to().pos().xy()};
			double di = s.distance(pos);

			// gradient direction
			UM::vec2 gdi = UM::vec2(-s.vector().y, s.vector().x).normalized();


			UM::vec2 gmi = smin.dsmin(di, distance, k);


			gradient = gdi * gmi.x + gradient * gmi.y;
			distance = smin.smin(distance, di, k);

		}

		return gradient;
	}


};


#endif // !SDF_SMOOTHMIN_HPP

