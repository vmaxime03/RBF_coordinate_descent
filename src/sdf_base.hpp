#ifndef SDF_BASE_HPP__
#define SDF_BASE_HPP__

#include "rbf.hpp"
#include "ultimaille/algebra/vec.h"
#include <concepts>
#include <memory>
#include <nlohmann/json.hpp>
#include <tuple>
#include <vector>


struct FunctionBase {
	UM::vec2 point;
	double alpha;
	UM::vec2 beta;
};

template<typename Function>
requires std::derived_from<Function, FunctionBase>
struct SDF_Base {
	using function_type = Function;

	std::vector<Function> fonctions;
	std::vector<bool> active;
	std::unique_ptr<RBF> rbf;

	virtual inline double fi(const size_t i, const UM::vec2 pos) const = 0;
	virtual inline UM::vec2 gi(const size_t i, const UM::vec2 pos) const = 0;

	explicit SDF_Base(std::unique_ptr<RBF> _rbf) : rbf(std::move(_rbf)) {};
	virtual ~SDF_Base() = default;

	virtual double distance(const UM::vec2 pos) {
		double t = 0.;
#pragma omp parallel for reduction(+:t)
		for (size_t i = 0; i < fonctions.size(); ++i) {
			if (!active[i]) continue;
			t += fi(i, pos);
		}
		return t;
	}
	virtual UM::vec2 gradient(const UM::vec2 pos) {
		double gx = 0.;
		double gy = 0.;
#pragma omp parallel for reduction(+:gx,gy)
		for (size_t i = 0; i < fonctions.size(); ++i) {
			if (!active[i]) continue;
			auto grad = gi(i, pos);	
			gx += grad.x;
			gy += grad.y;
		
		}
		return {gx, gy};
	}



	virtual inline std::tuple<double, UM::vec2> evali(const size_t i, const UM::vec2 pos) const = 0;

	virtual std::tuple<double, UM::vec2> eval(const UM::vec2 pos) {

		double d = 0.;
		double gx = 0.;
		double gy = 0.;

		#pragma omp parallel for reduction(+:d,gx,gy)
		for (size_t i = 0; i < fonctions.size(); ++i) {
			if (!active[i]) continue;
			auto [di, gi] =  evali(i, pos);
			d += di;
			gx += gi.x;
			gy += gi.y;

		}
		return {d, {gx, gy}};
	}





	inline void add_func(const Function& f) {
		fonctions.push_back(f);
		active.push_back(true);
	}

	inline void delete_point(size_t i) {
		active.erase(active.begin() + i);
		fonctions.erase(fonctions.begin() + i);
	}


	virtual std::string to_string() const = 0;
	virtual nlohmann::json to_json() const = 0;

};

#endif
