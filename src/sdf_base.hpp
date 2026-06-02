#ifndef SDF_BASE_HPP__
#define SDF_BASE_HPP__

#include "rbf.hpp"
#include "ultimaille/algebra/vec.h"
#include <memory>
#include <nlohmann/json.hpp>
#include <vector>


template<typename Function>
struct SDF_Base {
	std::vector<Function> fonctions;
	std::vector<bool> active;
	std::unique_ptr<RBF> rbf;

	virtual inline double fi(size_t i, UM::vec2 pos) const = 0;
	virtual inline UM::vec2 gi(size_t i, UM::vec2 pos) const = 0;

	explicit SDF_Base(std::unique_ptr<RBF> _rbf) : rbf(std::move(_rbf)) {};
	virtual ~SDF_Base() = default;

	double distance(UM::vec2 pos) {
		double t = 0.;
		for (size_t i = 0; i < fonctions.size(); ++i) {
			if (!active[i]) continue;
			t += fi(i, pos);
		}
		return t;
	}
	UM::vec2 gradient(UM::vec2 pos) {
		UM::vec2 grad(0., 0.);
		for (size_t i = 0; i < fonctions.size(); ++i) {
			if (!active[i]) continue;
			grad += gi(i, pos);	
		
		}
		return grad;
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
