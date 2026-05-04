#ifndef RBF_HPP__
#define RBF_HPP__ 

#include <cmath>
#include <string>


struct RBF {
	virtual ~RBF() = default;
	virtual inline double f(double r, double e) const = 0;
	virtual inline double df(double r, double e) const = 0;
	virtual inline std::string name() const = 0;
};

/*
struct WendlandC2 : RBF { // TODO 
	inline double f(double r, double p) const override {
		return 0.;
	}

	inline double df(double r, double p) const override {
		return 0.;
	}
};
*/


struct Gaussian : RBF {
	inline double f(double r, double s) const override {
		return std::exp(-r*r / (2*s*s));
	}
	inline double df(double r, double s) const override {
		return (-r/(s*s)) * f(r, s);
	}
	inline std::string name() const override { return "gaussian"; }
};


#endif
