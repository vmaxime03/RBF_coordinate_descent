#ifndef RBF_HPP__
#define RBF_HPP__ 

#include <cmath>
#include <string>


struct RBF {
	virtual ~RBF() = default;
	virtual inline double f(double r, double e) const = 0;
	virtual inline double df(double r, double e) const = 0;
	virtual inline double ddf(double r, double e) const = 0;
	virtual inline std::string name() const = 0;
};


struct Gaussian : RBF {
	inline double f(double r, double s) const override {
		return std::exp(-r*r / (2*s*s));
	}
	inline double df(double r, double s) const override {
		return (-r/(s*s)) * f(r, s);
	}
	inline double ddf(double r, double s) const override {
		return ((r*r)/(s*s*s*s) - 1./(s*s)) * f(r, s);
	}
	inline std::string name() const override { return "gaussian"; }
};


struct WendlandC2 : RBF {
    inline double f(double r, double e) const override {
        double t = r / e;
        if (t >= 1.0) return 0.0;
        double s = 1.0 - t;
        return s*s*s*s * (4.0*t + 1.0);
    }

    inline double df(double r, double e) const override {
        double t = r / e;
        if (t >= 1.0) return 0.0;
        double s = 1.0 - t;
        return (s*s*s * (-4.0*(4.0*t + 1.0)) + s*s*s*s * 4.0) / e;
    }

    inline double ddf(double r, double e) const override {
        double t = r / e;
        if (t >= 1.0) return 0.0;
        double s = 1.0 - t;
        return (s*s * (20.0 * (4.0*t + 1.0) - 8.0*s*4.0 - 4.0*s*4.0)) / (e*e);
    }

    inline std::string name() const override { return "wendlandC2"; }
};

struct Pow3 : RBF {
	inline double f(double r, double e) const override {
		return r * r * r * e;
	}
	inline double df(double r, double e) const override {
		return 3 * r * r * e; 
	}
	inline double ddf(double r, double e) const override {
		return 6 * r * e;
	}

	inline std::string name() const override { return "pow3"; }
};

struct InverseMultiquadric : RBF {
	inline double f(double r, double e) const override {
		return 1. / std::sqrt(1 + e * e * r * r);
	}
	inline double df(double r, double e) const override {
		return -e * e * r / std::pow(1 + e * e * r * r, 1.5); 
	}
	inline double ddf(double r, double e) const override {
		return e * e * (2 * e * e * r * r - 1) / std::pow(1 + e * e * r * r, 2.5); 
	}

	inline std::string name() const override { return "inverse-multiquadric"; }
};




#endif
