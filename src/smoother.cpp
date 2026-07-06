#include <algorithm>
#include <chrono>
#include <print>
#include "ultimaille/all.h"


#include "algo_base.hpp"
#include "debug_macros.hpp"
#include "sdf_elliptic.hpp"
#include "ultimaille/polyline.h"
#include <utility>
#include "rbf.hpp"
#include "sdf.hpp"
#include "sdf_elliptic.hpp"
#include "samples.hpp"
#include "algo_base_elliptic.hpp"
#include "algo_ellipse_imp.hpp"
#include "test_init.hpp"
#include <cfloat>
#include <csignal>
#include <cstddef>
#include <filesystem>
#include <numbers>
#include <utility>
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"
#include "algo_imp1.hpp"
#include "algo_imp2.hpp"
#include "output.hpp"
#include "samples.hpp"
#include "algo_base_elliptic.hpp"
#include "algo_ellipse_imp.hpp"

#include "nlohmann/json.hpp"

using namespace UM;
using Clock = std::chrono::high_resolution_clock;
using namespace std::literals::chrono_literals;

inline double chi(double eps, double det) {
    if (det>0)
        return (det + std::sqrt(eps*eps + det*det))*.5;
    return .5*eps*eps / (std::sqrt(eps*eps + det*det) - det);
}
inline double chi_deriv(double eps, double det) {
    return .5+det/(2.*std::sqrt(eps*eps + det*det));
}

int main(int argc, char** argv) {
	
	const std::string input_dir = INPUT_DIR;

	std::string plfilename = "bord.obj";
	std::string quadfilename = "quad_mesh.obj";

	int polyline_segmant_nsample = 5;
	double target_function_width = 0.075;
	int target_interpolant_neighbors = 6;
	double lambda_distance = 1.;
	double smoother_w = 3.0;
	bool export_field = true;
	int neighborhood_size = target_interpolant_neighbors * 2;

	if (argc > 1) {
		std::ifstream config_file(argv[1]);

		if (!config_file) {
			std::println(stderr, "Error: Could not open or find configuration file '{}'", argv[1]);
			return 1;
		}

		nlohmann::json config;
		try {
			config_file >> config;
		} catch (const std::exception& e) {
			std::println(stderr, "Error: invalid JSON in '{}': {}", argv[1], e.what());
			return 1;
		}
		plfilename                    = config.value("bord", plfilename);
		quadfilename                  = config.value("quad", quadfilename);
		polyline_segmant_nsample      = config.value("nsample", polyline_segmant_nsample);
		target_function_width         = config.value("target_function_width", target_function_width);
		target_interpolant_neighbors  = config.value("target_interpolant_neighbors", target_interpolant_neighbors);
		lambda_distance               = config.value("lambda_distance", lambda_distance);
		smoother_w                    = config.value("smoother_w", smoother_w);
		export_field                  = config.value("export_field", export_field);
		std::println("Loaded configuration from '{}'", argv[1]);
	}
    std::println(
		"Parameters :\n"
        "  bord                         = {}\n"
        "  quads                        = {}\n"
        "  nsample                      = {}\n"
        "  target_function_width        = {}\n"
        "  target_interpolant_neighbors = {}\n"
        "  lambda_distance              = {}\n"
        "  smoother_w                   = {}\n"
        "  export_field                 = {}\n",
        plfilename, quadfilename, polyline_segmant_nsample, target_function_width, target_interpolant_neighbors, lambda_distance, smoother_w, export_field);


	std::string plfile = INPUT_DIR + plfilename;
	std::string quadfile = INPUT_DIR + quadfilename;


	const std::string output_dir = OUTPUT_DIR;
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	std::filesystem::create_directories(output_dir + "smoother/");
	std::filesystem::create_directories(output_dir + "test/");

	// SDF INIT
	auto rbf = std::make_unique<WendlandC2>();
	UDF_Elliptic sdf(std::move(rbf));
	PolyLine pl;
	read_by_extension(plfile, pl);
	pl.connect();

	write_by_extension(output_dir + "smoother/bord.obj", pl);


	

	auto NORMAL = [](const vec2& v) -> vec2 { return {-v.y, v.x}; };
	SDFPointInit::shape(sdf, pl, [&](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, NORMAL(n) * r*0.5, r*0.5*n};});
	auto samples = samples::compute_edges_samples_normals(pl, polyline_segmant_nsample);
	auto algo = sdffitting_elliptic::TODONAMEFitter(sdf, samples);
	algo.fit_ellipses_radius(target_function_width, target_interpolant_neighbors);

	TIME_DEBUG("fit ellipse");
	algo.lambda_distance = lambda_distance;
	algo.resolve_ls();
	sdf.optimize();
	TIME_DEBUG("sdf optimization");

	sdf.neighborhood_size = neighborhood_size;

	if (export_field) {
		// output::export_samples(samples, output_dir + "samples.csv");
		output::export_polyline(pl, output_dir + "test/polyline.csv");
		output::export_sdf(sdf, output_dir + "test/sdf_params.csv");
		double samples_offset = target_function_width; 
		output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "test/sdf.csv", 1000);
		// std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;
		TIME_DEBUG("export udf");

	}

	// LISSEUR
    Quads m;
    read_by_extension(quadfile, m);
    m.connect();

    constexpr double theta = .25; // the energy is (1-theta)*(shape energy) + theta*(area energy)
    constexpr int bfgs_maxiter  = 100000;
    constexpr double bfgs_threshold  = 1e-14;

    CornerAttribute<mat<3,2>> reference(m);  // desired tri geometry
    CornerAttribute<double> area(m);
    std::vector<double> X(m.nverts()*2, 0.);
    for (int v : m.iter_vertices())
        for (int d : {0,1})
            X[2*v+d] = m.points[v][d];

    double total_area = 0;
    for (auto f : m.iter_facets())
        total_area += Quad3(f).unsigned_area();

    for (auto h : m.iter_halfedges()) {
        double e = std::sqrt(total_area/m.nfacets());
        vec2 A = {0,0}, B = {e, 0},C = {0, e};
        mat<2,2> ST = {{B-A, C-A}};
        reference[h] = mat<3,2>{{ {-1,-1},{1,0},{0,1} }}*ST.invert_transpose();
        area[h] = e*e;
    }

    const auto getJ = [&reference](const std::vector<double>& X, Surface::Halfedge h)->mat<2, 2> {
        int verts[3] = {h.from(), h.to(), h.prev().from() };
        mat<2, 2> J = {};
        for (int i : {0, 1, 2})
            for (int d : {0, 1})
                J[d] += reference[h][i] * X[verts[i]*2 + d];
        return J;
    };

    auto starting_time = Clock::now();
    std::vector<SpinLock> spin_locks(X.size());

    double mindet = 0.;
    int ninverted = 0;
    for (auto h : m.iter_halfedges()) {
        const mat<2, 2> J = getJ(X, h);
        double det = J.det();
        mindet = std::min(mindet, det);
        ninverted += (det<=0);
    }

    std::println("mindet: {}, inverted corners: {}", mindet, ninverted);

    constexpr double e0 = 1e-4;
    double eps = mindet>0 ? e0 : std::sqrt(e0*e0 + 0.004*mindet*mindet);
    std::cerr << "eps: " << eps << std::endl;

    const STLBFGS::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;
        ninverted = 0;
        mindet = std::numeric_limits<double>::max();
#pragma omp parallel for reduction(+:F) reduction(min:mindet) reduction(+:ninverted)
        for (int c=0; c<m.ncorners(); c++) {
            Surface::Halfedge h = m.halfedge(c);
            const mat<2, 2> J = getJ(X, m.halfedge(h));
            const mat<2, 2> K = { {{ +J[1].y, -J[1].x }, { -J[0].y, +J[0].x }} };

            const double det = J[0]*K[0];

            mindet = std::min(mindet, det);
            ninverted += (det<=0);

            const double c1 = chi(eps, det);
            const double c2 = chi_deriv(eps, det);

            const double f = J.sumsqr()/(2.*c1);
            const double g = (1+det*det)/(2.*c1);
            F += ((1-theta)*f + theta*g) * area[h];

            for (int d : {0, 1}) {
                const vec2& a = J[d];
                const vec2& b = K[d];
                const vec2 dfda = (a - b*f*c2)/c1;
                const vec2 dgda = b*(det - g*c2)/c1;
                for (int i : {0, 1, 2}) {
                    int verts[3] = {m.halfedge(h).from(), m.halfedge(h).to(), m.halfedge(h).prev().from() };
                    const int v = verts[i];
                    spin_locks[v*2+d].lock();
                    G[v*2+d] += (dfda*(1.-theta) + dgda*theta) * area[h] * reference[h][i];
                    spin_locks[v*2+d].unlock();
                }
            }

#if 1
            constexpr double w = 5e+3;
            int v[3] = { h.from(), h.to(), h.prev().from() };
            vec2 a = vec2{X[v[1]*2+0], X[v[1]*2+1]} - vec2{X[v[0]*2+0], X[v[0]*2+1]};
            vec2 b = vec2{X[v[2]*2+0], X[v[2]*2+1]} - vec2{X[v[0]*2+0], X[v[0]*2+1]};
            double s = a*b;
            F += w*s*s;

            for (int d : {0, 1}) {
                    spin_locks[v[1]*2+d].lock();
                    G[v[1]*2+d] += w*2*s*b[d];
                    spin_locks[v[1]*2+d].unlock();
                    spin_locks[v[2]*2+d].lock();
                    G[v[2]*2+d] += w*2*s*a[d];
                    spin_locks[v[2]*2+d].unlock();
                    spin_locks[v[0]*2+d].lock();
                    G[v[0]*2+d] += w*-2*s*(a+b)[d];
                    spin_locks[v[0]*2+d].unlock();
            }
#endif
        }

#pragma omp parallel for reduction(+:F)
        for (int c=0; c<m.ncorners(); c++) {

            Surface::Halfedge h = m.halfedge(c);
            if (h.opposite().active()) continue;

            vec2 a = {X[2*h.from()], X[2*h.from()+1]};
            vec2 b = {X[2*h.to()  ], X[2*h.to()  +1]};

            double w = smoother_w; // TODO
            for (double t=0; t<=1; t+=1e-2) {
                vec2 p = a*(1-t) + t*b;

                auto [val, grd] = sdf.eval(p); 
				
                F += w * val;

                for (int d : {0, 1}) {
                    spin_locks[h.from()*2+d].lock();
                    G[h.from()*2+d] += w * grd[d] * (1-t);
                    spin_locks[h.from()*2+d].unlock();

                    spin_locks[h.to()*2+d].lock();
                    G[h.to()*2+d] += w * grd[d] * t;
                    spin_locks[h.to()*2+d].unlock();
                }
            }
        }
    };

    double E_prev, E;
    std::vector<double> trash(X.size());
    func(X, E_prev, trash);

    STLBFGS::Optimizer opt{func};
    opt.ftol = opt.gtol = bfgs_threshold;
    opt.maxiter = bfgs_maxiter;
    opt.run(X);

    func(X, E, trash);
    std::println("E: {} --> {}, eps: {}, min det: {}", E_prev, E, eps, mindet);
    std::println("Running time: {} seconds", (Clock::now()-starting_time)/1.s);

    for (int v : m.iter_vertices())
        m.points[v] = {X[2*v+0], X[2*v+1]};
    write_by_extension(output_dir + "smoother/smoothed.obj", m);
    return 0;
}

