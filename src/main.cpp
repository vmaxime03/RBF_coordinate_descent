#include "algo_base.hpp"
#include "debug_macros.hpp"
#include "sdf_elliptic.hpp"
#include "ultimaille/polyline.h"
#include <algorithm>
#include <cfloat>
#include <csignal>
#include <cstddef>
#include <filesystem>
#include <memory>
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

#include "sdf_smoothmin.hpp"


// #include "test_eigen.hpp"

using namespace UM;

int main(int argc, char** argv) {

	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;

	TIME_DEBUG_INIT();

	PolyLine polyline;
	auto rbf = std::make_unique<WendlandC2>();
	
	
	// int n = 8; PolyLineGenerator::regular_polygon(pl, n);
	// PolyLineGenerator::read_from_file(polyline, input_dir + "duck.geogram");

	// PolyLineGenerator::read_from_file(pl, input_dir + "lapinpluslisseplusdense_quad_mesh.obj");
	
	// PolyLineGenerator::read_from_file(pl, input_dir + "cerf.geogram");
	// PolyLineGenerator::read_from_file(polyline, input_dir + "o.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o1.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o2.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "ooo.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "ooo2.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "e.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");

	// PolyLineGenerator::isoceles_triangle(pl, std::numbers::pi / 4);
	

	PolyLineGenerator::read_polyline_directly(polyline, input_dir + "lapinpluslisseplusdense_bord.obj");


	// SDF sdf(std::move(rbf));

	TIME_DEBUG("read pl, nedges: " << polyline.nedges());

	//
	//
	// auto pls = PolyLineGenerator::extract_subpolylines(polyline);
	// //
	// DEBUG(pls.size());
	// PolyLine& pl = *pls[0];

	PolyLine& pl = polyline;

	/* 
	bool flip_normals = false;
	double default_sigma = 0.05; // wendland
	//double default_sigma = 0.05; // gaussian

	// 1 point
	//auto e = pl.iter_edges().begin().h ; sdf.add_func((e.from().pos().xy() + e.to().pos().xy())/2, 1., {0., 0.}, 1.);
	

	// SDFPointInit::circle(sdf, 10, [](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, r};});
	// SDFPointInit::shape(sdf, pl, [](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, r/2};});
	// SDFPointInit::shape_multiple(sdf, pl, 10, [](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, r};});
	// SDFPointInit::equaly_spaced_shape(sdf, pl, 10, [](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, r};});

	// auto algo = sdffitting::TestElliptic(sdf, samples);
	auto samples = samples::compute_edges_samples_normals(pl, 100, flip_normals);

	//auto algo = sdffitting::AlphaBetaOnlyAddFitter(sdf, samples);
	auto algo = sdffitting::ClusteringFitter(sdf, samples);
	
	algo.K = 2;
	algo.margin = 1.;

	algo.angular_threshold = std::numbers::pi / 18;
	algo.min_distance_points = 0.05;
	algo.cluster_min_size = 2;
	algo.error_threshold = 0.1;

	algo.ADD_POINT_ERR_THRESHOLD = 1.;
	algo.default_sigma_add = default_sigma;

	algo.fix_alpha_zero = false;
	algo.fix_beta_zero = false;

	algo.MIN_IMPROVEMENT = 1e-8;

	algo.lambda_distance = 1.;
	algo.lambda_gradient = 1.;

	int IT = 10;
	// algo.fit(IT, IT, output_dir);

	algo.resolve_ls();

	*/

	// TEST Elliptic
	
	
	SDF_Elliptic sdf(std::move(rbf));
	// UDF_Elliptic sdf(std::move(rbf));


	auto NORMAL = [](const vec2& v) -> vec2 { return {-v.y, v.x}; };
	SDFPointInit::shape(sdf, pl, [&](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, NORMAL(n) * r*0.5, r*0.5*n};});
	// SDFPointInit::shape_multiple(sdf, pl, 10, [&](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, NORMAL(n) * r, n * r/4 };});
	// SDFPointInit::equaly_spaced_shape(sdf, pl, 10, [&](auto p, auto n, auto r) -> FunctionElliptic { return {p, 0., n, NORMAL(n) * r, n * r/2 };});
	//
	//

	TIME_DEBUG("sdf points");


	auto samples = samples::compute_edges_samples_normals(pl, 5);

	TIME_DEBUG("samples");

	// auto samples = samples::compute_edges_samples_normals_include_end(pl, 2);
	// auto samples = samples::compute_equally_spaced_samples_normals(pl, 400);


	// auto algo = sdffitting_elliptic::TestEllipse(sdf, samples);
	auto algo = sdffitting_elliptic::TODONAMEFitter(sdf, samples);

	//
	algo.adapt_minor(1.0);
	TIME_DEBUG("adapt minor");
	//
	algo.decimate(1e-4, false);
	//
	TIME_DEBUG("decimate");

	// algo.adapt_minor(1.0);
	

	//
	// algo.adapt_minor(1.); 
	// algo.adjust_major(20.0); 
	// //
	// algo.adjust_minor(10.);
	//
	algo.fit_ellipses_radius(0.05, 2);
	TIME_DEBUG("fit ellipse");

	double max_area = 0;
	double mean_area = 0;
	for (const auto& f : sdf.fonctions) {
		
		double area = f.ellipse_major.norm() * f.ellipse_minor.norm();

		max_area = std::max(max_area, area);
		mean_area += area;
	}
	mean_area /= static_cast<double>(sdf.fonctions.size());

	TIME_DEBUG("mean : " << mean_area << " max : " << max_area);


	
	algo.lambda_distance = 100.;

	algo.resolve_ls();
	//
	
	// DEBUG("init eigen sovler");
	// TestEigen::test(sdf, samples);
	// TestEigen::testBFGS(sdf, samples);
	// TestEigen::BFGS_LS::testBFGS_LS(sdf, samples);
	// DEBUG("end eigen solver");
	
	DEBUG("optimization");
	TIME_DEBUG_INIT();

	sdf.optimize();

	sdf.max_neighbors = 20;

	TIME_DEBUG("sdf optimization");

	

	
	// OUTPUT
	output::export_samples(samples, output_dir + "samples.csv");
	TIME_DEBUG("export_samples");
	output::export_polyline(pl, output_dir + "polyline.csv");
	TIME_DEBUG("export_polyline");
	output::export_sdf(sdf, output_dir + "sdf_params.csv");
	TIME_DEBUG("export_sdf");
	// output::export_samples_error(samples, sdf, output_dir + "samples_error.csv", algo.lambda_distance, algo.lambda_gradient);
	// output::sample_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");
	//
	double samples_offset = 0.05; // = default_sigma
	output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "sdf.csv", 1000);
	TIME_DEBUG("sample_sdf");
	// std::ofstream(output_dir + "last.json") << sdf.to_json().dump(2);

	// std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;



	// EXPSMIN smin;
	// RFUNC_MIN smin;
	// SMOOTHMIN smin(
	// 		[](double x) { 
	// 		if (x <= -1.0) return 0.0;
	// 		if (x >= 1.0)  return x;
	// 		return (x * (2.0 + x)+1)/4.0 ;
	// 		},
	// 		[](double x) { 
	// 		if (x <= -1.0) return 0.0;
	// 		if (x >= 1.0)  return 1.0;
	// 		return (1.0+x)/2.0 ;
	// 		});
	//
	// SDF_Smoothmin sdf(pl, smin);
	//
	//

	// auto samples = samples::compute_edges_samples_normals(pl, 1);
	// output::export_samples_error(samples, sdf, output_dir + "samples_error.csv");

	// output::export_polyline(pl, output_dir + "polyline.csv");
	// double samples_offset = 0.2; // = default_sigma
	// output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "sdf.csv", 500);
	//

	return 0;




}

