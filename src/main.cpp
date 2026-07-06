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
	
	
	// int n = 16; PolyLineGenerator::regular_polygon(polyline, n);
	// PolyLineGenerator::demi_circle(polyline, 10, 20);
	// PolyLineGenerator::read_from_file(polyline, input_dir + "duck.geogram");

	// PolyLineGenerator::read_from_file(pl, input_dir + "lapinpluslisseplusdense_quad_mesh.obj");
	
	// PolyLineGenerator::read_from_file(polyline, input_dir + "cerf.geogram");
	// PolyLineGenerator::read_from_file(polyline, input_dir + "o.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o1.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o2.obj");
	// PolyLineGenerator::read_from_file(polyline, input_dir + "ooo.obj");
	// PolyLineGenerator::read_from_file(polyline, input_dir + "ooo2.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "e.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");

	// PolyLineGenerator::isoceles_triangle(pl, std::numbers::pi / 4);
	

	PolyLineGenerator::read_polyline_directly(polyline, input_dir + "lapinpluslisseplusdense_bord.obj");
	// PolyLineGenerator::read_polyline_directly(polyline, input_dir + "rocket.geogram");

	// PolyLineGenerator::read_polyline_directly(polyline, input_dir + "patte.geogram");
	
	// PolyLineGenerator::read_from_file(polyline, input_dir + "square_bord.obj");


	TIME_DEBUG("read pl, nedges: " << polyline.nedges());


	auto pls = PolyLineGenerator::extract_subpolylines_corner(polyline);
	// //
	TIME_DEBUG("found subpolylines :" << pls.size());
	int k = 0;
	std::filesystem::create_directories(output_dir + "subpolylines/");
	for (auto& subpl : pls) {
		write_by_extension(output_dir + "subpolylines/" + std::to_string(k) + ".geogram", *subpl);

		DEBUG("subpl " << k << " npoints: " << subpl->nverts() << " nedges: " << subpl->nedges());
		++k;



	}
	// PolyLine& pl = *pls[0];

	PolyLine& pl = polyline;

	// double scale = PolyLineGenerator::normalize_pl(pl);
	double scale = 1.;



	SDF_Elliptic sdf(std::move(rbf));
	// UDF_Elliptic sdf(std::move(rbf));

	int polyline_segmant_nsample = 6;
	double target_function_width = 0.1 * scale;
	int target_interpolant_neighbors = 6;
	double lambda_distance = 10.;

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
	sdf.neighborhood_size = target_interpolant_neighbors + 4;






	double max_area = 0;
	double mean_area = 0;
	for (const auto& f : sdf.fonctions) {
		
		double area = f.ellipse_major.norm() * f.ellipse_minor.norm();

		max_area = std::max(max_area, area);
		mean_area += area;
	}
	mean_area /= static_cast<double>(sdf.fonctions.size());

	TIME_DEBUG("mean : " << mean_area << " max : " << max_area);

	

	
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
	double samples_offset = 0.2; // = default_sigma
	output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "sdf.csv", 500);
	TIME_DEBUG("sample_sdf");
	// std::ofstream(output_dir + "last.json") << sdf.to_json().dump(2);

	// std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;
	//

	return 0;




}

