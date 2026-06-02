#include "algo_base.hpp"
#include "sdf_elliptic.hpp"
#include "ultimaille/polyline.h"
#include <cfloat>
#include <csignal>
#include <cstddef>
#include <filesystem>
#include <utility>
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"
#include "algo_imp1.hpp"
#include "algo_imp2.hpp"
#include "output.hpp"
#include "samples.hpp"
#include "algo_base_elliptic.hpp"

using namespace UM;

int main(int argc, char** argv) {

	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;


	PolyLine pl;
	auto rbf = std::make_unique<WendlandC2>();
	

	// SDF sdf(std::move(rbf));

	// SDF_Elliptic sdf(std::move(rbf));

	// int n = 4; PolyLineGenerator::regular_polygon(pl, n);
	//PolyLineGenerator::random_polygon(pl, n);

	// PolyLineGenerator::read_from_file(pl, input_dir + "duck.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "cerf.geogram");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "e.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");
	// SDFPointInit::shape(sdf, pl);
	

	/* 
	
	bool flip_normals = false;
	double default_sigma = 0.05; // wendland
	//double default_sigma = 0.05; // gaussian

	//SDFPointInit::shape(sdf, pl, 0., default_sigma, flip_normals);
	
	SDFPointInit::shape_sigma_edge_length(sdf, pl, 0., 0.5, flip_normals);

	// SDFPointInit::shape_multiple(sdf, pl, 1, 0., 1.1, flip_normals);

	//SDFPointInit::equaly_spaced_shape(sdf, pl, 12, 0., default_sigma, flip_normals);
	
// 1 point
	//auto e = pl.iter_edges().begin().h ; sdf.add_func((e.from().pos().xy() + e.to().pos().xy())/2, 1., {0., 0.}, 1.);
	

	auto samples = samples::compute_edges_samples_normals(pl, 100, flip_normals);

	//auto samples = sdffitting::samples::compute_equally_spaced_samples_normals(pl, 150);

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
	// SDF_Elliptic sdf(std::move(rbf));
	SDF sdf(std::move(rbf));

	int n = 3; PolyLineGenerator::regular_polygon(pl, n);
	// PolyLineGenerator::read_from_file(pl, input_dir + "duck.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "cerf.geogram");
	// PolyLineGenerator::read_from_file(pl, input_dir + "o.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "e.obj");
	// PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");


	/*
	for (const auto& e : pl.iter_edges()) {
		auto v = e.to().pos() - e.from().pos();

		sdf.add_func(
				{
				(e.from().pos() + (v/2)).xy(), 
				0., 
				{0., 0.}, 
				v.xy()/1, 
				UM::vec2(-v.y, v.x)/1
				}
				);
	}
	*/

	SDFPointInit::shape_sigma_edge_length(sdf, pl, 0., 1., false);

	auto samples = samples::compute_edges_samples_normals(pl, 100);
	// auto algo = sdffitting_elliptic::TestEllipse(sdf, samples);
	auto algo = sdffitting::TestCircular(sdf, samples);

	algo.resolve_ls();


	// OUTPUT
	output::export_samples(samples, output_dir + "samples.csv");
	output::export_polyline(pl, output_dir + "polyline.csv");
	output::export_sdf(sdf, output_dir + "sdf_params.csv");
	output::export_samples_error(samples, sdf, output_dir + "samples_error.csv", algo.lambda_distance, algo.lambda_gradient);
	//output::sample_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");
	//
	double samples_offset = 0.2; // = default_sigma
	output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "sdf.csv", 500);
	std::ofstream(output_dir + "last.json") << sdf.to_json().dump(2);

	std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;
	return 0;




}
