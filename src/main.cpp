#include "ultimaille/polyline.h"
#include <cfloat>
#include <csignal>
#include <cstddef>
#include <utility>
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"
#include "algo.hpp"
#include "types.hpp"
#include "output.hpp"
#include "samples.hpp"

using namespace UM;

int main(int argc, char** argv) {


	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;


	PolyLine pl;
	auto rbf = std::make_unique<WendlandC2>();
	

	SDF sdf(std::move(rbf));

	int n = 8; PolyLineGenerator::regular_polygon(pl, n);
	//PolyLineGenerator::random_polygon(pl, n);

	//PolyLineGenerator::read_from_file(pl, input_dir + "duck.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "cerf.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "e.obj");
	//PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");
	// SDFPointInit::shape(sdf, pl);
	

	bool flip_normals = false;
	double default_sigma = 1.; // wendland
	//double default_sigma = 0.05; // gaussian

	//SDFPointInit::shape(sdf, pl, 0., default_sigma, flip_normals);
	
	SDFPointInit::shape_sigma_edge_length(sdf, pl, 0., 0.6, flip_normals);

	//SDFPointInit::equaly_spaced_shape(sdf, pl, 10, 0., default_sigma, flip_normals);
	

	auto samples = samples::compute_edges_samples_normals(pl, 50, flip_normals);

	//auto samples = sdffitting::samples::compute_equally_spaced_samples_normals(pl, 150);

	auto algo = sdffitting::AlphaBetaOnlyAddFitter(sdf, samples);


	/*
	algo.lb_sigma = 0.1;
	algo.ub_sigma = 3.;

	algo.lb_point = -1.5;
	algo.ub_point = 1.5;
	*/

	algo.lb_sigma = 0;
	algo.ub_sigma = DBL_MAX;
	algo.lb_point = -DBL_MAX;
	algo.ub_point = DBL_MAX;


	algo.ADD_POINT_ERR_THRESHOLD = 1.;

	algo.default_sigma_add = default_sigma;

	algo.fix_alpha_zero = false;
	algo.fix_beta_zero = false;

	algo.MIN_IMPROVEMENT = 1e-8;


	int IT = 1;
	algo.fit(IT, IT, output_dir);



	// OUTPUT
	output::export_samples(samples, output_dir + "samples.csv");
	output::export_polyline(pl, output_dir + "polyline.csv");
	output::export_sdf(sdf, output_dir + "sdf_params.csv");
	output::export_samples_error(samples, sdf, output_dir + "samples_error.csv");
	//output::sample_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");
	//
	double samples_offset = 0.2; // = default_sigma
	output::sample_sdf(sdf, -1. - samples_offset, -1. - samples_offset, 1. + samples_offset, 1. + samples_offset, output_dir + "sdf.csv", 100);
	std::ofstream(output_dir + "last.json") << sdf.to_json().dump(2);

	std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;
	return 0;
}
