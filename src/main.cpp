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

using namespace UM;

int main(int argc, char** argv) {


	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;


	PolyLine pl;
	auto rbf = std::make_unique<WendlandC2>();
	

	SDF sdf(std::move(rbf));

	int n = 3;
	//PolyLineGenerator::regular_polygon(pl, n);
	//PolyLineGenerator::random_polygon(pl, n);

	PolyLineGenerator::read_from_file(pl, input_dir + "duck.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "cerf.geogram");
	//PolyLineGenerator::read_from_file(pl, input_dir + "u.obj");
	// SDFPointInit::shape(sdf, pl);
	

	bool flip_normals = false;
	double default_sigma = 0.2; // wendland
	//double default_sigma = 0.05; // gaussian

	//SDFPointInit::shape(sdf, pl, 0., default_sigma, flip_normals);
	

	//SDFPointInit::equaly_spaced_shape(sdf, pl, 50, 0., default_sigma, flip_normals);
	

	auto samples = sdffitting::samples::compute_edges_samples_normals(pl, 5, flip_normals);

	//auto samples = sdffitting::samples::compute_equally_spaced_samples_normals(pl, 150);

	auto algo = sdffitting::Fitter(sdf, samples);
	//auto algo = sdffitting::Fitter(sdf, samples);;

	algo.step_sigma = 0.0;

	algo.step_point = 0.1;


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

	algo.MIN_IMPROVMENT = 1e-8;



	//algo.fit(100, 100, output_dir);
	//algo.fit_moving_points(500, 100, output_dir);
	//algo.fit_add_point(5000, 100, output_dir);
	//algo.fit_moving_add_points(5000, 500, output_dir);

	//algo.fit_adaptative_step(500, 500, output_dir);
	//algo.fit_adaptative_step_add_points(500, 50, output_dir);
	//algo.fit_adaptative_step_moving_points(500, 100, output_dir);
	//algo.fit_adaptative_step_moving_add_points(1000, 200, output_dir);
	
	//algo.fit_alpha_beta_add(0, 0, output_dir);

	algo.resolve_ls();


	// OUTPUT
	output::export_samples(samples, output_dir + "samples.csv");
	output::export_polyline(pl, output_dir + "polyline.csv");
	output::export_sdf(sdf, output_dir + "sdf_params.csv");
	//output::sample_sdf(sdf, -3, -3, 3, 3, output_dir + "sdf.csv");
	//
	//
	output::sample_sdf(sdf, -1. - default_sigma, -1. - default_sigma, 1. + default_sigma, 1. + default_sigma, output_dir + "sdf.csv", 200);
	std::ofstream(output_dir + "last.json") << sdf.to_json().dump(2);

	std::cout << "FINAL SDF:\n" << sdf.to_string() << std::endl;
	return 0;
}
