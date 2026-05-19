#include "ultimaille/polyline.h"
#include <cfloat>
#include <csignal>
#include <utility>
#include "rbf.hpp"
#include "test_init.hpp"
#include "sdf.hpp"
#include "algo.hpp"

using namespace UM;

int main(int argc, char** argv) {


	const std::string output_dir = OUTPUT_DIR + std::string("test/");
	std::filesystem::remove_all(output_dir);
	std::filesystem::create_directories(output_dir);
	const std::string input_dir = INPUT_DIR;


	PolyLine pl;
	auto rbf = std::make_unique<InverseMultiquadric>();
	SDF sdf(std::move(rbf));



	//int n = 3;
	//PolyLineGenerator::regular_polygon(pl, n);
	
	//SDFPointInit::shape(sdf, pl, 1., 1.5);
	// SDFPointInit::circle(sdf, n, 1., std::numbers::pi / n);

	//PolyLineGenerator::random_polygon(pl, n);

	/*
	int n = 10;
	PolyLineGenerator::random_polygon(pl, n);
	SDFPointInit::shape(sdf, pl);
	*/


	//SDFPointInit::circle(sdf, 3, 1., std::numbers::pi / 3);

	
	PolyLineGenerator::read_from_file(pl, input_dir + "duck.geogram");
	// SDFPointInit::shape(sdf, pl);
	



	auto algo = sdffitting::Fitter(pl, sdf, sdffitting::compute_samples_normals(pl, 1));

	algo.step_sigma = 0.01;

	algo.step_point = 0.1;


	/*
	algo.lb_sigma = 0.1;
	algo.ub_sigma = 3.;

	algo.lb_point = -1.5;
	algo.ub_point = 1.5;
	*/

	algo.lb_sigma = -DBL_MAX;
	algo.ub_sigma = DBL_MAX;
	algo.lb_point = -DBL_MAX;
	algo.ub_point = DBL_MAX;


	algo.ADD_POINT_ERR_THRESHOLD = 1.;

	algo.default_sigma_add = 3.;


	//algo.fit(100, 100, output_dir);
	//algo.fit_moving_points(500, 100, output_dir);
	//algo.fit_add_point(5000, 100, output_dir);
	//algo.fit_moving_add_points(5000, 500, output_dir);

	//algo.fit_adaptative_step(500, 500, output_dir);
	//algo.fit_adaptative_step_moving_points(500, 100, output_dir);
	algo.fit_adaptative_step_moving_add_points(4000, 200, output_dir);
	return 0;
}
