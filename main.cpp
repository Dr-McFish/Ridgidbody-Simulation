
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <polyscope/polyscope.h>

#include "rendering.h"
#include "physics_engine.h"


int main() {
	polyscope::init();

	Eigen::Vector3f a(1,2,3);


	struct ridgidbody test_bodyl;
	initialize_rigidbody(test_bodyl, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(1, 1, 0));
	struct rendering::mesh test_cube = rendering::create_square("test_cube");	

	// if desired, set up a userCallback to add ImGui UI elements

	while(true) {

		full_integration_step1(test_bodyl, 1./60.);
		rendering::update_mesh(test_cube, test_bodyl.x, test_bodyl.Q);
		//std::cout << test_bodyl.x << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
