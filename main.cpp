
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <polyscope/polyscope.h>

#include "collision.h"
#include "rendering.h"
#include "physics_engine.h"


// int main() {
// 	polyscope::init();

// 	Eigen::Vector3f a(1,2,3);


// 	struct ridgidbody test_bodyl;
// 	initialize_rigidbody(test_bodyl, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(1, 1, 0));
// 	struct rendering::mesh test_cube = rendering::create_square("test_cube");	

// 	// if desired, set up a userCallback to add ImGui UI elements

// 	while(true) {

// 		full_integration_step1(test_bodyl, 1./60.);
// 		rendering::update_mesh(test_cube, test_bodyl.x, test_bodyl.Q);
// 		//std::cout << test_bodyl.x << std::endl;
// 		polyscope::frameTick(); // renders one UI frame, returns immediately
// 	}
// }

int main() {
	polyscope::init();

	Eigen::Vector3f a(1,2,3);


	struct physics_system system;
	system.ridgidbody_count = 0;
	system.ridgidbodyies = NULL;
	system.colliders = NULL;
	system.timestep_seconds = 1./60.;

	struct collider csphere;
	csphere.pos = NULL;
	csphere.type = COLIDER_SPHERE;
	csphere.u.sphere_colider = {.radius = 0.5f};

	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(1, 1, 0), csphere);
	struct rendering::mesh cube = rendering::create_cube("cube1");
	system.ridgidbodyies[0].mesh = &cube;

	struct collider cplane;
	cplane.pos = NULL;
	cplane.type = COLIDER_HALF_SPACE;
	cplane.u.half_space_colider = { .normal = Eigen::Vector3f(0,1,0) };
	add_imovablbe(system, Eigen::Vector3f(0, 0, 0), cplane);

	struct rendering::mesh square = rendering::create_square("plane1");
	square.scale = 10.f;
	system.ridgidbodyies[1].mesh = &square;
	
	//if desired, set up a userCallback to add ImGui UI elements

	while(true) {

		//full_integration_step1(test_bodyl, 1./60.);
		full_integration_step1(system.ridgidbodyies[0], 1./60.);
		physys_render_update(system);
		auto contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		if(NULL != contacts) {
			printf("4\n");
		}
		visualise_collisions(system, contacts);
		free(contacts);
		
		//rendering::update_mesh(test_cube, test_bodyl.x, test_bodyl.Q);
		// std::cout << system.ridgidbodyies[0].x << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
