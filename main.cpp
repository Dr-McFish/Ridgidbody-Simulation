
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <polyscope/polyscope.h>

#include "collision.h"
#include "rendering.h"
#include "physics_engine.h"


int main() {
	polyscope::init();

	Eigen::Vector3f a(1,2,3);


	struct physics_system system;// = initialise_system();

	struct collider csphere;
	csphere.type = COLIDER_SPHERE;
	csphere.u.sphere_colider = {.radius = 0.5f};

	struct rendering::mesh cube = rendering::create_cube("cube1");

	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(1, 1, 0), csphere, cube);

	struct collider cplane;
	cplane.type = COLIDER_HALF_SPACE;
	cplane.u.half_space_colider = { .normal = Eigen::Vector3f(0,1,0) };

	struct rendering::mesh square = rendering::create_square("plane1");
	square.scale = 5.f;

	add_imovablbe(system, Eigen::Vector3f(0, 0, 0), cplane, square);


	//if desired, set up a userCallback to add ImGui UI elements

	
	
	while(true) {

		//full_integration_step1(test_bodyl, 1./60.);
		//full_integration_step1(system.ridgidbodyies[0], 1./30.);
		integration_step(system);
		physys_render_update(system);
		printf("step\n");
	
		//auto contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		//visualise_collisions(system, contacts);
		//free(contacts);

		//std::cout << system.ridgidbodyies[0].x << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
