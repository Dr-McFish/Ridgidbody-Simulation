
#include <cassert>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <polyscope/polyscope.h>

#include "collision.h"
#include "rendering.h"
#include "physics_engine.h"

/*spring scenario*/

int main() {
	polyscope::init();

	Eigen::Vector3f a(1,2,3);


	struct physics_system system;// = initialise_system();

	struct collider csphere1;
	csphere1.type = COLIDER_SPHERE;
	csphere1.u.sphere_colider = {.radius = 0.5f};

	struct rendering::mesh cube1 = rendering::create_cube("cube");
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(+3, 0.52, 0), csphere1, cube1);\

	struct collider cplane;
	cplane.type = COLIDER_HALF_SPACE;
	cplane.u.half_space_colider = { .normal = Eigen::Vector3f(0,1,0) };


	struct rendering::mesh square = rendering::create_square("plane1");
	square.scale = 5.f;


	// add_imovablbe(system, Eigen::Vector3f(-4, 0, 0), cplaneN, square);
	// add_imovablbe(system, Eigen::Vector3f(+4, 0, 0), cplaneS, square);
	// add_imovablbe(system, Eigen::Vector3f(0, 0, -2), cplaneE, square);
	// add_imovablbe(system, Eigen::Vector3f(0, 0, +2), cplaneW, square);
	add_imovablbe(system, Eigen::Vector3f(0, 0, 0), cplane, square);

	//if desired, set up a userCallback to add ImGui UI elements

	system.base_timestep_seconds = 1.f/30.f;

	system.springs.push_back((struct spring) {
		.body_i = 0,
		.body_j = 1,
		.r_i = Eigen::Vector3f(0.f, -0.f, 0.f),
		.r_j = Eigen::Vector3f(0.f, 1.5f, 0.f),
		.k = 2.f,
		.l0 = 0.5,
	});
	start_spring_visualisation(system);

	while(true) {

		//full_integration_step1(test_bodyl, 1./60.);
		//full_integration_step1(system.ridgidbodyies[0], 1./30.);
		integration_step(system);
		physys_render_update(system);
		//printf("step\n");
	
		//auto contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		//visualise_collisions(system, contacts);
		//free(contacts);

		printf("%3.3f, %3.3f\n", potential_energy(system), kinetic_energy(system));

		//std::cout << ith_v(system, 0) << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
