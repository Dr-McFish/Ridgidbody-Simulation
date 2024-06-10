
#include <cassert>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <polyscope/polyscope.h>

#include "collision.h"
#include "rendering.h"
#include "physics_engine.h"

/* Stack Scenario */

int main() {
	polyscope::init();

	Eigen::Vector3f a(1,2,3);


	struct physics_system system;// = initialise_system();

	struct collider csphere1;
	csphere1.type = COLIDER_SPHERE;
	csphere1.u.sphere_colider = {.radius = 0.5f};

	struct collider csphere2 = csphere1;


	char noms[60];
	for(int i = 0; i < 12; i++) {
		noms[2*i] = 'a' + i;
		noms[2*i +1] = '\0';
		struct rendering::mesh cube1 = rendering::create_sphere(&noms[2*i]);
		add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(-2 + 1.2*(i/ 5), 1.2*(i%5) + 1, 0), csphere1, cube1);\
	}

	struct collider cplane;
	cplane.type = COLIDER_HALF_SPACE;
	cplane.u.half_space_colider = { .normal = Eigen::Vector3f(0,1,0) };


	struct collider cplaneN;
	cplaneN.type = COLIDER_HALF_SPACE;
	cplaneN.u.half_space_colider = { .normal = Eigen::Vector3f(1,0,0) };


	struct collider cplaneS;
	cplaneS.type = COLIDER_HALF_SPACE;
	cplaneS.u.half_space_colider = { .normal = Eigen::Vector3f(-1,0,0) };


	struct collider cplaneE;
	cplaneE.type = COLIDER_HALF_SPACE;
	cplaneE.u.half_space_colider = { .normal = Eigen::Vector3f(0, 0,1) };

	struct collider cplaneW;
	cplaneW.type = COLIDER_HALF_SPACE;
	cplaneW.u.half_space_colider = { .normal = Eigen::Vector3f(0, 0,-1) };

	struct rendering::mesh square = rendering::create_square("plane1");
	square.scale = 5.f;


	// add_imovablbe(system, Eigen::Vector3f(-4, 0, 0), cplaneN, square);
	// add_imovablbe(system, Eigen::Vector3f(+4, 0, 0), cplaneS, square);
	// add_imovablbe(system, Eigen::Vector3f(0, 0, -2), cplaneE, square);
	// add_imovablbe(system, Eigen::Vector3f(0, 0, +2), cplaneW, square);
	add_imovablbe(system, Eigen::Vector3f(0, 0, 0), cplane, square);

	//if desired, set up a userCallback to add ImGui UI elements

	system.base_timestep_seconds = 1.f/30.f;

	int i = 0;
	while(true) {
		//full_integration_step1(test_bodyl, 1./60.);
		//full_integration_step1(system.ridgidbodyies[0], 1./30.);
		if( i % 300 < 15) {
		integration_step(system);
		physys_render_update(system);
		printf("step\n");
		}
	
		//auto contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		//visualise_collisions(system, contacts);
		//free(contacts);

		//std::cout << ith_v(system, 0) << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
		//getchar();
		i++;
	}
}
