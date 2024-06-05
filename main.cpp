
#include <cassert>
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

	struct collider csphere1;
	csphere1.type = COLIDER_SPHERE;
	csphere1.u.sphere_colider = {.radius = 0.5f};

	struct collider csphere2 = csphere1;


	struct rendering::mesh cube1 = rendering::create_cube("cube1");
	struct rendering::mesh cube2 = rendering::create_cube("cube2");
	struct rendering::mesh cube3 = rendering::create_cube("cube3");
	struct rendering::mesh cube4 = rendering::create_cube("cube4");
	struct rendering::mesh cube5 = rendering::create_cube("cube5");
	struct rendering::mesh cube6 = rendering::create_cube("cube6");
	struct rendering::mesh cube7 = rendering::create_cube("cube7");

	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(2, 1, 0), csphere1, cube1);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(-2, 1, 0), csphere2, cube2);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(0, 1, 0), csphere1, cube3);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(0, 3, 0), csphere1, cube4);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(0, 5, 0), csphere1, cube5);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(0, 7, 0), csphere1, cube6);
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(0, 9, 0), csphere1, cube7);

	

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


	add_imovablbe(system, Eigen::Vector3f(-4, 0, 0), cplaneN, square);
	add_imovablbe(system, Eigen::Vector3f(+4, 0, 0), cplaneS, square);
	add_imovablbe(system, Eigen::Vector3f(0, 0, -4), cplaneE, square);
	add_imovablbe(system, Eigen::Vector3f(0, 0, +4), cplaneW, square);
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

		//std::cout << ith_v(system, 0) << std::endl;
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
