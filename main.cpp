
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


	struct physics_system system;

	struct collider csphere1;
	csphere1.type = COLIDER_SPHERE;
	csphere1.u.sphere_colider = {.radius = 0.5f};

	struct rendering::mesh cube1 = rendering::create_sphere("cube");
	add_body(system, 10, Eigen::Matrix3f::Identity(), Eigen::Vector3f(+3, 0.50, 0), csphere1, cube1);\

	struct collider cplane;
	cplane.type = COLIDER_HALF_SPACE;
	cplane.u.half_space_colider = { .normal = Eigen::Vector3f(0,1,0) };


	struct rendering::mesh square = rendering::create_square("plane1");
	square.scale = 5.f;


	add_imovablbe(system, Eigen::Vector3f(0, 0, 0), cplane, square);

	system.base_timestep_seconds = 1.f/30.f;

	system.springs.push_back((struct spring) {
		.body_i = 0,
		.body_j = 1,
		.r_i = Eigen::Vector3f(0.f, -0.f, 0.f),
		.r_j = Eigen::Vector3f(0.f, 1.5f, 0.f),
		.k = 10.f,
		.l0 = 0.5,
	});
	start_spring_visualisation(system);

	printf("E_c,E_p,E_m\n");
	
	int i =0;
	for (int step = 0; step < 20000; step++) {
		i++;
		integration_step(system);
		physys_render_update(system);
	
		float E_p = potential_energy(system);
		float E_c = kinetic_energy(system);
		printf("%3.3f,%3.3f,%3.3f\n", E_p, E_c, E_p +E_c);
		polyscope::frameTick(); // renders one UI frame, returns immediately
	}
}
