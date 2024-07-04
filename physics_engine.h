#ifndef PHYSICS_ENGINE_H
# define PHYSICS_ENGINE_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <vector>

#include "collision.h"
#include "polyscope/curve_network.h"
#include "rendering.h"

struct spring {
	int body_i;
	int body_j;

	// points de application dans le refferencielle des solides
	Eigen::Vector3f r_i;
	Eigen::Vector3f r_j;

	float k;
	float l0;
};

struct physics_system {
	size_t ridgidbody_count {0};

	// constatnts
	std::vector<Eigen::Matrix3f> Ibody_inv;
	std::vector<float> invmass;

	// 7N generalized position vector
	Eigen::VectorXf s;

	// 6N Generalised Velocity
	Eigen::VectorXf u;

	// 6N Generalized force(incl. torque)
	Eigen::VectorXf force;
	
	float base_timestep_seconds {1./30.};
	float remaining_seconds_until_next_timestep {0.f};
	float minimum_nonpentration_velocity {penetration_epsilon_m * 30};

	std::vector<struct collider> colliders;
	std::vector<struct rendering::mesh> mesh;
	std::vector<struct spring> springs;
	polyscope::CurveNetwork* spring_visualisation {NULL};
};
Eigen::Vector3f s_ith_x(Eigen::VectorXf& s, int i);
Eigen::Quaternionf s_ith_Q(Eigen::VectorXf& s, int i);

Eigen::Vector3f ith_x(struct physics_system& system, int i);
Eigen::Quaternionf ith_Q(struct physics_system& system, int i);
Eigen::Vector3f ith_v(struct physics_system& system, int i);
Eigen::Vector3f ith_omega(struct physics_system& system, int i);
Eigen::Vector3f ith_force(struct physics_system& system, int i);

void start_spring_visualisation(struct physics_system& system);
void physys_render_update(struct physics_system& system);
void physys_render_update(struct physics_system& system);
void visualise_collisions(struct physics_system& system, struct contact_list* contacts);

void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, 
			  Eigen::Vector3f x_initial, struct collider& c, struct rendering::mesh& m);
void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, 
					collider& c, struct rendering::mesh& mesh);

void integration_step(struct physics_system& system);

//Donnés :

float potential_energy(struct physics_system& system);
float kinetic_energy(struct physics_system& system);

#endif
