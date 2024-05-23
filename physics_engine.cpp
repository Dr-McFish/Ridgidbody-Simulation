#include "physics_engine.h"

#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include "collision.h"
#include "quaternion_helper.h"


// Gravity
static const Eigen::Vector3f g = 9.81 * Eigen::Vector3f(0, -1, 0);

// In pace update the derived quantities
void compute_derived1(struct ridgidbody& body) {
	Eigen::Matrix3f R = body.Q.matrix();
	body.Iinv = R * body.Ibody_inv * R.transpose();
	body.v = body.inv_mass * body.p;
	body.omega = body.Iinv * body.L;
}

// Only accumulates gravity for now
void compute_force_and_torque1(struct ridgidbody& body, float delta_t) {
	body.force = g;
	body.torque = Eigen::Vector3f::Zero();
}

void integration_only_step1(struct ridgidbody& body, float delta_t) {
	// Semi-implicit / "Symplectic" Euler
	body.p += delta_t * body.force;
	body.L += delta_t * body.torque;
	body.x += delta_t * body.v;
	body.Q = body.Q + q_sclale(delta_t, (quaternionf_from_vector_part(body.omega) * body.Q));
	body.Q.normalize();
}

void full_integration_step1(struct ridgidbody& body, float delta_t) {
	compute_force_and_torque1(body, delta_t);
	compute_derived1(body);
	integration_only_step1(body, delta_t);
}

void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial, struct collider collider) {
	system.ridgidbody_count += 1;

	system.colliders = (struct collider*)realloc(system.colliders, system.ridgidbody_count * sizeof(struct collider));
	system.ridgidbodyies = (struct ridgidbody*)realloc(system.ridgidbodyies, system.ridgidbody_count * sizeof(struct collider));

	system.colliders[system.ridgidbody_count -1].type = collider.type;
	copy_collider(collider, system.colliders[system.ridgidbody_count-1]);

	initialize_rigidbody(system.ridgidbodyies[system.ridgidbody_count -1], 
							mass_kg, Inertia_body, x_initial);
}

void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, collider collider) {
	system.ridgidbody_count += 1;
	
	system.colliders = (struct collider*)realloc(system.colliders, system.ridgidbody_count * sizeof(struct collider));
	system.ridgidbodyies = (struct ridgidbody*)realloc(system.ridgidbodyies, system.ridgidbody_count * sizeof(struct collider));

	system.colliders[system.ridgidbody_count -1].type = collider.type;
	copy_collider(collider, system.colliders[system.ridgidbody_count-1]);
	
	struct ridgidbody& last_ridgidbody = system.ridgidbodyies[system.ridgidbody_count -1];
	initialize_rigidbody(last_ridgidbody, 1.0, Eigen::Matrix3f::Identity(), x_initial);
	// immovable
	last_ridgidbody.Ibody_inv = Eigen::Matrix3f::Zero();
	last_ridgidbody.inv_mass = 0;
}

void integration_step(struct physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody body = system.ridgidbodyies[i];
		compute_force_and_torque1(body, system.timestep_seconds);
		compute_derived1(body);
	}
	// TODO : Collision detction
	// TODO : Contact resolution

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody body = system.ridgidbodyies[i];
		integration_only_step1(body, system.timestep_seconds);
	}
}

void initialize_rigidbody(struct ridgidbody& r, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial) {
	r.inv_mass = 1/mass_kg;
	r.Ibody_inv = Inertia_body.inverse();

	// State variables
	r.x = x_initial;
	r.Q = Eigen::Quaternionf::Identity();
	r.p = Eigen::Vector3f(1, 10, 5); // Quantité de mouvement
	r.L = Eigen::Vector3f(1, 5, 0); // moment cinetique

	// Derived quantities (auxiliary variables)
	compute_derived1(r);

	//  Computed quantities
	r.force =  Eigen::Vector3f::Zero();
	r.torque = Eigen::Vector3f::Zero();
}