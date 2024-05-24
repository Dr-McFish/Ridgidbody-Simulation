#ifndef PHYSICS_ENGINE_H
# define PHYSICS_ENGINE_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include "collision.h"
#include "rendering.h"

struct ridgidbody {
	// Constant quantities
	float inv_mass;
	Eigen::Matrix3f Ibody_inv;

	// State variables
	Eigen::Vector3f x;
	Eigen::Quaternionf Q;
	Eigen::Vector3f p; // Quantité de mouvement
	Eigen::Vector3f L; // moment cinetique

	// Derived quantities (auxiliary variables)
	Eigen::Vector3f v;
	Eigen::Vector3f omega;
	Eigen::Matrix3f Iinv;

	//  Computed quantities
	Eigen::Vector3f force;
	Eigen::Vector3f torque;

	//optional
	struct rendering::mesh* mesh;

};

void initialize_rigidbody(struct ridgidbody& r, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial);

struct physics_system {
	size_t ridgidbody_count;
	struct ridgidbody* ridgidbodyies;
	struct collider* colliders; // `rigidbodies[k]` is assosiated with collider `colliders[k]`
	float timestep_seconds;

	// TODO: champs de forces
};

void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial, struct collider& c);
void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, collider& c);
void full_integration_step1(struct ridgidbody& body, float delta_t);
void integration_step(struct physics_system& system);
void physys_render_update(struct physics_system& system);


void physys_render_update(struct physics_system& system);
void visualise_collisions(struct physics_system& system, struct contact_list* contacts);
// struct contact_list* detect_collisions(struct physics_system);
// void integration_step(struct physics_system& system);

// void physics_step(struct physics_system& system);



#endif
