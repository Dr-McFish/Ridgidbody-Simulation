#include "physics_engine.h"

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include "algebra_utils.h"
#include "collision.h"
#include "quaternion_helper.h"
#include "rendering.h"


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

void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial, struct collider& collider) {
	system.ridgidbody_count += 1;

	{
		// struct collider* tmp = (struct collider*)calloc(system.ridgidbody_count, sizeof(struct collider));
		// memcpy(tmp, system.colliders, (system.ridgidbody_count -1) * sizeof(struct collider*));
		// free(system.colliders);
		// system.colliders = tmp;
		//assert(NULL != system.colliders);
		system.colliders = (struct collider*)realloc(system.colliders, system.ridgidbody_count * sizeof(struct collider));
	}

	//system.ridgidbodyies = (struct ridgidbody*)realloc(system.ridgidbodyies, system.ridgidbody_count * sizeof(struct ridgidbody));
	{
		// struct ridgidbody* tmp = (struct ridgidbody*)calloc(system.ridgidbody_count, sizeof(struct ridgidbody));
		// memcpy(tmp, system.ridgidbodyies, (system.ridgidbody_count -1) * sizeof(struct ridgidbody*));
		// free(system.ridgidbodyies);
		// system.ridgidbodyies = tmp;
		//assert(NULL != system.ridgidbodyies);
		system.ridgidbodyies = (struct ridgidbody*)realloc(system.ridgidbodyies ,system.ridgidbody_count* sizeof(struct ridgidbody));
	}

	// TODO pas malin
	for (size_t i = 0 ; i < system.ridgidbody_count; i++) {
		system.colliders[i].pos = &system.ridgidbodyies[i].x;
		// printf("%ld\n", i);
		// assert(system.colliders[i].pos == &system.ridgidbodyies[i].x);
		// assert(*system.colliders[i].pos == system.ridgidbodyies[i].x);
	}

	struct collider& last_collider = system.colliders[system.ridgidbody_count -1];
	last_collider.type = collider.type;
	copy_collider(&collider, &last_collider);

	struct ridgidbody& last_ridgidbody = system.ridgidbodyies[system.ridgidbody_count -1];
	last_collider.pos = &last_ridgidbody.x;

	initialize_rigidbody(last_ridgidbody, mass_kg, Inertia_body, x_initial);
}

void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, collider& collider) {
	add_body(system, 1.0, Eigen::Matrix3f::Identity(), x_initial, collider);

	struct ridgidbody& last_ridgidbody = system.ridgidbodyies[system.ridgidbody_count -1];
	// immovable
	last_ridgidbody.Ibody_inv = Eigen::Matrix3f::Zero();
	last_ridgidbody.inv_mass = 0;
}

// blue - no collisions
// yellow - in contact
// red - penetration (bad theoreticly)
void visualise_collisions(struct physics_system& system, struct contact_list* contacts) {
	const glm::vec3 blue = glm::vec3(0.1, 0.1, 0.6);
	const glm::vec3 yellow = glm::vec3(0.6, 0.6, 0.0);
	const glm::vec3 red = glm::vec3(0.8, 0.1, 0.1);
	for (size_t i = 0; i < system.ridgidbody_count; i++) {
		if(NULL != system.ridgidbodyies[i].mesh) {
			system.ridgidbodyies[i].mesh->mesh->setSurfaceColor(blue);
		}
	}

	while (contacts != NULL) {
		const glm::vec3 coloring = contacts->penetration ? red : yellow;
		//if(contacts->penetration) {printf("red\n");} else {printf("yellow\n");} 
		if(NULL != system.ridgidbodyies[contacts->solid1_id].mesh){
			system.ridgidbodyies[contacts->solid1_id].mesh->mesh->setSurfaceColor(coloring);
			system.ridgidbodyies[contacts->solid2_id].mesh->mesh->setSurfaceColor(coloring);
		}
		
		contacts = contacts->next;
	}
}

void calculate_contact_velocity(struct physics_system& system, struct contact_list* contacts) {
	if(contacts == NULL) {
		return;
	}

	struct ridgidbody& body1 = system.ridgidbodyies[contacts->solid1_id];
	struct ridgidbody& body2 = system.ridgidbodyies[contacts->solid2_id];

	Eigen::Vector3f relative_pos1 = contacts->contact_pos - body1.x;
	Eigen::Vector3f relative_pos2 = contacts->contact_pos - body2.x;

	Eigen::Vector3f contact_velocity1 = body1.v + cross_product(body1.omega, relative_pos1);
	Eigen::Vector3f contact_velocity2 = body2.v + cross_product(body2.omega, relative_pos2);

	contacts->contact_velocity = contact_velocity2 - contact_velocity1;
	contacts->colliding_contact = dot_product(contacts->contact_velocity, contacts->contact_normal) < -system.minimum_nonpentration_velocity;

	calculate_contact_velocity(system, contacts->next);
}

void existsts_penetration_colliding_contact(struct contact_list* contacts, bool* E_colliding_contact, bool* E_penetration) {
	assert(NULL != E_colliding_contact);
	assert(NULL != E_penetration);\

	*E_colliding_contact = false;
	*E_penetration = false;
	
	while (NULL != contacts && !(*E_colliding_contact && *E_penetration)) {
		
		*E_colliding_contact =	*E_colliding_contact || contacts->colliding_contact;
		*E_penetration = 		*E_penetration 		 || contacts->penetration;

		contacts = contacts->next;
	}
}

void bisective_integration_step(struct physics_system& system, float sub_timestep, int maxdepth) {
	printf("depth: %d\n", maxdepth);
	if(system.stoped) {
		return;
	}
	else if(maxdepth <= 0) {
		fprintf(stderr, "Warning: maximum bisection depth reached! consder augmenting collision epsilon or decresing speed\n");
		system.stoped = true;
		return;
	}

	struct ridgidbody* old_state = (struct ridgidbody*)calloc(system.ridgidbody_count, sizeof(struct ridgidbody));
	memcpy(old_state, system.ridgidbodyies, system.ridgidbody_count * sizeof(struct ridgidbody));

	// Integration

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		compute_force_and_torque1(body, sub_timestep);
		compute_derived1(body);
	}


	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		integration_only_step1(body, sub_timestep);
	}

	// Collision/Penetration detction
	{
		struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		bool E_colliding_contact, E_penetration;
		existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
		
		if (NULL != contacts ){
			printf("penetration_depth : %2.3f\n", contacts->penetration_depth);
		} else {
			printf("no contact!\n");
			std::cout << system.ridgidbodyies[0].x << std::endl;
		}
		
		if(E_penetration)
		{
			// take one substep back
	
			free(system.ridgidbodyies); // restore old state
			system.ridgidbodyies = old_state;
			// AGH!
			for(size_t i = 0; i < system.ridgidbody_count; i++) {
				system.colliders[i].pos = &system.ridgidbodyies[i].x;
			}

			bisective_integration_step(system, 0.5 * sub_timestep, maxdepth -1);
		} 
		else if(E_colliding_contact)
		{
			// We are done
			system.stoped = true;
			free(old_state);
			visualise_collisions(system, contacts);
		}
		else
		{ 
			// take one substep forward
			bisective_integration_step(system, 0.5 * sub_timestep, maxdepth -1);
			free(old_state);
		}

		//visualise_collisions(system, contacts);
		free_contact_list(contacts);
	}

}

void integration_step0(struct physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		compute_force_and_torque1(body, system.base_timestep_seconds);
		compute_derived1(body);
	}
	// Collision detction
	struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
	visualise_collisions(system, contacts);
	free_contact_list(contacts);
	
	// TODO : Contact resolution

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		integration_only_step1(body, system.base_timestep_seconds);
	}
}

void integration_step(struct physics_system& system) {
	struct ridgidbody* old_state = (struct ridgidbody*)calloc(system.ridgidbody_count, sizeof(struct ridgidbody));
	memcpy(old_state, system.ridgidbodyies, system.ridgidbody_count * sizeof(struct ridgidbody));

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		compute_force_and_torque1(body, system.base_timestep_seconds);
		compute_derived1(body);
	}


	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		integration_only_step1(body, system.base_timestep_seconds);
	}

	
	// Collision detction
	struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
	visualise_collisions(system, contacts);

	bool E_colliding_contact, E_penetration;
	existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
	if (E_penetration)
	{
		
		//restore old state:
		free(system.ridgidbodyies);
		system.ridgidbodyies = old_state;
		/// AGH!
		for(size_t i = 0; i < system.ridgidbody_count; i++) {
			system.colliders[i].pos = &system.ridgidbodyies[i].x;
		}
		
		bisective_integration_step(system, 0.5 * system.base_timestep_seconds, 32);
		assert(system.stoped);

		free_contact_list(contacts);
		contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		visualise_collisions(system, contacts);
		free_contact_list(contacts);
	}
	else 
	{
		visualise_collisions(system, contacts);
		free_contact_list(contacts);
		free(old_state);
		old_state = NULL;
	}
	
	// TODO : Contact resolution

}


void initialize_rigidbody(struct ridgidbody& r, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial) {
	r.inv_mass = 1/mass_kg;
	r.Ibody_inv = Inertia_body.inverse();

	// State variables
	r.x = x_initial;
	r.Q = Eigen::Quaternionf::Identity();
	r.p = Eigen::Vector3f(1, 10, 5); // Quantité de mouvement
	r.L = Eigen::Vector3f(.1, 0.5, 0); // moment cinetique

	// Derived quantities (auxiliary variables)
	compute_derived1(r);

	//  Computed quantities
	r.force =  Eigen::Vector3f::Zero();
	r.torque = Eigen::Vector3f::Zero();

	r.mesh = NULL;
}

void physys_render_update(struct physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& ridgidbody = system.ridgidbodyies[i];
		rendering::update_mesh(*ridgidbody.mesh, ridgidbody.x, ridgidbody.Q);
	}
}

struct physics_system initialise_system() {
	struct physics_system system;
	system.ridgidbody_count = 0;
	system.ridgidbodyies = NULL;
	system.colliders = NULL;
	system.base_timestep_seconds = 1.f/30.f;
	system.remaining_seconds_until_next_timestep = 0.f;

	// 1/8 de distance nessesaire pour brisser le contact en 1 pas.
	system.minimum_nonpentration_velocity = (0.125f * penetration_epsilon_m) / system.base_timestep_seconds;

	//system.contact_count = 0;
	//system.contacts = NULL;
	return system;
}
