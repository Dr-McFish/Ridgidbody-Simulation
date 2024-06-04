#include "physics_engine.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <eigen3/Eigen/src/SparseCore/SparseMatrix.h>
#include <iostream>
#include "algebra_utils.h"
#include "collision.h"
#include "glm/fwd.hpp"
#include "quaternion_helper.h"
#include "rendering.h"
#include "complementarity_solver.h"

// Gravity
static const Eigen::Vector3f g = 9.81 * Eigen::Vector3f(0, -1, 0);
static const float restitution_factor = 0.09;

// In pace update the derived quantities
void compute_derived1(struct ridgidbody& body) {
	Eigen::Matrix3f R = body.Q.matrix();
	body.Iinv = R * body.Ibody_inv * R.transpose();
	body.v = body.inv_mass * body.p;
	body.omega = body.Iinv * body.L;
}

void compute_derived_all(physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		compute_derived1(system.ridgidbodyies[i]);
	}
}

Eigen::SparseMatrix<float> Jmatrix(physics_system& system, contact_list* contact_table, int num_contacts){
	const int K = num_contacts;
	const int N = system.ridgidbody_count;

	Eigen::SparseMatrix<float> J(3*K, 6*N);
	std::vector<Eigen::Triplet<float>> tripletListJ;
	tripletListJ.reserve(K*2*3*3);

	for (int k = 0; k < K; k++) {
		const contact_list& contact_k = contact_table[k];


		// if it goes the wrong direction, n_k might be wrong way
		const Eigen::Vector3f n_k  = contact_k.contact_normal;
		const Eigen::Vector3f t_k1 = contact_k.contact_tangent1;
		const Eigen::Vector3f t_k2 = contact_k.contact_tangent2;

		const Eigen::Vector3f r_ki = contact_k.contact_pos - system.ridgidbodyies[contact_k.bodyi_id].x;
		const Eigen::Vector3f r_ki_x_n  = cross_product(r_ki, n_k);
		const Eigen::Vector3f r_ki_x_t1 = cross_product(r_ki, t_k1);
		const Eigen::Vector3f r_ki_x_t2 = cross_product(r_ki, t_k2);

		const Eigen::Vector3f r_kj = contact_k.contact_pos - system.ridgidbodyies[contact_k.bodyj_id].x;		
		const Eigen::Vector3f r_kj_x_n  = cross_product(r_kj, n_k);
		const Eigen::Vector3f r_kj_x_t1 = cross_product(r_kj, t_k1);
		const Eigen::Vector3f r_kj_x_t2 = cross_product(r_kj, t_k2);
		
		Eigen::Matrix3f Ji_klin;
		Ji_klin << -n_k(0) , -n_k(1) , -n_k(2) ,
				   -t_k1(0), -t_k1(1), -t_k1(2),
				   -t_k2(0), -t_k2(1), -t_k2(2);

		Eigen::Matrix3f Jj_klin;
		Jj_klin << n_k(0) , n_k(1) , n_k(2) ,
				   t_k1(0), t_k1(1), t_k1(2),
				   t_k2(0), t_k2(1), t_k2(2);
		
		Eigen::Matrix3f Ji_kang;
		Ji_kang << -r_ki_x_n (0), -r_ki_x_n (1), -r_ki_x_n (2),
				   -r_ki_x_t1(0), -r_ki_x_t1(1), -r_ki_x_t1(2),
				   -r_ki_x_t2(0), -r_ki_x_t2(1), -r_ki_x_t2(2);

		Eigen::Matrix3f Jj_kang;
		Jj_kang << r_kj_x_n (0), r_kj_x_n (1), r_kj_x_n (2),
				   r_kj_x_t1(0), r_kj_x_t1(1), r_kj_x_t1(2),
				   r_kj_x_t2(0), r_kj_x_t2(1), r_kj_x_t2(2);
		
		for (int matx = 0; matx < 3; matx++) {
			for (int maty = 0; maty < 3; maty++) {
				tripletListJ.push_back(Eigen::Triplet<float>(
					// i´, j´ J_i´,j´
					3*k + matx, 6*contact_k.bodyi_id + maty, Ji_klin(matx,maty)
				));

				tripletListJ.push_back(Eigen::Triplet<float>(
					// i´, j´ J_i´,j´
					3*k + matx, 6*contact_k.bodyj_id + maty, Jj_klin(matx,maty)
				));

				tripletListJ.push_back(Eigen::Triplet<float>(
					// i´, j´ J_i´,j´
					3*k + matx, 6*contact_k.bodyi_id + 3 + maty, Ji_kang(matx,maty)
				));

				tripletListJ.push_back(Eigen::Triplet<float>(
					// i´, j´ J_i´,j´
					3*k + matx, 6*contact_k.bodyj_id + 3 + maty, Jj_kang(matx,maty)
				));
			}
		}
	}
	J.setFromTriplets(tripletListJ.begin(), tripletListJ.end());

	return J;
}

Eigen::SparseMatrix<float> generalizedMass_matrix_inv(physics_system& system) {
	const int N = system.ridgidbody_count;

	std::vector<Eigen::Triplet<float>> tripletListMinv;
	tripletListMinv.reserve(N*12);

	Eigen::SparseMatrix<float> Minv(6*N, 6*N);

	for(int body_idx = 0; body_idx < N; body_idx++) {
		
		// diagonal block
		for (int repeat = 0; repeat < 3; repeat++) {
			tripletListMinv.push_back(Eigen::Triplet<float>(
				// i´, j´ J_i´,j´
				6*body_idx + repeat, 6*body_idx + repeat, system.ridgidbodyies[body_idx].inv_mass
			));
		}

		// Iinv block
		for (int matx = 0; matx < 3; matx++) {
			for (int maty = 0; maty < 3; maty++) {
				tripletListMinv.push_back(Eigen::Triplet<float>(
				// i´, j´ J_i´,j´
				6*body_idx + 3 + matx, 6*body_idx + 3 + maty, system.ridgidbodyies[body_idx].Ibody_inv(matx, maty)
			));
			}
		}
	}

	Minv.setFromTriplets(tripletListMinv.begin(), tripletListMinv.end());
	return Minv;
}

Eigen::VectorXf u_generlised_velocity(physics_system& system) {
	const int N = system.ridgidbody_count;

	Eigen::VectorXf u(6 * N);

	for (int i = 0; i < N; i++) {
		ridgidbody& body_i = system.ridgidbodyies[i];

		u(6*i     + 0) = body_i.x(0);
		u(6*i     + 1) = body_i.x(1);
		u(6*i     + 2) = body_i.x(2);
		u(6*i + 3 + 0) = body_i.omega(0);
		u(6*i + 3 + 1) = body_i.omega(1);
		u(6*i + 3 + 2) = body_i.omega(2);
	}

	return u;
}

Eigen::VectorXf forces_ext_generlised(physics_system& system) {
	const int N = system.ridgidbody_count;

	Eigen::VectorXf forces(6 * N);

	for (int i = 0; i < N; i++) {
		ridgidbody& body_i = system.ridgidbodyies[i];

		forces(6*i     + 0) = body_i.force(0);
		forces(6*i     + 1) = body_i.force(1);
		forces(6*i     + 2) = body_i.force(2);
		forces(6*i + 3 + 0) = body_i.torque(0);
		forces(6*i + 3 + 1) = body_i.torque(1);
		forces(6*i + 3 + 2) = body_i.torque(2);
	}

	return forces;
}

// calculates and adds on the contact forces based on (wat) solved LPC
void lambda_to_forces(physics_system& system, contact_list* contact_table, int K, Eigen::VectorXf lambda) {
	///std::cout << lambda << std::endl;

	for(int idx = 0; idx < K; idx++) {
		ridgidbody& body_i = system.ridgidbodyies[contact_table[idx].bodyi_id];
		ridgidbody& body_j = system.ridgidbodyies[contact_table[idx].bodyj_id];

		Eigen::Vector3f ri = contact_table[idx].contact_pos - body_i.x;
		Eigen::Vector3f rj = contact_table[idx].contact_pos - body_j.x;
		
		// body i on  body j
		Eigen::Vector3f force_ij = lambda(3*idx   ) * contact_table[idx].contact_normal
								 + lambda(3*idx +1) * contact_table[idx].contact_tangent1;
								 + lambda(3*idx +2) * contact_table[idx].contact_tangent2;
		// force_ji = - force_ij

		Eigen::Vector3f torque_i = cross_product(ri, -force_ij);
		Eigen::Vector3f torque_j = cross_product(rj,  force_ij);

		body_i.force  -= -force_ij;
		body_i.torque -= torque_i;

		body_j.force  -= force_ij;
		body_j.torque -= torque_j;
	}
}

// Call once external forces are known and stored in the rigidbodiesr
// each ellement of the contact table is one long
void compute_contact_forces(physics_system& system, struct contact_list* contact_list, float timestep) {
	
	int K; //N
	struct contact_list* contact_table = list_to_array(contact_list, &K);
	
	Eigen::SparseMatrix<float> Minv = generalizedMass_matrix_inv(system);
	Eigen::SparseMatrix<float> J = Jmatrix(system, contact_table, K);

	// printf("Jmatrix ar = \n");
	// std::cout << J.toDense() << std::endl;

	struct linear_complementarity_problem lpc;

	lpc.A = timestep * J * Minv * J.transpose();

	auto u = u_generlised_velocity(system);
	auto f_ext = forces_ext_generlised(system);

	lpc.b = J * (u + timestep * (Minv * f_ext) );


	const float mu_friction_coef = 0;
	lpc.lamba_min = Eigen::VectorXf(3* K);
	lpc.lamba_max = Eigen::VectorXf(3* K);

	for (int j = 0; j < K; j++) {
		lpc.lamba_min(j + 0) = 0;
		lpc.lamba_min(j + 1) = -mu_friction_coef * 0; // friction at some point maybe
		lpc.lamba_min(j + 2) = -mu_friction_coef * 0; // friction at some point maybe

		lpc.lamba_max(j + 0) = +INFINITY;
		lpc.lamba_max(j + 1) = mu_friction_coef * 0; // friction at some point maybe
		lpc.lamba_max(j + 2) = mu_friction_coef * 0; // friction at some point maybe
	}

	Eigen::VectorXf lambda = pgs_solve(&lpc, 10);

	// printf("Amatrix = \n");
	// std::cout << lpc.A.toDense() << std::endl;
	// printf("b_vec = \n");
	// std::cout << lpc.b << std::endl;

	// printf("lambda = \n");
	// std::cout << lambda << std::endl;

	printf("force_before = \n");
	std::cout << system.ridgidbodyies[0].force << std::endl;
	
	lambda_to_forces(system, contact_table, K, lambda);
	
	printf("force_after = \n");
	std::cout << system.ridgidbodyies[0].force << std::endl;

	//assert(false); // crash

	free(contact_table);
}

// Only accumulates gravity for now
void compute_force_and_torque1(struct ridgidbody& body) {
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
	compute_force_and_torque1(body);
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
	const glm::vec3 orange  = glm::vec3(0.9f, 0.65f, 0.f);
	for (size_t i = 0; i < system.ridgidbody_count; i++) {
		if(NULL != system.ridgidbodyies[i].mesh) {
			system.ridgidbodyies[i].mesh->mesh->setSurfaceColor(blue);
		}
	}

	while (contacts != NULL) {
		glm::vec3 coloring;
		if(contacts->penetration) {
			coloring = red; 
		} else if (contacts->colliding_contact) {
			coloring = orange;
		} else {
			coloring = yellow;
		}

		//if(contacts->penetration) {printf("red\n");} else {printf("yellow\n");} 
		if(NULL != system.ridgidbodyies[contacts->bodyi_id].mesh){
			system.ridgidbodyies[contacts->bodyi_id].mesh->mesh->setSurfaceColor(coloring);
			system.ridgidbodyies[contacts->bodyj_id].mesh->mesh->setSurfaceColor(coloring);
		}
		
		contacts = contacts->next;
	}
}

void calculate_contact_velocity(struct physics_system& system, struct contact_list* contacts) {
	if(contacts == NULL) {
		return;
	}

	struct ridgidbody& body_i = system.ridgidbodyies[contacts->bodyi_id];
	struct ridgidbody& body_j = system.ridgidbodyies[contacts->bodyj_id];

	Eigen::Vector3f relative_pos_i = contacts->contact_pos - body_i.x;
	Eigen::Vector3f relative_pos_j = contacts->contact_pos - body_j.x;

	Eigen::Vector3f contact_velocity_i = body_i.v + cross_product(body_i.omega, relative_pos_i);
	Eigen::Vector3f contact_velocity_j = body_j.v + cross_product(body_j.omega, relative_pos_j);

	contacts->contact_velocity = contact_velocity_j - contact_velocity_i;
	contacts->contact_velocity_projected = dot_product(contacts->contact_velocity, contacts->contact_normal);
	contacts->colliding_contact = contacts->contact_velocity_projected > system.minimum_nonpentration_velocity;

	calculate_contact_velocity(system, contacts->next);
}

void existsts_penetration_colliding_contact(struct contact_list* contacts, bool* E_colliding_contact, bool* E_penetration) {
	assert(NULL != E_colliding_contact);
	assert(NULL != E_penetration);

	*E_colliding_contact = false;
	*E_penetration = false;
	
	while (NULL != contacts && !(*E_colliding_contact && *E_penetration)) {
		
		*E_colliding_contact =	*E_colliding_contact || contacts->colliding_contact;
		*E_penetration = 		*E_penetration 		 || contacts->penetration;

		//printf("penetration %f", contacts->penetration_depth);

		contacts = contacts->next;
	}
}

void bisective_integration_step(struct physics_system& system, float sub_timestep, int maxdepth) {
	//printf("depth: %d\n", maxdepth);
	if(system.stoped) {
		return;
	}
	else if(maxdepth <= 0) {
		fprintf(stderr, "Warning: maximum bisection depth reached! consder augmenting collision epsilon or decresing speed\n");
		//system.stoped = true;
		return;
	}

	struct ridgidbody* old_state = (struct ridgidbody*)calloc(system.ridgidbody_count, sizeof(struct ridgidbody));
	memcpy(old_state, system.ridgidbodyies, system.ridgidbody_count * sizeof(struct ridgidbody));

	// Integration

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		compute_derived1(body);
		compute_force_and_torque1(body);
	}

	// normal forces
	{
	struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
	calculate_contact_velocity(system, contacts);

	if(contacts != NULL) {
		compute_contact_forces(system, contacts, sub_timestep);
	}
	}

	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		struct ridgidbody& body = system.ridgidbodyies[i];
		integration_only_step1(body, sub_timestep);
	}

	// Collision/Penetration detction
	{
		struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		calculate_contact_velocity(system, contacts);
		bool E_colliding_contact, E_penetration;
		existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
		
		if(E_colliding_contact)
		{
			// We are done
			system.remaining_seconds_until_next_timestep -= sub_timestep;
			//system.stoped = true;
			free(old_state);
			//visualise_collisions(system, contacts);
		}
		else if(E_penetration)
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
		else
		{ 
			// take one substep forward
			system.remaining_seconds_until_next_timestep -= sub_timestep;
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
		compute_force_and_torque1(body);
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

void resolve_colliding_contacts(struct physics_system& system, struct contact_list* contacts) {
	if(NULL == contacts) {
		return;
	} else if(contacts->contact_velocity_projected >= -system.minimum_nonpentration_velocity) {
		// not colliding

		resolve_colliding_contacts(system, contacts->next);
	}

	struct ridgidbody& body_a = system.ridgidbodyies[contacts->bodyi_id];
	struct ridgidbody& body_b = system.ridgidbodyies[contacts->bodyj_id];

	Eigen::Vector3f relative_pos_a = contacts->contact_pos - body_a.x;
	Eigen::Vector3f relative_pos_b = contacts->contact_pos - body_b.x;

	//float restitution = (contacts->contact_velocity_projected <= -system.minimum_nonpentration_velocity*4) ?  0.f : restitution_factor;

	float j = -((1.f + restitution_factor) * contacts->contact_velocity_projected)/ 
			(body_b.inv_mass + body_a.inv_mass
			 + dot_product(contacts->contact_normal, cross_product(body_a.Iinv * (cross_product(relative_pos_a, contacts->contact_normal)),
																   relative_pos_a))
			 + dot_product(contacts->contact_normal, cross_product(body_b.Iinv * (cross_product(relative_pos_b, contacts->contact_normal)),
																   relative_pos_b))
			);
	
	Eigen::Vector3f impulse = j * contacts->contact_normal;

	body_a.p -= impulse;
	body_b.p += impulse;

	body_a.L -= cross_product(relative_pos_a, impulse);
	body_b.L += cross_product(relative_pos_a, impulse);

	compute_derived1(body_a);
	compute_derived1(body_b);

	resolve_colliding_contacts(system, contacts->next);
}

void integration_step(struct physics_system& system) {
	if(system.stoped) {
		return;
	}
	system.remaining_seconds_until_next_timestep = system.base_timestep_seconds;

	while(system.remaining_seconds_until_next_timestep > 1E-6)
	{
		struct ridgidbody* old_state = (struct ridgidbody*)calloc(system.ridgidbody_count, sizeof(struct ridgidbody));
		memcpy(old_state, system.ridgidbodyies, system.ridgidbody_count * sizeof(struct ridgidbody));

		for(size_t i = 0; i < system.ridgidbody_count; i++) {
			struct ridgidbody& body = system.ridgidbodyies[i];
			compute_derived1(body);
			compute_force_and_torque1(body);
		}

		{
		struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		calculate_contact_velocity(system, contacts);

		if(contacts != NULL) {
			compute_contact_forces(system, contacts, system.remaining_seconds_until_next_timestep);
		}
		}
		
		for(size_t i = 0; i < system.ridgidbody_count; i++) {
			struct ridgidbody& body = system.ridgidbodyies[i];
			integration_only_step1(body, system.remaining_seconds_until_next_timestep);
		}
		compute_derived_all(system);

		// Collision detction
		struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
		calculate_contact_velocity(system, contacts);
		

		bool E_colliding_contact, E_penetration;
		existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
		if (E_penetration)
		{
			//printf(" boop\n");
			
			//restore old state:
			free(system.ridgidbodyies);
			system.ridgidbodyies = old_state;
			/// AGH!
			for(size_t i = 0; i < system.ridgidbody_count; i++) {
				system.colliders[i].pos = &system.ridgidbodyies[i].x;
			}
			
			bisective_integration_step(system, 0.5 * system.base_timestep_seconds, 10);
			//assert(system.stoped);

			//debug
			struct contact_list* contacts = collision_detectoion(system.ridgidbody_count, system.colliders);
			calculate_contact_velocity(system, contacts);
			visualise_collisions(system, contacts);

			bool E_colliding_contact, E_penetration;
			existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
			assert(!E_penetration);
			assert(E_colliding_contact);
			
			// Bouce ressolution
			resolve_colliding_contacts(system, contacts);
			calculate_contact_velocity(system, contacts);
			
			existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
			assert(!E_penetration);
			assert(!E_colliding_contact);

			free_contact_list(contacts);
			//debug

		}
		else 
		{
			system.remaining_seconds_until_next_timestep = 0.f;
			visualise_collisions(system, contacts);
			if(E_colliding_contact) {
				resolve_colliding_contacts(system, contacts);
				calculate_contact_velocity(system, contacts);

				bool E_colliding_contact, E_penetration;
				existsts_penetration_colliding_contact(contacts, &E_colliding_contact, &E_penetration);
				assert(!E_penetration);
				assert(!E_colliding_contact);

			}


			free(old_state);
			old_state = NULL;
		}
		
		free_contact_list(contacts);
	}	

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

	system.stoped = false;
	// 1/8 de distance nessesaire pour brisser le contact en 1 pas.
	system.minimum_nonpentration_velocity = (0.5f * penetration_epsilon_m) / system.base_timestep_seconds;

	//system.contact_count = 0;
	//system.contacts = NULL;
	return system;
}
