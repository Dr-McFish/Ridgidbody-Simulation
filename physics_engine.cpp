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
#include <eigen3/Eigen/src/Core/ArithmeticSequence.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <eigen3/Eigen/src/SparseCore/SparseMatrix.h>
#include <iostream>
#include "algebra_utils.h"
#include "collision.h"
#include "quaternion_helper.h"
#include "rendering.h"
#include "complementarity_solver.h"

// Gravity
static const Eigen::Vector3f g = 9.81 * Eigen::Vector3f(0, -1, 0);
static const float restitution_factor = 0.090;

Eigen::Quaternionf ith_Q(struct physics_system& system, int i) {
	return Eigen::Quaternionf(system.s(7*i +3), system.s(7*i +4), system.s(7*i +5), system.s(7*i +6) );
}
Eigen::Vector3f ith_x(struct physics_system& system, int i) {
	return system.s(Eigen::seq(7*i, 7*i +3));
}

Eigen::Vector3f ith_v(struct physics_system& system, int i) {
	return system.u(Eigen::seq(6*i, 6*i +2));
}
Eigen::Vector3f ith_omega(struct physics_system& system, int i) {
	return system.u(Eigen::seq(6*i+3, 6*i +5));
}

Eigen::Vector3f ith_force(struct physics_system& system, int i) {
	auto f =  system.force(Eigen::seq(6*i, 6*i +2));
	return f;
}



Eigen::Matrix3f ith_Iinv(struct physics_system& system, int i) {
	auto Q = ith_Q(system, i);
	Eigen::Matrix3f R = Q.matrix();

	return R * system.Ibody_inv[i] * R.transpose();
}



struct ridgidbody ith_body(struct physics_system& system, int i) {
	struct ridgidbody r;
	r.Ibody_inv = system.Ibody_inv[i];
	r.inv_mass = system.invmass[i];

	r.x = ith_x(system, i);
	r.Q = ith_Q(system, i);

	r.v = ith_v(system, i);
	r.omega = ith_omega(system, i);

	r.Iinv = ith_Iinv(system, i);

	return r;
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
		const Eigen::Vector3f n_k = contact_k.contact_normal;
		const Eigen::Vector3f t_k1 = contact_k.contact_tangent1;
		const Eigen::Vector3f t_k2 = contact_k.contact_tangent2;

		const Eigen::Vector3f r_ki = contact_k.contact_pos - ith_x(system, contact_k.bodyi_id);
		const Eigen::Vector3f r_ki_x_n  = cross_product(r_ki, n_k);
		const Eigen::Vector3f r_ki_x_t1 = cross_product(r_ki, t_k1);
		const Eigen::Vector3f r_ki_x_t2 = cross_product(r_ki, t_k2);

		const Eigen::Vector3f r_kj = contact_k.contact_pos - ith_x(system, contact_k.bodyj_id);		
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
				6*body_idx + repeat, 6*body_idx + repeat, system.invmass[body_idx]
			));
		}

		// Iinv block
		Eigen::Matrix3f Iinv = ith_Iinv(system, body_idx);
		for (int matx = 0; matx < 3; matx++) {
			for (int maty = 0; maty < 3; maty++) {
				tripletListMinv.push_back(Eigen::Triplet<float>(
				// i´, j´ J_i´,j´
				6*body_idx + 3 + matx, 6*body_idx + 3 + maty, Iinv(matx, maty)
			));
			}
		}
	}

	Minv.setFromTriplets(tripletListMinv.begin(), tripletListMinv.end());
	return Minv;
}

// Eigen::VectorXf u_generlised_velocity(physics_system& system) {
// 	const int N = system.ridgidbody_count;

// 	// Eigen::VectorXf u(6 * N);

// 	// for (int i = 0; i < N; i++) {
// 	// 	ridgidbody& body_i = system.ridgidbodyies[i];

// 	// 	u(6*i     + 0) = body_i.x(0);
// 	// 	u(6*i     + 1) = body_i.x(1);
// 	// 	u(6*i     + 2) = body_i.x(2);
// 	// 	u(6*i + 3 + 0) = body_i.omega(0);
// 	// 	u(6*i + 3 + 1) = body_i.omega(1);
// 	// 	u(6*i + 3 + 2) = body_i.omega(2);
// 	// }

// 	return system.u;
// }

// Eigen::VectorXf forces_ext_generlised(physics_system& system) {
// 	const int N = system.ridgidbody_count;

// 	Eigen::VectorXf forces(6 * N);

// 	for (int i = 0; i < N; i++) {
// 		ridgidbody& body_i = system.ridgidbodyies[i];

// 		forces(6*i     + 0) = body_i.force(0);
// 		forces(6*i     + 1) = body_i.force(1);
// 		forces(6*i     + 2) = body_i.force(2);
// 		forces(6*i + 3 + 0) = body_i.torque(0);
// 		forces(6*i + 3 + 1) = body_i.torque(1);
// 		forces(6*i + 3 + 2) = body_i.torque(2);
// 	}

// 	return forces;
// }

// // calculates and adds on the contact forces based on (wat) solved LPC
// void lambda_to_forces(physics_system& system, contact_list* contact_table, int K, Eigen::VectorXf lambda) {
// 	///std::cout << lambda << std::endl;

// 	for(int idx = 0; idx < K; idx++) {
// 		Eigen::Vector3f ri = contact_table[idx].contact_pos - ith_x(system, contact_table->bodyi_id);
// 		Eigen::Vector3f rj = contact_table[idx].contact_pos - ith_x(system, contact_table->bodyj_id);
		
// 		// body i on  body j
// 		Eigen::Vector3f force_ij = lambda(3*idx   ) * contact_table[idx].contact_normal
// 								 + lambda(3*idx +1) * contact_table[idx].contact_tangent1;
// 								 + lambda(3*idx +2) * contact_table[idx].contact_tangent2;
// 		// force_ji = - force_ij

// 		Eigen::Vector3f torque_i = cross_product(ri, -force_ij);
// 		Eigen::Vector3f torque_j = cross_product(rj,  force_ij);

// 		int i = contact_table->bodyi_id;
// 		int j = contact_table->bodyj_id;
// 		system.force(Eigen::seq(6*i  , 6*i +2)) += -force_ij;
// 		system.force(Eigen::seq(6*i+3, 6*i +5)) += torque_i;

// 		system.force(Eigen::seq(6*j  , 6*j +2)) += force_ij;
// 		system.force(Eigen::seq(6*j+3, 6*j +5)) += torque_j;
// 	}
// }


void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial, struct collider& collider, struct rendering::mesh& mesh) {
	system.ridgidbody_count += 1;
	system.invmass.push_back(1.f/mass_kg);
	system.Ibody_inv.push_back(Inertia_body);
	system.s.conservativeResize(7 * system.ridgidbody_count);
	system.s(Eigen::seq(system.ridgidbody_count - 8, system.ridgidbody_count-8 +3)) = x_initial;
	system.s(Eigen::seq(system.ridgidbody_count - 4, system.ridgidbody_count-1)) = Eigen::Vector4f(1,0,0,0);

	system.u.conservativeResize(6 * system.ridgidbody_count);

	system.force.conservativeResize(6 * system.ridgidbody_count);
	system.u(Eigen::seq(system.ridgidbody_count - 8, system.ridgidbody_count-1)) << 0,0,0 , 0,0,0;


	system.colliders.push_back(collider);
	system.mesh.push_back(mesh);
}

void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, collider& collider, struct rendering::mesh& mesh) {
	system.ridgidbody_count += 1;
	system.invmass.push_back(0.f); // 1/m = 0
	system.Ibody_inv.push_back(Eigen::Matrix3f::Zero()); // I⁻¹ = 0
	system.s.conservativeResize(7 * system.ridgidbody_count);
	system.s(Eigen::seq( 7 *system.ridgidbody_count - 8, 7 * system.ridgidbody_count-8 +3)) = x_initial;
	system.s(Eigen::seq( 7 *system.ridgidbody_count - 4, 7 * system.ridgidbody_count-1)) = Eigen::Vector4f(1,0,0,0);

	system.u.conservativeResize(6 * system.ridgidbody_count);

	system.force.conservativeResize(6 * system.ridgidbody_count);
	system.u(Eigen::seq(6* system.ridgidbody_count - 6, system.ridgidbody_count-1)) << 0,0,0 , 0,0,0;


	system.colliders.push_back(collider);
	system.mesh.push_back(mesh);
}

// blue - no collisions
// yellow - in contact
// red - penetration (bad theoreticly)
void visualise_collisions(struct physics_system& system, struct contact_list* contacts) {
	const glm::vec3 blue = glm::vec3(0.1, 0.1, 0.6);
	const glm::vec3 yellow = glm::vec3(0.6, 0.6, 0.0);
	const glm::vec3 red = glm::vec3(0.8, 0.1, 0.1);
	for (size_t i = 0; i < system.ridgidbody_count; i++) {
		system.mesh[i].mesh->setSurfaceColor(blue);
	}

	while (contacts != NULL) {
		const glm::vec3 coloring = contacts->penetration ? red : yellow;
		//if(contacts->penetration) {printf("red\n");} else {printf("yellow\n");} 
		system.mesh[contacts->bodyi_id].mesh->setSurfaceColor(coloring);
		system.mesh[contacts->bodyj_id].mesh->setSurfaceColor(coloring);
		
		contacts = contacts->next;
	}
}


// void initialize_rigidbody(struct ridgidbody& r, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial) {
// 	r.inv_mass = 1/mass_kg;
// 	r.Ibody_inv = Inertia_body.inverse();

// 	// State variables
// 	r.x = x_initial;
// 	r.Q = Eigen::Quaternionf::Identity();
// 	// Derived quantities (auxiliary variables)

// 	//  Computed quantities
// 	r.force =  Eigen::Vector3f::Zero();
// 	r.torque = Eigen::Vector3f::Zero();

// 	r.mesh = NULL;
// }

void physys_render_update(struct physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		// todo
		Eigen::Quaternionf Q(system.s(7*i +3), system.s(7*i +4), system.s(7*i +5), system.s(7*i +6) );
		rendering::update_mesh(system.mesh[i], system.s(Eigen::seq(7*i, 7*i +2)), Q);
	}
}

// struct physics_system initialise_system() {
// 	assert(false); // todo
// 	return system;
// }


void sove_for_contact_forces(struct physics_system& system) {
	
}
