#include "physics_engine.h"


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
#include <vector>
#include "algebra_utils.h"
#include "collision.h"
#include "glm/fwd.hpp"
#include "polyscope/curve_network.h"
#include "rendering.h"
#include "complementarity_solver.h"

// Gravity
static const Eigen::Vector3f g = 9.81 * Eigen::Vector3f(0, -1, 0);
static const float restitution_factor = 0.6f;
const float mu_friction_coef = 2.f;


Eigen::Quaternionf s_ith_Q(Eigen::VectorXf& s, int i) {
	return Eigen::Quaternionf(s(7*i +3), s(7*i +4), s(7*i +5), s(7*i +6) );
}

Eigen::Quaternionf ith_Q(struct physics_system& system, int i) {
	return s_ith_Q(system.s, i);
}
Eigen::Vector3f s_ith_x(Eigen::VectorXf& s, int i) {
	return s(Eigen::seq(7*i, 7*i +2));
}

Eigen::Vector3f ith_x(struct physics_system& system, int i) {
	return system.s(Eigen::seq(7*i, 7*i +2));
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



void compute_spring_forces(physics_system& system) {
	for (size_t j = 0; j < system.springs.size();  j ++) {
		struct spring& spring = system.springs[j];
		int body_i = spring.body_i ;
		int body_j = spring.body_j ;

		Eigen::Vector3f r_i_world = ith_Q(system, body_i)._transformVector( spring.r_i );
		Eigen::Vector3f r_j_world = ith_Q(system, body_j)._transformVector( spring.r_j );

		Eigen::Vector3f spring_end_i = ith_x(system, body_i) + r_i_world;
		Eigen::Vector3f spring_end_j = ith_x(system, body_j) + r_j_world;

		float l = (spring_end_i - spring_end_j).norm();
		Eigen::Vector3f u = (1.f/l) * (spring_end_i - spring_end_j);
		float force_magnitude = spring.k * (l - spring.l0);

		auto f_i = system.force(Eigen::seq(6*body_i, 6*body_i +2));
		auto f_j = system.force(Eigen::seq(6*body_j, 6*body_j +2));

		auto torque_i = system.force(Eigen::seq(6*body_i + 3, 6*body_i +5));
		auto torque_j = system.force(Eigen::seq(6*body_j + 3, 6*body_j +5));

		Eigen::VectorXf spirng_force_i = - force_magnitude * u;

		f_i += spirng_force_i;
		f_j += - spirng_force_i;

		torque_i += cross_product(r_i_world, spirng_force_i);
		torque_j += cross_product(r_j_world, -spirng_force_i);
	}
}


void compute_forces(physics_system& system) {
	for (size_t i = 0; i < system.ridgidbody_count; i++) {
		auto force_i  = system.force(Eigen::seq(6*i    , 6*i + 2));
		auto torque_i = system.force(Eigen::seq(6*i + 3, 6*i + 5));

		force_i = g;


		Eigen::Vector3f omega = ith_omega(system, i);
		torque_i = - cross_product(omega, ith_Iinv(system, i) * omega);
	}

	compute_spring_forces(system);
}

Eigen::SparseMatrix<float> Jmatrix(physics_system& system, std::vector<struct contact_list> contact_table){
	const int K = contact_table.size();
	const int N = system.ridgidbody_count;

	Eigen::SparseMatrix<float> J(3*K, 6*N);
	std::vector<Eigen::Triplet<float>> tripletListJ;
	tripletListJ.reserve(K*2*3*3);

	for (int k = 0; k < K; k++) {
		const contact_list& contact_k = contact_table[k];

		const Eigen::Vector3f n_k = -contact_k.contact_normal;
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

//
Eigen::SparseMatrix<float> S_position_derivative_matrix(physics_system& system) {
	const int N = system.ridgidbody_count;

	std::vector<Eigen::Triplet<float>> tripletListS;
	tripletListS.reserve(N*12);

	Eigen::SparseMatrix<float> S(7*N, 6*N);

	for(int body_idx = 0; body_idx < N; body_idx++) {

		// diagonal block
		for (int repeat = 0; repeat < 3; repeat++) {
			tripletListS.push_back(Eigen::Triplet<float>(
				// i´, j´ J_i´,j´
				7*body_idx + repeat, 6*body_idx + repeat, 1
			));
		}

		Eigen::Quaternionf Q = ith_Q(system, body_idx);

		Eigen::Matrix<float, 4, 3> Q_d ;
		Q_d << -Q.x(), -Q.y(), -Q.z(),
			    Q.w(),  Q.z(), -Q.y(),
			   -Q.z(),  Q.w(),  Q.x(),
			    Q.y(), -Q.x(),  Q.w();
		
		Q_d *= 0.5;

		for (int matx = 0; matx < 4; matx++) {
			for (int maty = 0; maty < 3; maty++) {
				tripletListS.push_back(Eigen::Triplet<float>(
				// i´, j´ J_i´,j´
				7*body_idx + 3 + matx, 6*body_idx + 3 + maty, Q_d(matx, maty)
			));
			}
		}
	}

	S.setFromTriplets(tripletListS.begin(), tripletListS.end());
	return S;
}


void add_body(struct physics_system& system, float mass_kg, Eigen::Matrix3f Inertia_body, Eigen::Vector3f x_initial, struct collider& collider, struct rendering::mesh& mesh) {
	system.ridgidbody_count += 1;
	system.invmass.push_back(1.f/mass_kg);
	system.Ibody_inv.push_back(Inertia_body);
	system.s.conservativeResize(7 * system.ridgidbody_count);


	system.s(Eigen::seq(7*system.ridgidbody_count - 7    , 7*system.ridgidbody_count -7 +2)) = x_initial;
	system.s(Eigen::seq(7*system.ridgidbody_count - 7 + 3, 7*system.ridgidbody_count -1)) = Eigen::Vector4f(1,0,0,0);

	system.u.conservativeResize(6 * system.ridgidbody_count);

	system.force.conservativeResize(6 * system.ridgidbody_count);
	system.u(Eigen::seq(6*system.ridgidbody_count - 6, 6*system.ridgidbody_count-1)) << 0,0,0 , 0,0,0;


	system.colliders.push_back(collider);
	system.mesh.push_back(mesh);
}

void add_imovablbe(struct physics_system& system, Eigen::Vector3f x_initial, collider& collider, struct rendering::mesh& mesh) {
	system.ridgidbody_count += 1;
	system.invmass.push_back(0.f); // 1/m = 0
	system.Ibody_inv.push_back(Eigen::Matrix3f::Zero()); // I⁻¹ = 0
	system.s.conservativeResize(7 * system.ridgidbody_count);
	system.s(Eigen::seq( 7 *system.ridgidbody_count - 7    , 7 * system.ridgidbody_count-7 +2)) = x_initial;
	system.s(Eigen::seq( 7 *system.ridgidbody_count - 7 + 3, 7 * system.ridgidbody_count-1)) = Eigen::Vector4f(1,0,0,0);

	system.u.conservativeResize(6 * system.ridgidbody_count);
	system.force.conservativeResize(6 * system.ridgidbody_count);

	system.u(Eigen::seq(6* system.ridgidbody_count - 6, 6*system.ridgidbody_count-1)) << 0,0,0 , 0,0,0;


	system.colliders.push_back(collider);
	system.mesh.push_back(mesh);
}


// returns the lambda vector
Eigen::VectorXf compute_contact_impulses(int K, physics_system& system, Eigen::SparseMatrix<float> J, Eigen::SparseMatrix<float> Minv, float timestep)
{
	struct linear_complementarity_problem lcp;
	lcp.A = J * Minv * J.transpose();

	// bounce
	Eigen::VectorXf b_bouce = restitution_factor * J * system.u;
	for (int i = 0; i < K; i++) {
		b_bouce(3*i + 1) = 0.f;
		b_bouce(3*i + 2) = 0.f;
	}

	lcp.b = J * (system.u + timestep * (Minv * system.force) ) + b_bouce;

	lcp.lambda_min = Eigen::VectorXf(3* K);
	lcp.lambda_max = Eigen::VectorXf(3* K);

	for (int j = 0; j < K; j++) {
		lcp.lambda_min(3*j + 0) = 0;
		lcp.lambda_min(3*j + 1) = -mu_friction_coef * 0; // friction at some point maybe
		lcp.lambda_min(3*j + 2) = -mu_friction_coef * 0; // friction at some point maybe

		lcp.lambda_max(3*j + 0) = +INFINITY;
		lcp.lambda_max(3*j + 1) = mu_friction_coef * 0; // friction at some point maybe
		lcp.lambda_max(3*j + 2) = mu_friction_coef * 0; // friction at some point maybe
	}

	Eigen::VectorXf lambda = pgs_solve(&lcp, 10, mu_friction_coef);

	return lambda;
}

void normalize_quaternions(physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++)
	{
		auto Q = system.s(Eigen::seq(7*i +3, 7*i +6) );

		Q.normalize();
	}
}

Eigen::VectorXf force_impulses(struct physics_system& system, float timestep_seconds) {
	compute_forces(system);

	contact_list* contacts = collision_detectoion(system.colliders, system.s);
	visualise_collisions(system, contacts);
	std::vector<struct contact_list> contact_table = list_to_array(contacts);

	Eigen::SparseMatrix<float> Minv = generalizedMass_matrix_inv(system);
	Eigen::SparseMatrix<float> J = Jmatrix(system, contact_table);
	const int K = contact_table.size();


	Eigen::VectorXf lambda = compute_contact_impulses(K, system, J, Minv, timestep_seconds);

	return  (Minv * J.transpose() * lambda)
		  + (timestep_seconds * Minv * system.force);
}

void integration_step_symplectic(struct physics_system& system, float timestep_seconds)
{
	// Velocity update
	system.u += force_impulses(system, timestep_seconds);

	// Position update
	Eigen::SparseMatrix<float> S = S_position_derivative_matrix(system);
	system.s = system.s + timestep_seconds * S * system.u;

	normalize_quaternions(system);
}

void integration_step_explicit(struct physics_system& system, float timestep_seconds)
{
	Eigen::VectorXf impulses = force_impulses(system, timestep_seconds);

	// Position update
	Eigen::SparseMatrix<float> S = S_position_derivative_matrix(system);
	system.s = system.s + timestep_seconds * S * system.u;
	
	normalize_quaternions(system);

	// Velocity update
	system.u += impulses;

}

void integration_step_midpoint(struct physics_system& system, float timestep_seconds)
{
	Eigen::VectorXf u_old = Eigen::VectorXf(system.u);
	Eigen::VectorXf s_old = Eigen::VectorXf(system.s);

	integration_step_explicit(system, timestep_seconds/2);

	compute_forces(system);

	contact_list* contacts2 = collision_detectoion(system.colliders, system.s);
	visualise_collisions(system, contacts2);
	std::vector<struct contact_list> contact_table = list_to_array(contacts2);

	Eigen::SparseMatrix<float> Minv2 = generalizedMass_matrix_inv(system);
	Eigen::SparseMatrix<float> J2 = Jmatrix(system, contact_table);
	const int K2 = contact_table.size();

	linear_complementarity_problem lcp;

	// bounce
	Eigen::VectorXf b_bouce = restitution_factor * J2 * system.u;
	for (int i = 0; i < K2; i++) {
		b_bouce(3*i + 1) = 0.f;
		b_bouce(3*i + 2) = 0.f;
	}
	lcp.A = J2 * Minv2 * J2.transpose();
	lcp.b = J2 * (u_old + timestep_seconds * (Minv2 * system.force) ) + b_bouce;


	lcp.lambda_min = Eigen::VectorXf(3* K2);
	lcp.lambda_max = Eigen::VectorXf(3* K2);

	for (int j = 0; j < K2; j++) {
		lcp.lambda_min(3*j + 0) = 0;
		lcp.lambda_min(3*j + 1) = -mu_friction_coef * 0; // friction at some point maybe
		lcp.lambda_min(3*j + 2) = -mu_friction_coef * 0; // friction at some point maybe

		lcp.lambda_max(3*j + 0) = +INFINITY;
		lcp.lambda_max(3*j + 1) = mu_friction_coef * 0; // friction at some point maybe
		lcp.lambda_max(3*j + 2) = mu_friction_coef * 0; // friction at some point maybe
	}

	Eigen::VectorXf lambda = pgs_solve(&lcp, 10, mu_friction_coef);

	Eigen::VectorXf w = lcp.A * lambda + lcp.b;
	

	Eigen::SparseMatrix<float> S = S_position_derivative_matrix(system);
	system.s = s_old + timestep_seconds * S * system.u;
	
	// Velocity update
	system.u = u_old + (Minv2 * J2.transpose() * lambda)
			 		 + (timestep_seconds * Minv2 * system.force);	
	normalize_quaternions(system);
}


void integration_step_midpoint_II(struct physics_system& system, float timestep_seconds)
{
	//

	Eigen::VectorXf u_old = Eigen::VectorXf(system.u);
	Eigen::VectorXf s_old = Eigen::VectorXf(system.s);

	integration_step_explicit(system, timestep_seconds/2);

	compute_forces(system);

	contact_list* contacts2 = collision_detectoion(system.colliders, system.s);
	visualise_collisions(system, contacts2);
	std::vector<struct contact_list> contact_table = list_to_array(contacts2);
	//assert(!exists_big_penetration(contact_table));

	Eigen::SparseMatrix<float> Minv2 = generalizedMass_matrix_inv(system);
	Eigen::SparseMatrix<float> J2 = Jmatrix(system, contact_table);
	const int K2 = contact_table.size();

	Eigen::VectorXf lambda = compute_contact_impulses(K2, system, J2, Minv2, timestep_seconds);
	Eigen::SparseMatrix<float> S = S_position_derivative_matrix(system);
	
	system.s = s_old + timestep_seconds * S * system.u;
	system.u = u_old + (Minv2 * J2.transpose() * lambda)
			 		 + (timestep_seconds * Minv2 * system.force);	
	normalize_quaternions(system);
}

void integration_step(struct physics_system& system)
{
	// Symplectic
	// integration_step_symplectic(system, system.base_timestep_seconds);
	// ----------

	// Explicit
	// integration_step_explicit(system, system.base_timestep_seconds);
	// --------

	// Midpoint method
	integration_step_midpoint(system, system.base_timestep_seconds);
	//
}

// blue - no collisions
// yellow - in contact
// red - penetration (bad)
void visualise_collisions(struct physics_system& system, struct contact_list* contacts) {
	const glm::vec3 blue = glm::vec3(0.1, 0.1, 0.6);
	const glm::vec3 yellow = glm::vec3(0.6, 0.6, 0.0);
	const glm::vec3 red = glm::vec3(0.8, 0.1, 0.1);
	for (size_t i = 0; i < system.ridgidbody_count; i++) {
		system.mesh[i].mesh->setSurfaceColor(blue);
	}

	while (contacts != NULL) {
		const glm::vec3 coloring = contacts->penetration ? red : yellow;
		system.mesh[contacts->bodyi_id].mesh->setSurfaceColor(coloring);
		system.mesh[contacts->bodyj_id].mesh->setSurfaceColor(coloring);

		contacts = contacts->next;
	}
}

void physys_render_update(struct physics_system& system) {
	for(size_t i = 0; i < system.ridgidbody_count; i++) {
		// todo
		Eigen::Quaternionf Q(system.s(7*i +3), system.s(7*i +4), system.s(7*i +5), system.s(7*i +6) );
		rendering::update_mesh(system.mesh[i], system.s(Eigen::seq(7*i, 7*i +2)), Q);
	}

	if (NULL != system.spring_visualisation) {
		std::vector<Eigen::Vector3f> nodes;

		for (size_t j = 0; j < system.springs.size();  j ++) {
				struct spring& spring = system.springs[j];
				int body_i = spring.body_i ;
				int body_j = spring.body_j ;

				Eigen::Vector3f r_i_world = ith_Q(system, body_i)._transformVector( spring.r_i );
				Eigen::Vector3f r_j_world = ith_Q(system, body_j)._transformVector( spring.r_j );

				Eigen::Vector3f spring_end_i = ith_x(system, body_i) + r_i_world;
				Eigen::Vector3f spring_end_j = ith_x(system, body_j) + r_j_world;

				nodes.push_back(spring_end_i);
				nodes.push_back(spring_end_j);
		}

		system.spring_visualisation->updateNodePositions(nodes);
	}
}

void start_spring_visualisation(struct physics_system& system) {

	std::vector<Eigen::Vector3f> nodes;
	std::vector<std::array<size_t, 2>> edges;

	for (size_t j = 0; j < system.springs.size();  j ++) {
			nodes.push_back(Eigen::Vector3f::Zero());
			nodes.push_back(Eigen::Vector3f::Zero());
			edges.push_back({2*j, 2*j+1});
	}


	system.spring_visualisation = polyscope::registerCurveNetwork("Spring Curve network", nodes, edges);
	system.spring_visualisation->setColor(glm::vec3(0,0,0));
}


float potential_energy(struct physics_system& system)
{
	float sum = 0.f;

	// E_pp, Gravité
	for (size_t i = 0; i < system.ridgidbody_count; i++ ) {
		if(system.invmass[i] != 0) { // exclude immovable of zero energy
			sum += (0.5/system.invmass[i]) * s_ith_x(system.s, i).y();
		}
	}

	// Ressorts: 
	for (size_t j = 0; j < system.springs.size(); j++) {
		struct spring& spring = system.springs[j];
		int body_i = spring.body_i ;
		int body_j = spring.body_j ;

		Eigen::Vector3f r_i_world = ith_Q(system, body_i)._transformVector( spring.r_i );
		Eigen::Vector3f r_j_world = ith_Q(system, body_j)._transformVector( spring.r_j );

		Eigen::Vector3f spring_end_i = ith_x(system, body_i) + r_i_world;
		Eigen::Vector3f spring_end_j = ith_x(system, body_j) + r_j_world;

		float l = (spring_end_i - spring_end_j).norm();
		sum += 0.5 * spring.k * (l - spring.l0) * (l - spring.l0); 
	}
	return sum;
}

float kinetic_energy(struct physics_system& system)
{
	float sum = 0.f;

	for (size_t i = 0; i < system.ridgidbody_count; i++ ) {
		if(system.invmass[i] != 0) { // exclude immovable of zero energy
			sum += (0.5/system.invmass[i]) * ith_v(system, i).squaredNorm();

			Eigen::Vector3f omega = ith_omega(system, i);
			sum += 0.5 * dot_product(omega, ith_Iinv(system, i).inverse() * omega);
		}
	}
	return sum;
}

