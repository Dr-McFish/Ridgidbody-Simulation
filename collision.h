#ifndef COLLISION_H
#define COLLISION_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

const float penetration_epsilon_m = 1e-2;

struct sphere_colider {
	float radius;
};

struct half_space_colider {
	Eigen::Vector3f normal; // unitary vector!
};

struct convex_hull_collider {
	int number_points;
	
	float** points; // dans le referentiel du solide
	float** points_world; // dans le refferentiell galliean
};

struct convex_hull_collider make_cube_collider();

// add more when needed
enum colider_type {COLIDER_SPHERE, COLIDER_HALF_SPACE, COLIDER_CONVEX_HULL};


struct collider{
	enum colider_type type;
	union U{
		struct sphere_colider sphere_colider; 
		struct half_space_colider half_space_colider;
		struct convex_hull_collider convex_hull_collider;

		// https://stackoverflow.com/questions/70428563/unions-default-constructor-is-implicitly-deleted
		constexpr U() {}
		bool dummy{false};
	} u;

	collider(const collider&);
	collider();
};

struct contact_list{
	int bodyi_id;
	int bodyj_id;
	// invarient: solid1 =/= solid2
	bool penetration;
	Eigen::Vector3f contact_pos;
	Eigen::Vector3f contact_normal; // from 1 toward 2
	Eigen::Vector3f contact_tangent1;
	Eigen::Vector3f contact_tangent2;

	float penetration_depth;

	Eigen::Vector3f contact_velocity;
	float contact_velocity_projected;
	bool colliding_contact; // true if colliding contat, false if sliding contact

	// NULL terminated list
	struct contact_list* next;
};
void free_contact_list(struct contact_list* c);

//borrows an array of colliders and returns a list of contacts
struct contact_list* collision_detectoion(std::vector<struct collider> coliders, Eigen::VectorXf& s);

// destroys the list
std::vector<struct contact_list> list_to_array(struct contact_list* list);

#endif
