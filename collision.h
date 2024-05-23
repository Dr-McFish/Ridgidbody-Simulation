#ifndef COLLISION_H
#define COLLISION_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

const float penetration_epsilon = 1e-5;

struct sphere_colider {
	Eigen::Vector3f pos;
	float radius;
};

struct half_space_colider {
	Eigen::Vector3f pos;
	Eigen::Vector3f normal; // unitary vector!
};

// add more when needed
enum colider_type {COLIDER_SPHERE, COLIDER_HALF_SPACE};

struct collider{
	enum colider_type type;
	union{
		struct sphere_colider sphere_colider; 
		struct half_space_colider half_space_colider;
	} u;
};

struct contact_list{
	int solid1_id;
	int solid2_id;
	// invarient: solid1 =/= solid2
	bool penetration;
	Eigen::Vector3f contact_pos;
	Eigen::Vector3f contact_normal;
	float penetration_depth;

	// NULL terminated list
	struct contact_list* next;
};

//borrows an array of colliders and returns a list of contacts
struct contact_list* collision_detectoion(struct collider* coliders);
void copy_collider(struct collider& source, struct collider& dest);

#endif
