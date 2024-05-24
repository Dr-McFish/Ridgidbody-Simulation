#ifndef COLLISION_H
#define COLLISION_H

#include <cstring>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

const float penetration_epsilon_m = 1e-2;

struct sphere_colider {
	float radius;
};

struct half_space_colider {
	Eigen::Vector3f normal; // unitary vector!
};

// add more when needed
enum colider_type {COLIDER_SPHERE, COLIDER_HALF_SPACE};


struct collider{
	Eigen::Vector3f* pos;
	enum colider_type type;
	union U{
		struct sphere_colider sphere_colider; 
		struct half_space_colider half_space_colider;

		// https://stackoverflow.com/questions/70428563/unions-default-constructor-is-implicitly-deleted
		constexpr U() {}
		bool dummy{false};
	} u;
};

struct contact_list{
	int solid1_id;
	int solid2_id;
	// invarient: solid1 =/= solid2
	bool penetration;
	Eigen::Vector3f contact_pos;
	Eigen::Vector3f contact_normal; // from 1 toward 2
	float penetration_depth;

	Eigen::Vector3f contact_velocity;
	bool colliding_contact; // true if colliding contat, false if sliding contact

	// NULL terminated list
	struct contact_list* next;
};
void free_contact_list(struct contact_list* c);

//borrows an array of colliders and returns a list of contacts
struct contact_list* collision_detectoion(int count, struct collider* coliders);
void copy_collider(struct collider* source, struct collider* dest);

#endif
