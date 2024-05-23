#include "collision.h"
#include <assert.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include "algebra_utils.h"

struct contact_list* sphere_sphere_collision_detectoion(struct sphere_colider sphere1, struct sphere_colider sphere2);
struct contact_list* sphere_half_space_collision_detectoion(struct sphere_colider sphere, struct half_space_colider half_space);

//borrows a list of colliders and returns a list of contacts
struct contact_list* collision_detectoion(int count, struct collider* colliders) {
	struct contact_list* ret = NULL;

	// Very naive O(n^2)
	for(int i =0; i < count; i++) {
		for(int j =0; j < i; j++) {
			struct contact_list* collision_maybe = NULL;

			if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_SPHERE)
			{
				collision_maybe = sphere_sphere_collision_detectoion(colliders[i].u.sphere_colider, colliders[j].u.sphere_colider);
			} 
			else if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[i].u.sphere_colider, colliders[j].u.half_space_colider);
			}
			else if(colliders[j].type == COLIDER_SPHERE && colliders[i].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[j].u.sphere_colider, colliders[i].u.half_space_colider);
			} // else detection not supported.

			// add collision to list
			if(collision_maybe != NULL) {
				collision_maybe->next = ret;
				ret = collision_maybe;
			}
		}
	}

	return ret;
}


// returns NULL if there is no contact,
// 		   a single ellement list in case of contactr
// /!\ solid1_id and solid2_id unintialized
struct contact_list* sphere_sphere_collision_detectoion(struct sphere_colider sphere1, struct sphere_colider sphere2) {
	const float collision_distance = sphere1.radius - sphere2.radius;

	if( (sphere1.pos - sphere2.pos).squaredNorm() > collision_distance*collision_distance) {
		// No collision
		return NULL;
	}
	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	//avrege
	ret->contact_pos = (1.0 / (sphere1.radius + sphere2.radius))*(sphere1.pos*sphere1.radius  + sphere2.pos*sphere2.radius);
	ret->contact_normal = sphere1.pos - sphere2.pos;
	ret->contact_normal.normalize();
	ret->penetration_depth = fabsf(collision_distance) - (sphere1.pos - sphere2.pos).norm();
	assert(ret->penetration_depth >= 0.);
	ret->penetration = ret->penetration_depth > penetration_epsilon;
	
	// ret->solid1_id;
	// ret->solid2_id;

	return ret; 
}

// returns NULL if there is no contact,
// 		   a single ellement list in case of contactr
// /!\ solid1_id and solid2_id unintialized
struct contact_list* sphere_half_space_collision_detectoion(struct sphere_colider sphere, struct half_space_colider half_space) {
	const Eigen::Vector3f relative_pos = sphere.pos - half_space.pos;
	const float normal_projection = dot_product(relative_pos, half_space.normal);


	if ( normal_projection > 0 ) {
		// No collision
		return NULL;
	}
	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	ret->contact_normal = half_space.normal;
	ret->contact_pos = sphere.pos - sphere.radius*half_space.normal;

	ret->penetration_depth = -normal_projection;	
	ret->penetration = ret->penetration_depth > penetration_epsilon;
	// ret->solid1_id;
	// ret->solid2_id;

	return ret; 
}

void copy_collider(struct collider& source, struct collider& dest) {
	switch(source.type) {
	  
	  case COLIDER_SPHERE:
		dest.u.sphere_colider = source.u.sphere_colider;
		break;
	  case  COLIDER_HALF_SPACE:
		dest.u.sphere_colider = source.u.sphere_colider;
		break;
	  default:
		assert(false); // NOT IMPLEMENTED
		break;
	}
}