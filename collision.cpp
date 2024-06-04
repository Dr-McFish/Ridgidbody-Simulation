#include "collision.h"
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include "algebra_utils.h"


struct contact_list* sphere_sphere_collision_detectoion(struct sphere_colider& sphere1, Eigen::Vector3f sphere1_pos,
														struct sphere_colider& sphere2, Eigen::Vector3f sphere2_pos);

struct contact_list* sphere_half_space_collision_detectoion(struct sphere_colider& sphere,		  	Eigen::Vector3f sphere_pos, 
															struct half_space_colider& half_space,	Eigen::Vector3f half_space_pos);

void free_contact_list(struct contact_list* c) {
	if(c != NULL) {
		free_contact_list(c->next);
		free(c);
	}
}

//borrows a list of colliders and returns a list of contacts
struct contact_list* collision_detectoion(int count, struct collider* colliders) {
	struct contact_list* ret = NULL;

	auto add_collision_to_list = [&](int idx1, int idx2, struct contact_list* collision_maybe) {
		if(collision_maybe != NULL) {
			collision_maybe->bodyi_id = idx1;
			collision_maybe->bodyj_id = idx2;
			
			collision_maybe->next = ret;
			ret = collision_maybe;
		}
	};

	// Very naive O(n^2)
	for(int i =0; i < count; i++) {
		for(int j =0; j < i; j++) {
			struct contact_list* collision_maybe = NULL;

			if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_SPHERE)
			{
				collision_maybe = sphere_sphere_collision_detectoion(colliders[i].u.sphere_colider, *colliders[i].pos,
																	 colliders[j].u.sphere_colider, *colliders[j].pos);
				add_collision_to_list(i,j, collision_maybe);
			} 
			else if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[i].u.sphere_colider,		*colliders[i].pos,
																		 colliders[j].u.half_space_colider,*colliders[j].pos);
				add_collision_to_list(i,j, collision_maybe);
			}
			else if(colliders[j].type == COLIDER_SPHERE && colliders[i].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[j].u.sphere_colider,		*colliders[j].pos,
																		 colliders[i].u.half_space_colider,*colliders[i].pos);
				// NOTE the inversion of j and i
				add_collision_to_list(j,i, collision_maybe); 
			} 
			// else detection not supported.

			if (NULL != collision_maybe) {
				Eigen::Matrix3f rotaion_90x;
				rotaion_90x <<	1, 0, 0,
								0, 0, 1,
								0,-1, 0;
				
				collision_maybe->contact_tangent1 = rotaion_90x * collision_maybe->contact_normal;
				collision_maybe->contact_tangent1.normalize();

				collision_maybe->contact_tangent2 = cross_product(collision_maybe->contact_normal, collision_maybe->contact_tangent1)
											.normalized();
			}
		}
	}

	return ret;
}


// returns NULL if there is no contact,
// 		   a single ellement list in case of contactr
// /!\ solid1_id and solid2_id unintialized
struct contact_list* sphere_sphere_collision_detectoion(struct sphere_colider& sphere1, Eigen::Vector3f sphere1_pos,
														struct sphere_colider& sphere2, Eigen::Vector3f sphere2_pos) {
	const float collision_distance = sphere1.radius + sphere2.radius;

	if( (sphere1_pos - sphere2_pos).squaredNorm() > collision_distance*collision_distance) {
		// No collision
		return NULL;
	}
	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	//avrege
	ret->contact_pos = (1.0 / (sphere1.radius + sphere2.radius))*(sphere1_pos*sphere1.radius  + sphere2_pos*sphere2.radius);
	ret->contact_normal = sphere1_pos - sphere2_pos;
	ret->contact_normal.normalize();
	ret->penetration_depth = fabsf(collision_distance) - (sphere1_pos - sphere2_pos).norm();
	assert(ret->penetration_depth >= 0.);
	ret->penetration = ret->penetration_depth > penetration_epsilon_m;
	
	// ret->solid1_id;
	// ret->solid2_id;
	// ret->tangent1;
	// ret->tangent2;

	return ret; 
}

// returns NULL if there is no contact,
// 		   a single ellement list in case of contactr
// /!\ solid1_id and solid2_id unintialized
struct contact_list* sphere_half_space_collision_detectoion(struct sphere_colider& sphere,		  	Eigen::Vector3f sphere_pos, 
															struct half_space_colider& half_space,	Eigen::Vector3f half_space_pos){
	const Eigen::Vector3f relative_pos = sphere_pos - half_space_pos;
	const float normal_projection = dot_product(relative_pos, half_space.normal);
	const float penetration_depth = normal_projection - sphere.radius;

	if ( penetration_depth > penetration_epsilon_m) {
		// No collision
		return NULL;
	}
	// std::cout << sphere_pos << std::endl;
	// std::cout << half_space_pos << std::endl;
	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	ret->contact_normal = half_space.normal;
	ret->contact_pos = sphere_pos - sphere.radius*half_space.normal;

	ret->penetration_depth = - penetration_depth;	
	ret->penetration = penetration_depth < -penetration_epsilon_m;
	// invarient: penetration_depth \in [-penetration_epsilon, +penetration_epsilon]

	// ret->solid1_id;
	// ret->solid2_id;
	// ret->tangent1;
	// ret->tangent2;
	return ret; 
}

void copy_collider(struct collider* source, struct collider* dest) {
	memcpy(dest, source, sizeof(struct collider));
}

int list_length(struct contact_list* list) {
	if(NULL == list) {
		return 0;
	} else {
		return 1 + list_length(list->next);
	}
}

void parcours_coppie(struct contact_list* list, struct contact_list* array, int pos) {
	if(pos < 0) {
		return;
	} else {
		assert(NULL != list);
		array[pos] = *list;
		array[pos].next = NULL;
		parcours_coppie(list->next, array, pos - 1);
		free(list);
	}
}

struct contact_list* list_to_array(struct contact_list* list, int* adr_taille) {
	*adr_taille = list_length(list);

	struct contact_list* array = (struct contact_list*)calloc(*adr_taille, sizeof(struct contact_list));
	parcours_coppie(list, array, *adr_taille -1);
	
	return array;
}

