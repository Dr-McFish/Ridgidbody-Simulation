#include "collision.h"
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <openGJK/openGJK.h>
#include "algebra_utils.h"
#include "physics_engine.h"

struct convex_hull_collider make_cube_collider() {
	struct convex_hull_collider ret;

	ret.number_points = 8;

	ret.points = (float**)calloc(ret.number_points, sizeof(float*));
	ret.points_world = (float**)calloc(ret.number_points, sizeof(float*));

	for (int i = 0; i < ret.number_points; i++) {
		ret.points[i] = (float*)calloc(3, sizeof(float));
		ret.points_world[i] = (float*)calloc(3, sizeof(float));
	}

	ret.points[0][0] = +0.5f; ret.points[0][1] = +0.5f; ret.points[0][2] = +0.5f;
	ret.points[1][0] = +0.5f; ret.points[1][1] = +0.5f; ret.points[1][2] = -0.5f;
	ret.points[2][0] = +0.5f; ret.points[2][1] = -0.5f; ret.points[2][2] = +0.5f;
	ret.points[3][0] = +0.5f; ret.points[3][1] = -0.5f; ret.points[3][2] = -0.5f;
	ret.points[4][0] = -0.5f; ret.points[4][1] = +0.5f; ret.points[4][2] = +0.5f;
	ret.points[5][0] = -0.5f; ret.points[5][1] = +0.5f; ret.points[5][2] = -0.5f;
	ret.points[6][0] = -0.5f; ret.points[6][1] = -0.5f; ret.points[6][2] = +0.5f;
	ret.points[7][0] = -0.5f; ret.points[7][1] = -0.5f; ret.points[7][2] = -0.5f;

	return ret;
}

void update_relative_positions(struct convex_hull_collider& hull, Eigen::Vector3f translation, Eigen::Quaternionf rotaion){
	for (int i = 0; i < hull.number_points; i++) {
		Eigen::Vector3f original = Eigen::Vector3f(hull.points[i][0], hull.points[i][1], hull.points[i][2]);
		Eigen::Vector3f transformed = translation + rotaion._transformVector(original);

		hull.points_world[i][0] = transformed.x();
		hull.points_world[i][1] = transformed.y();
		hull.points_world[i][2] = transformed.z();
	}
}


struct contact_list* sphere_sphere_collision_detectoion(struct sphere_colider& sphere1, Eigen::Vector3f sphere1_pos,
														struct sphere_colider& sphere2, Eigen::Vector3f sphere2_pos);

struct contact_list* sphere_half_space_collision_detectoion(struct sphere_colider& sphere,		  	Eigen::Vector3f sphere_pos, 
															struct half_space_colider& half_space,	Eigen::Vector3f half_space_pos);

struct contact_list* sphere_convex_hull_collision_detectoion(struct sphere_colider& sphere, Eigen::Vector3f sphere_pos,
														struct convex_hull_collider& hull, Eigen::Vector3f hull_pos, Eigen::Quaternionf hull_Q);

struct contact_list* half_space_convex_hull_collision_detectoion(struct half_space_colider& half_space,	Eigen::Vector3f half_space_pos,
														struct convex_hull_collider& hull, Eigen::Vector3f hull_pos, Eigen::Quaternionf hull_Q);

struct contact_list* convex_hull_convex_hull_collision_detectoion
													(struct convex_hull_collider& hull1, Eigen::Vector3f hull_pos1, Eigen::Quaternionf hull_Q1,
													 struct convex_hull_collider& hull2, Eigen::Vector3f hull_pos2, Eigen::Quaternionf hull_Q2);



void free_contact_list(struct contact_list* c) {
	if(c != NULL) {
		free_contact_list(c->next);
		free(c);
	}
}
struct contact_list* collision_detectoion(std::vector<struct collider> colliders, Eigen::VectorXf& s) {
	int count = colliders.size();
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
				collision_maybe = sphere_sphere_collision_detectoion(colliders[i].u.sphere_colider, s_ith_x(s, i),
																	 colliders[j].u.sphere_colider, s_ith_x(s, j));
				add_collision_to_list(i,j, collision_maybe);
			} 
			else if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[i].u.sphere_colider,		s_ith_x(s, i),
																		 colliders[j].u.half_space_colider,s_ith_x(s, j));
				add_collision_to_list(i,j, collision_maybe);
			}
			else if(colliders[j].type == COLIDER_SPHERE && colliders[i].type == COLIDER_HALF_SPACE)
			{
				collision_maybe = sphere_half_space_collision_detectoion(colliders[j].u.sphere_colider,		s_ith_x(s, j),
																		 colliders[i].u.half_space_colider,s_ith_x(s, i));
				// NOTE the inversion of j and i
				add_collision_to_list(j,i, collision_maybe); 
			}
			else if(colliders[i].type == COLIDER_SPHERE && colliders[j].type == COLIDER_CONVEX_HULL) {
				collision_maybe = sphere_convex_hull_collision_detectoion(colliders[i].u.sphere_colider,     s_ith_x(s, i),
																		  colliders[j].u.convex_hull_collider, s_ith_x(s, j), s_ith_Q(s, j) );
				add_collision_to_list(i, j, collision_maybe); 
			}
			else if(colliders[j].type == COLIDER_SPHERE && colliders[i].type == COLIDER_CONVEX_HULL) {
				collision_maybe = sphere_convex_hull_collision_detectoion(colliders[j].u.sphere_colider,     s_ith_x(s, j),
																		  colliders[i].u.convex_hull_collider, s_ith_x(s, i), s_ith_Q(s, i) );

				add_collision_to_list(j,i, collision_maybe); 
			}
			else if(colliders[i].type == COLIDER_HALF_SPACE && colliders[j].type == COLIDER_CONVEX_HULL) {
				collision_maybe = half_space_convex_hull_collision_detectoion(colliders[i].u.half_space_colider,     s_ith_x(s, i),
																		  colliders[j].u.convex_hull_collider, s_ith_x(s, j), s_ith_Q(s, j) );
				add_collision_to_list(i, j, collision_maybe); 
			}
			else if(colliders[j].type == COLIDER_HALF_SPACE && colliders[i].type == COLIDER_CONVEX_HULL) {
				collision_maybe = half_space_convex_hull_collision_detectoion(colliders[j].u.half_space_colider,     s_ith_x(s, j),
																		  colliders[i].u.convex_hull_collider, s_ith_x(s, i), s_ith_Q(s, i) );
				add_collision_to_list(j, i, collision_maybe); 
			}
			else if(colliders[i].type == COLIDER_CONVEX_HULL && colliders[j].type == COLIDER_HALF_SPACE) {
				collision_maybe = convex_hull_convex_hull_collision_detectoion(colliders[i].u.convex_hull_collider, s_ith_x(s, i), s_ith_Q(s, i),
																			   colliders[j].u.convex_hull_collider, s_ith_x(s, j), s_ith_Q(s, j));
				add_collision_to_list(i, j, collision_maybe); 
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
	//assert(ret->penetration_depth >= 0.);
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

	//assert(false);
	return ret; 
}


struct contact_list* sphere_convex_hull_collision_detectoion(struct sphere_colider& sphere, Eigen::Vector3f sphere_pos,
														 struct convex_hull_collider& hull, Eigen::Vector3f hull_pos, Eigen::Quaternionf hull_Q)
{
	update_relative_positions(hull, hull_pos, hull_Q);
	float sphere_center[3] = {sphere_pos.x(), sphere_pos.y(), sphere_pos.z()};
	float* _single_point = sphere_center;
	float** single_point = &_single_point;

	gkSimplex s;
	/* Initialise simplex as empty */
	s.nvrtx = 0;

	gkPolytope sphere_point = {.numpoints = 1, .coord = single_point};
	gkPolytope hull_tope = {.numpoints = hull.number_points, .coord = hull.points_world};

	float distance = 5.f; //compute_minimum_distance(sphere_point, hull_tope, &s);

	const float penetration_depth = distance - sphere.radius;
	if ( penetration_depth > penetration_epsilon_m) {
		// No collision
		return NULL;
	}

	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	assert(false); // no way to get contact normal or pos.
	ret->contact_normal;
	ret->contact_pos;

	ret->penetration_depth = - penetration_depth;	
	ret->penetration = penetration_depth < -penetration_epsilon_m;
	// invarient: penetration_depth \in [-penetration_epsilon, +penetration_epsilon]

	// ret->solid1_id;
	// ret->solid2_id;
	// ret->tangent1;
	// ret->tangent2;

	//assert(false);
	return ret; 
}


struct contact_list* half_space_convex_hull_collision_detectoion(struct half_space_colider& half_space,	Eigen::Vector3f half_space_pos,
														struct convex_hull_collider& hull, Eigen::Vector3f hull_pos, Eigen::Quaternionf hull_Q)
{
	update_relative_positions(hull, hull_pos, hull_Q);
	float min_distance_to_plane = +INFINITY;
	int min_idx;

	for (int i = 0; i < hull.number_points; i++) {
		Eigen::Vector3f pos = Eigen::Vector3f(hull.points_world[i][0], hull.points_world[i][1], hull.points_world[i][2]);
		float distance = dot_product(pos - half_space_pos, half_space.normal);
		if (min_distance_to_plane > distance) {
			min_distance_to_plane = distance;
			min_idx = i;
		}
	}
	const float penetration_depth = min_distance_to_plane;
	if ( penetration_depth > penetration_epsilon_m) {
		// No collision
		return NULL;
	}
	struct contact_list* ret = (struct contact_list*)malloc(sizeof(struct contact_list));

	ret->contact_normal = -half_space.normal;
	ret->contact_pos = Eigen::Vector3f(hull.points_world[min_idx][0], hull.points_world[min_idx][1], hull.points_world[min_idx][2]);

	ret->penetration_depth = - penetration_depth;	
	ret->penetration = penetration_depth < -penetration_epsilon_m;
	// invarient: penetration_depth \in [-penetration_epsilon, +penetration_epsilon]

	// ret->solid1_id;
	// ret->solid2_id;
	// ret->tangent1;
	// ret->tangent2;

	//assert(false);
	return ret;
}

struct contact_list* convex_hull_convex_hull_collision_detectoion
													(struct convex_hull_collider& hull1, Eigen::Vector3f hull_pos1, Eigen::Quaternionf hull_Q1,
													 struct convex_hull_collider& hull2, Eigen::Vector3f hull_pos2, Eigen::Quaternionf hull_Q2)
{
	update_relative_positions(hull2, hull_pos2, hull_Q2);
	update_relative_positions(hull1, hull_pos1, hull_Q1);
	//
	return NULL;
}


int list_length(struct contact_list* list) {
	if(NULL == list) {
		return 0;
	} else {
		return 1 + list_length(list->next);
	}
}

void parcours_coppie(struct contact_list* list, std::vector<struct contact_list>& array) {
	if(NULL == list) {
		return;
	} else {
		assert(NULL != list);
		array.push_back(*list);
		array[array.size() - 1].next = NULL;
		parcours_coppie(list->next, array);
		free(list);
	}
}

std::vector<struct contact_list> list_to_array(struct contact_list* list)
{
	int taille = list_length(list);

	std::vector<struct contact_list> array;
	array.reserve(taille);
	
	parcours_coppie(list, array);
	
	return array;
}

collider::collider(const collider& ref) {
	this->type = ref.type;

	switch (this->type) {
		case COLIDER_SPHERE:
			this->u.sphere_colider = ref.u.sphere_colider;
			break;

		case COLIDER_HALF_SPACE:
			this->u.half_space_colider = ref.u.half_space_colider;
			break;
		case COLIDER_CONVEX_HULL:
			this->u.convex_hull_collider = ref.u.convex_hull_collider;
			break;
		default:
			assert(false);	
	}
}
collider::collider() {
	this->u.dummy = false;
}