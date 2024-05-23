#ifndef RENDERING_H
# define RENDERING_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <polyscope/polyscope.h>
#include "polyscope/surface_mesh.h"

namespace rendering {


struct mesh{
	Eigen::MatrixXf vertices; // V x 3
	Eigen::MatrixXf faces; // F x 4
	polyscope::SurfaceMesh* mesh;
	float scale;
};

struct mesh create_cube(std::string name);
void update_mesh(rendering::mesh& mesh, Eigen::Vector3f position, Eigen::Quaternionf rotaion);

}

#endif
