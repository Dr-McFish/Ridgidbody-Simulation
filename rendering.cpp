#include "rendering.h"
#include "glm/fwd.hpp"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <string>

struct rendering::mesh rendering::create_cube(std::string name) {
	struct rendering::mesh cube;

	cube.scale = 0.5;
	cube.vertices = Eigen::MatrixXf(8, 3);
	cube.vertices << -1.0,  1.0,  1.0,
					  1.0,  1.0,  1.0, 
					 -1.0, -1.0,  1.0,
					  1.0, -1.0,  1.0, 
					 -1.0,  1.0, -1.0,
					  1.0,  1.0, -1.0, 
					 -1.0, -1.0, -1.0,
					  1.0, -1.0, -1.0;

	cube.faces = Eigen::MatrixXf(6, 4);
	cube.faces << 1,0,2,3,
				  5,4,6,7,
	              2,0,4,6,
				  3,1,5,7,
	              1,0,4,5,
				  3,2,6,7 ;

	cube.mesh = polyscope::registerSurfaceMesh(name, cube.vertices, cube.faces);
	cube.mesh->setSurfaceColor(glm::vec3(0.1, 0.5, 0.6));
	return cube;
}

struct rendering::mesh rendering::create_square(std::string name) {
	struct rendering::mesh square;

	square.scale = 0.5;
	square.vertices = Eigen::MatrixXf(4, 3);
	square.vertices <<  1.0,  0.0,  1.0,  
					   -1.0,  0.0,  1.0,
					   -1.0,  0.0, -1.0,
					    1.0,  0.0, -1.0, 

	square.faces = Eigen::MatrixXf(1, 4);
	square.faces << 1, 2, 3, 4;

	square.mesh = polyscope::registerSurfaceMesh(name, square.vertices, square.faces);
	square.mesh->setSurfaceColor(glm::vec3(0.1, 0.5, 0.6));
	return square;
}


void rendering::update_mesh(rendering::mesh& mesh, Eigen::Vector3f position, Eigen::Quaternionf rotaion) {
Eigen::MatrixXf new_vertices(8, 3);
		for(int i = 0; i < 8; i++) {
			Eigen::Vector3f vertex(mesh.vertices(i, 0), mesh.vertices(i, 1), mesh.vertices(i, 2));
			Eigen::Vector3f new_vertex = position + mesh.scale* rotaion._transformVector(vertex);

			new_vertices(i, Eigen::all) << new_vertex(0), new_vertex(1), new_vertex(2);
		}
		mesh.mesh->updateVertexPositions(new_vertices);
}