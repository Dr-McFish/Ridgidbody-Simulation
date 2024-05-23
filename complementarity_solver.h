#ifndef COMPLEMENTARITY_SOLVER_H
#define COMPLEMENTARITY_SOLVER_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>

struct linear_complementarity_problem {
	// Aλ + b >= 0
	// λ_min <= λ <= λ_max

	//invarient: A is n x n
	//			 λ, λ_min, λ_max are n x 1
	//           b is n x 1

	// ussually n = 3K 

	Eigen::MatrixXf A;

	Eigen::VectorXf b;
	Eigen::VectorXf lamba_min;
	Eigen::VectorXf lamba_max;
};

//returns lambda
Eigen::VectorXf pgs_solve(struct linear_complementarity_problem* problem, int iterations);

#endif

