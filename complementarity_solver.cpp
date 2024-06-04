#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>
#include <math.h>
#include "complementarity_solver.h"

float window(float minval, float val, float maxval) {
	return fminf(fmaxf(minval, val), maxval);
}

//returns lambda
Eigen::VectorXf pgs_solve(struct linear_complementarity_problem* problem, int iterations) {
	const int system_size = problem->lamba_max.size();

	Eigen::VectorXf lambda = Eigen::VectorXf(system_size);
	// auto A_Lower = problem->A.triangularView<Eigen::StrictlyLower>();
	// auto A_Upper = problem->A.triangularView<Eigen::StrictlyUpper>();

	for (int iteration = 0; iteration < iterations; iteration++) {
		// Gauss-Siedel
		for (int i = 0; i < system_size; i++) {
			lambda(i) = -problem->b(i);
			//sum
			for(int j = 0; j < i; j++) {
				lambda(i) -= problem->A.coeff(i, j) * lambda(j);
			}
			for(int j = i+1; j < system_size; j++) {
				lambda(i) -= problem->A.coeff(i, j) * lambda(j);
			}

			lambda(i) = lambda(i) / problem->A.coeff(i, i);
	
		}
		// Projection
		for (int i = 0; i < system_size; i++) {
			lambda(i) = window(problem->lamba_min(i), lambda(i), problem->lamba_max(i));
		}

	}
	
	return lambda;
}