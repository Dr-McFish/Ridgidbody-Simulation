#ifndef ALGEBRA_UTILS_H
#define ALGEBRA_UTILS_H
#include <eigen3/Eigen/Dense>

inline float dot_product(Eigen::Vector3f a, Eigen::Vector3f b) {
	return a.dot(b);
}

inline Eigen::Vector3f cross_product(Eigen::Vector3f a, Eigen::Vector3f b) {
	return a.cross(b);
}
#endif

