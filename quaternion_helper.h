#ifndef QUATERNION_HELPER_H
#define QUATERNION_HELPER_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>

inline Eigen::Quaternionf quaternionf_from_vector_part(Eigen::Vector3f v) {
	return Eigen::Quaternionf(0, v[0], v[1], v[2]);
}

inline Eigen::Quaternionf operator+(Eigen::Quaternionf a, Eigen::Quaternionf b) {
	return Eigen::Quaternionf(a.w() + b.w(), a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}

inline Eigen::Quaternionf q_sclale(float lambda, Eigen::Quaternionf q) {
	return Eigen::Quaternionf(lambda * q.w(),lambda * q.x(), lambda * q.y(), lambda * q.z());
}

#endif
