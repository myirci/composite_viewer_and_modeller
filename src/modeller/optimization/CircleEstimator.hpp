#ifndef CIRCLE_ESTIMATOR_HPP
#define CIRCLE_ESTIMATOR_HPP

#include <Eigen/Dense>

class Circle3D;
class Ellipse2D;
class PersProjParam;

class CircleEstimator {

public:
    CircleEstimator() { }
    void estimate_3d_circles_with_fixed_radius_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_radius);
    void estimate_3d_circles_with_fixed_depth_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_depth);
    int estimate_3d_circles_with_fixed_radius(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_radius);
    int estimate_3d_circles_with_fixed_depth(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_depth);
    int estimate_unit_3d_circles(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp);
private:
    void construct_change_of_basis_matrix(Eigen::Matrix3d& mat, const Eigen::Vector3d& vec2);
    bool check_eigenvalue_constraints(const Eigen::Vector3d& eigenvalues);
    void estimate_unit_3d_circles_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp);
};

#endif // CIRCLE_ESTIMATOR_HPP
