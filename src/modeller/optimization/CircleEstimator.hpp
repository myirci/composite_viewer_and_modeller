#ifndef CIRCLE_ESTIMATOR_HPP
#define CIRCLE_ESTIMATOR_HPP

#include <Eigen/Dense>

class Circle3D;
class Ellipse2D;
class PersProjParam;

class CircleEstimator {

public:
    CircleEstimator() { }
    int estimate_3d_circles_with_fixed_radius(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_radius);
    int estimate_3d_circles_with_fixed_depth(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_depth);
    int estimate_unit_3d_circles(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp);
    void estimate_3d_circles_under_orthographic_projection(const Ellipse2D& ellipse, Circle3D& circle, double near);
    void estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(const Ellipse2D& ellipse, Circle3D& circle, double near, double desired_depth);
    void estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(const Ellipse2D& ellipse, Circle3D& circle, double near, double desired_depth);
    // void estimate_3d_circles_with_fixed_radius_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_radius);
    // void estimate_3d_circles_with_fixed_depth_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_depth);
private:
    bool check_eigenvalue_constraints(const Eigen::Vector3d& eigenvalues);
    // void estimate_unit_3d_circles_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp);
    // void construct_change_of_basis_matrix(Eigen::Matrix3d& mat, const Eigen::Vector3d& vec2);
};

#endif // CIRCLE_ESTIMATOR_HPP
