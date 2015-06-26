#ifndef CIRCLE3D_HPP
#define CIRCLE3D_HPP

#include <Eigen/Dense>
#include <osg/Array>
#include <vector>

class Circle3D {
public:

    double radius;
    Eigen::Vector3d center;
    Eigen::Vector3d normal;

    // constructors:
    Circle3D() : radius(0.0), center(Eigen::Vector3d::Zero()), normal(Eigen::Vector3d::Zero()) { }
    Circle3D(const Eigen::Vector3d& c, const Eigen::Vector3d& n, double r): center(c), normal(n.normalized()), radius(r) { }

    // copy constructor
    Circle3D(const Circle3D& rhs);

    //assignment operator
    Circle3D& operator=(const Circle3D& rhs);

    // data generators
    void generate_data(osg::ref_ptr<osg::Vec3Array>& vertices, osg::ref_ptr<osg::Vec3Array>& normals, int num_points, bool random = false) const;
    void generate_data(osg::ref_ptr<osg::Vec3Array>& vertices, int num_points, bool random = false) const;
    void generate_aligned_data(osg::ref_ptr<osg::Vec3Array>& vertices, osg::ref_ptr<osg::Vec3Array>& normals, int num_points, const Circle3D& circle) const;

    // A 3D circle can be represented by a dual quadrics
    void get_matrix_representation(Eigen::Matrix4d& mat) const;

    // Scaled normal
    void get_scaled_normal(Eigen::Vector3d& vec) const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    void find_orthonomal_basis(Eigen::Vector3d& e1, Eigen::Vector3d& e2) const;
    void find_random_orthonomal_basis(Eigen::Vector3d& e1, Eigen::Vector3d& e2) const;
};

std::ostream& operator<<(std::ostream& out, const Circle3D& circle);

#endif // CIRCLE3D_HPP
