#include <cassert>
#include <cmath>

#include "Circle3D.hpp"
#include "../utility/RandomNumberGenerator.hpp"
#include "../utility/Utility.hpp"
#include "../osg/OsgUtility.hpp"

#include <osg/Vec3d>

Circle3D::Circle3D(const Circle3D& rhs) {
    this->center = rhs.center;
    this->normal = rhs.normal;
    this->radius = rhs.radius;
}

Circle3D& Circle3D::operator =(const Circle3D& rhs) {
    if(this != &rhs) {
        this->center = rhs.center;
        this->normal = rhs.normal;
        this->radius = rhs.radius;
    }
    return *this;
}

void Circle3D::generate_data(osg::ref_ptr<osg::Vec3Array>& vertices, osg::ref_ptr<osg::Vec3Array>& normals, int num_points, bool random) const {

    Eigen::Vector3d e1, e2;
    if(random) find_random_orthonomal_basis(e1, e2);
    else find_orthonomal_basis(e1, e2);

    // Construct the data points for the given planar section
    const double step = TWO_PI/static_cast<double>(num_points);
    Eigen::Vector3d n, pos;

    for(double t = 0.0; t < TWO_PI; t += step) {
        n = e1*(radius*cos(t)) + e2*(radius*sin(t));
        pos = center + n;
        vertices->push_back(osg::Vec3(pos[0], pos[1], pos[2]));
        n.normalize();
        normals->push_back(osg::Vec3(n[0], n[1], n[2]));
    }
}

void Circle3D::generate_data(osg::ref_ptr<osg::Vec3Array>& vertices, int num_points, bool random) const {

    Eigen::Vector3d e1, e2;
    if(random) find_random_orthonomal_basis(e1, e2);
    else find_orthonomal_basis(e1, e2);

    const double step = TWO_PI/static_cast<double>(num_points);
    Eigen::Vector3d pos;

    for(double t = 0.0; t < TWO_PI; t+= step) {
        pos = center + (radius*cos(t))*e1 + (radius*sin(t))*e2;
        vertices->push_back(osg::Vec3(pos[0], pos[1], pos[2]));
    }
}

void Circle3D::find_orthonomal_basis(Eigen::Vector3d& e1, Eigen::Vector3d& e2) const {

    // Find 3 othonormal basis vectors for the planar section: 2 vectors
    // in the plane of the section and 1 vector along the direction of the normal
    if(normal[0] != 0 || normal[1] != 0) {
        e1 = Eigen::Vector3d(-normal[1], normal[0], 0);
        e1.normalize();
    }
    else {
        e1 = Eigen::Vector3d(1, 0, 0);
    }

    e2 = normal.cross(e1);

#ifdef DEBUG
    // Check that e1 and e2 is on the plane:
    Eigen::Vector4d plane(normal(0), normal(1), normal(2), - normal.dot(center));
    Eigen::Vector3d v1 = e1 + center;
    Eigen::Vector3d v2 = e2 + center;
    assert(plane.dot(Eigen::Vector4d(v1(0), v1(1), v1(2), 1)) < 0.00001);
    assert(plane.dot(Eigen::Vector4d(v2(0), v2(1), v2(2), 1)) < 0.00001);
#endif // DEBUG

}

void Circle3D::find_random_orthonomal_basis(Eigen::Vector3d& e1, Eigen::Vector3d& e2) const {

    // Find 3 othonormal basis vectors for the planar section: 2 vectors
    // in the plane of the section and 1 vector along the direction of the normal
    RandomNumberGenerator gen;
    gen.initialize_uniform_double_distributor(-2.0, 2.0);
    Eigen::Vector3d random_vec(gen.generate_double(),
                               gen.generate_double(),
                               gen.generate_double());
    e1 = normal.cross(random_vec);
    e1.normalize();
    e2 = normal.cross(e1);
    e2.normalize();
#ifdef DEBUG
    // Check that e1 is not a zero vector plane:
    assert(e1 != Eigen::Vector3d::Zero());
    // Check that e1 and e2 is on the plane:
    Eigen::Vector4d plane(normal(0), normal(1), normal(2), - normal.dot(center));
    Eigen::Vector3d v1 = e1 + center;
    Eigen::Vector3d v2 = e2 + center;
    assert(plane.dot(Eigen::Vector4d(v1(0), v1(1), v1(2), 1)) < 0.00001);
    assert(plane.dot(Eigen::Vector4d(v2(0), v2(1), v2(2), 1)) < 0.00001);
#endif // DEBUG
}

std::ostream& operator<<(std::ostream& out, const Circle3D& circle) {

    out << "center: " << circle.center[0] << "\t" << circle.center[1] << "\t" << circle.center[2] << std::endl;
    out << "normal: " << circle.normal[0] << "\t" << circle.normal[1] << "\t" << circle.normal[2] << std::endl;
    out << "radius: " << circle.radius << std::endl;
    out << "d: " << -(circle.center.dot(circle.normal)) << std::endl;
    return out;
}
