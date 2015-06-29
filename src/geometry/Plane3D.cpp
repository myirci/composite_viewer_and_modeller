#include "Plane3D.hpp"

#include <iostream>
#include <cmath>

#include "Ray3D.hpp"
#include "../utility/Utility.hpp"

Plane3D::Plane3D(const osg::Vec4d& vec) : plane(vec) { }
Plane3D::Plane3D(double a, double b, double c, double d) : plane(osg::Vec4d(a,b,c,d)) { }
Plane3D::Plane3D(const osg::Vec3d& normal, const osg::Vec3d& p) { plane = osg::Vec4d(normal.x(), normal.y(), normal.z(), -normal * p); }
Plane3D::Plane3D(const Plane3D& pl) { plane = pl.get_plane(); }
Plane3D::Plane3D(const osg::Vec3d& p1, const osg::Vec3d& p2, const osg::Vec3d& p3) : Plane3D((p2-p1)^(p3-p1), p3) { }
Plane3D::Plane3D(const Eigen::Vector3d& normal, const Eigen::Vector3d& pt) {
    plane.x() = normal[0];
    plane.y() = normal[1];
    plane.z() = normal[2];
    plane.w() = - normal.dot(pt);
}

void Plane3D::get_unit_normal(osg::Vec3d& normal) const {

    normal.x() = plane.x();
    normal.y() = plane.y();
    normal.z() = plane.z();
    normal.normalize();
}

void Plane3D::hessian_normal_form() {
    plane /= std::sqrt(plane.x()*plane.x() + plane.y()*plane.y() + plane.z()*plane.z());
}

void Plane3D::get_normal(osg::Vec3d& normal) const {

    normal.x() = plane.x();
    normal.y() = plane.y();
    normal.z() = plane.z();
}

void Plane3D::print() const { std::cout << "plane: " << plane; }

void Plane3D::translate(const osg::Vec3d& tvec) {

    osg::Vec3d normal;
    get_normal(normal);
    plane.w() += normal * tvec;
}

void Plane3D::rotate(const osg::Matrix3d& mat) {

    osg::Vec3d col1(mat(0,0), mat(1,0), mat(2,0));
    osg::Vec3d col2(mat(0,1), mat(1,1), mat(2,1));
    osg::Vec3d col3(mat(0,2), mat(1,2), mat(2,2));
    osg::Vec3d normal;
    get_normal(normal);
    plane.x() = normal * col1;
    plane.y() = normal * col2;
    plane.z() = normal * col3;
}

void Plane3D::transform(const osg::Matrixd& mat) {

    osg::Matrix3d rot;
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j)
            rot(i,j) = mat(i,j);
    }
    rotate(rot);

    osg::Vec3d tvec;
    tvec[0] = mat(0,3);
    tvec[1] = mat(1,3);
    tvec[2] = mat(2,3);

    translate(tvec);
}

bool Plane3D::is_on_plane(const osg::Vec3d& pt) const {
    double eval = pt.x() * plane.x() + pt.y() * plane.y() + pt.z() * plane.z() + plane.w();
    return  std::abs(eval) < EPSILON;
}
