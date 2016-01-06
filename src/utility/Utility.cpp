#include <iostream>
#include "Utility.hpp"

#include "../geometry/Ray3D.hpp"
#include "../geometry/Plane3D.hpp"

#include <Eigen/Dense>
#include <osg/Vec2d>
#include <osg/Vec3d>
#include <osg/Vec4d>

double deg2rad(double deg) {

    return (PI/180.0)*deg;
}

double rad2deg(double rad) {

    return (180.0/PI)*rad;
}

double calculate_angle_in_degrees(const osg::Vec3d& v1, const osg::Vec3d& v2) {

    return rad2deg(acos((v1 * v2) / (v1.length() * v2.length())));
}

double calculate_angle_in_radians(const osg::Vec3d& v1, const osg::Vec3d& v2) {

    return acos((v1 * v2) / (v1.length() * v2.length()));
}

std::ostream& operator<<(std::ostream& out, const osg::Vec2d& vec) {

    out << vec.x() << " " << vec.y();
    return out;
}

std::ostream& operator<<(std::ostream& out, const osg::Vec3d& vec) {

    out << vec.x() << " " << vec.y() << " " << vec.z();
    return out;
}

std::ostream& operator<<(std::ostream& out, const osg::Vec4d& vec) {

    out << vec.x() << " " << vec.y() << " " << vec.z() << " " << vec.w();
    return out;
}

bool is_parallel(const osg::Vec3d& vec1, const osg::Vec3d& vec2) {

    return (vec1 ^ vec2).length() < EPSILON ;
}

bool is_parallel(const Ray3D& ray, const Plane3D& plane) {

    osg::Vec3d normal;
    plane.get_unit_normal(normal);
    return std::abs(ray.direction_vector() * normal) < EPSILON;
}

bool is_parallel(const Plane3D& pl1, const Plane3D& pl2) {

    osg::Vec3d n1, n2;
    pl1.get_unit_normal(n1);
    pl2.get_unit_normal(n2);
    return is_parallel(n1, n2);
}

bool intersect_ray_and_plane(const Ray3D& ray, const Plane3D& plane, osg::Vec3d& intersection) {

    if(is_parallel(ray, plane)) {
        std::cout << "Ray and plane are parallel" << std::endl;
        return false;
    }
    else {
        osg::Vec3d normal;
        plane.get_normal(normal);
        double t = -(ray.start_point() * normal + plane.get_plane().w()) / (ray.direction_vector() * normal);

        if(t < 0) {
            std::cout << "Ray intersects with the plane in the reverse direction" << " t: " << t << std::endl;
            return false;
        }
        else if(t == 0) {
            std::cout << "Ray starting point is on the plane" << std::endl;
            intersection = ray.start_point();
        }
        else {
            intersection = ray.start_point() + ray.direction_vector() * t;
        }
    }
    return true;
}

bool is_intersecting(const Ray3D& ray, const Plane3D& plane) {

    if(is_parallel(ray, plane)) return false;
    else {
        osg::Vec3d normal;
        plane.get_normal(normal);
        double t = -(ray.start_point() * normal + plane.get_plane().w()) / (ray.direction_vector() * normal);
        if(t < 0) return false;
    }
    return true;
}

bool intersect_two_planes(const Plane3D& pl1, const Plane3D& pl2, osg::Vec3d& start, osg::Vec3d& dir) {

    if(is_parallel(pl1, pl2)) {
        std::cerr << "Planes are parallel!" << std::endl;
        return false;
    }

    Plane3D plane1(pl1);
    plane1.hessian_normal_form();
    Plane3D plane2(pl2);
    plane2.hessian_normal_form();

    Eigen::Matrix<double, 2, 3> M;
    M << plane1.get_plane().x(), plane1.get_plane().y(), plane1.get_plane().z(),
         plane2.get_plane().x(), plane2.get_plane().y(), plane2.get_plane().z();

    Eigen::Vector2d b;
    b << -plane1.get_plane().w(), -plane2.get_plane().w();
    Eigen::FullPivLU<Eigen::Matrix<double, 2, 3> > lu_decomp(M);
    Eigen::Vector3d kernel = lu_decomp.kernel();
    Eigen::Vector3d sol = lu_decomp.solve(b);
    start.x() = sol[0];
    start.y() = sol[1];
    start.z() = sol[2];
    dir.x() = kernel[0];
    dir.y() = kernel[1];
    dir.z() = kernel[2];
}
