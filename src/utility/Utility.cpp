#include "Utility.hpp"
#include <osg/Vec3d>

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
