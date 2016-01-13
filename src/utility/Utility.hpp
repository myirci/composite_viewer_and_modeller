#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <ostream>

const double PI = 3.141592653589793238462643383279502884197169399375105;
const double TWO_PI = 2 * PI;
const double HALF_PI = PI / 2.0;
const double EPSILON = 0.00000001;

double deg2rad(double deg);
double rad2deg(double rad);

namespace osg { class Vec2d; class Vec3d; class Vec4d; }
double calculate_angle_in_degrees(const osg::Vec3d& v1, const osg::Vec3d& v2);
double calculate_angle_in_radians(const osg::Vec3d& v1, const osg::Vec3d& v2);

std::ostream& operator<<(std::ostream& out, const osg::Vec2d& vec);
std::ostream& operator<<(std::ostream& out, const osg::Vec3d& vec);
std::ostream& operator<<(std::ostream& out, const osg::Vec4d& vec);

class Ray3D;
class Plane3D;
bool is_parallel(const osg::Vec3d& vec1, const osg::Vec3d& vec2);
bool is_parallel(const Plane3D& pl1, const Plane3D& pl2);
bool is_parallel(const Ray3D& ray, const Plane3D& plane);
bool intersect_ray_and_plane(const Ray3D& ray, const Plane3D& plane, osg::Vec3d& intersection);
bool intersect_ray_and_plane_2(const Ray3D& ray, const Plane3D& plane, osg::Vec3d& intersection);
bool is_intersecting(const Ray3D& ray, const Plane3D& plane);
bool intersect_two_planes(const Plane3D& pl1, const Plane3D& pl2, osg::Vec3d& start, osg::Vec3d& dir);

#endif // UTILITY_HPP
