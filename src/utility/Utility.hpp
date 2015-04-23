#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <ostream>

const double PI = 3.141592653589793238462643383279502884197169399375105;
const double TWO_PI = 2 * PI;
const double HALF_PI = PI / 2.0;

double deg2rad(double deg);
double rad2deg(double rad);

namespace osg { class Vec2d; class Vec3d; }
double calculate_angle_in_degrees(const osg::Vec3d& v1, const osg::Vec3d& v2);
double calculate_angle_in_radians(const osg::Vec3d& v1, const osg::Vec3d& v2);

std::ostream& operator<<(std::ostream& out, const osg::Vec2d& vec);
std::ostream& operator<<(std::ostream& out, const osg::Vec3d& vec);


#endif // UTILITY_HPP
