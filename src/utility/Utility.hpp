#ifndef UTILITY_HPP
#define UTILITY_HPP

const double PI = 3.141592653589793238462643383279502884197169399375105;
const double TWO_PI = 2 * PI;
const double HALF_PI = PI / 2.0;

double deg2rad(double deg);
double rad2deg(double rad);

class Ellipse2D;
void convert_ellipse_from_logical_device_coordinates_to_projected_point_coordinates(int w, int h, double n, double half_fovy, const Ellipse2D& elp_log, Ellipse2D& elp_prj);

#endif // UTILITY_HPP
