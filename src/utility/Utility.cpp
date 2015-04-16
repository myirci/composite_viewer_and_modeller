#include "Utility.hpp"
#include "../geometry/Ellipse2D.hpp"

double deg2rad(double deg) { return (PI/180.0)*deg; }
double rad2deg(double rad) { return (180.0/PI)*rad; }

void convert_ellipse_from_logical_device_coordinates_to_projected_point_coordinates(
        int w, int h, double n, double half_fovy, const Ellipse2D& elp_log, Ellipse2D& elp_prj) {

    // if coefficient array has not been defined before, allocate memory
    if(elp_prj.coeff == nullptr)
        elp_prj.coeff = new double[6];

    // rot angle (angle between the major axis and the positive x-axis) does not change
    elp_prj.rot_angle = elp_log.rot_angle;

    // find normalized device coordinates
    double consx = 2.0/static_cast<double>(w);
    double consy = 2.0/static_cast<double>(h);
    osg::Vec2d points_ndc[4];
    points_ndc[0] = osg::Vec2d(consx * elp_log.points[0].x() - 1, consy * elp_log.points[0].y() - 1);
    points_ndc[1] = osg::Vec2d(consx * elp_log.points[1].x() - 1, consy * elp_log.points[1].y() - 1);
    points_ndc[2] = osg::Vec2d(consx * elp_log.points[2].x() - 1, consy * elp_log.points[2].y() - 1);
    points_ndc[3] = osg::Vec2d(consx * elp_log.points[3].x() - 1, consy * elp_log.points[3].y() - 1);

    // find projected point coordinates
    double aspect = static_cast<double>(w)/static_cast<double>(h);
    double c1 = n * tan(half_fovy);
    double c2 = aspect * c1;
    elp_prj.points[0] = osg::Vec2d(points_ndc[0].x()*c2, points_ndc[0].y()*c1);
    elp_prj.points[1] = osg::Vec2d(points_ndc[1].x()*c2, points_ndc[1].y()*c1);
    elp_prj.points[2] = osg::Vec2d(points_ndc[2].x()*c2, points_ndc[2].y()*c1);
    elp_prj.points[3] = osg::Vec2d(points_ndc[3].x()*c2, points_ndc[3].y()*c1);

    // calculate the remaining parameters and the coefficients
    elp_prj.center.x() = (elp_prj.points[0].x() + elp_prj.points[1].x())/2.0;
    elp_prj.center.y() = (elp_prj.points[0].y() + elp_prj.points[1].y())/2.0;
    elp_prj.smj_axis = (elp_prj.points[0] - elp_prj.center).length();
    elp_prj.smn_axis = (elp_prj.points[2] - elp_prj.center).length();
    double as = elp_prj.smj_axis * elp_prj.smj_axis;
    double bs = elp_prj.smn_axis * elp_prj.smn_axis;

    elp_prj.coeff[0] = 0.5*(as + bs + cos(2*elp_prj.rot_angle)*(bs - as));
    elp_prj.coeff[1] = sin(2*elp_prj.rot_angle)*(bs - as);
    elp_prj.coeff[2] = 0.5*(as + bs - cos(2*elp_prj.rot_angle)*(bs - as));
    elp_prj.coeff[3] = -2*elp_prj.center.x()*elp_prj.coeff[0] - elp_prj.center.y()*elp_prj.coeff[1];
    elp_prj.coeff[4] = -2*elp_prj.center.y()*elp_prj.coeff[2] - elp_prj.center.x()*elp_prj.coeff[1];
    elp_prj.coeff[5] = elp_prj.center.x()*elp_prj.center.x()*elp_prj.coeff[0] +
                       elp_prj.center.x()*elp_prj.center.y()*elp_prj.coeff[1] +
                       elp_prj.center.y()*elp_prj.center.y()*elp_prj.coeff[2] - as*bs;
}
