#include "Ellipse2DLight.hpp"

Ellipse2DLight::Ellipse2DLight(double smj, double smn, double rot, const osg::Vec2d& center_pt) :
    smj_axis(smj), smn_axis(smn), rot_angle(rot), center(center_pt) { }

Ellipse2DLight& Ellipse2DLight::operator=(const Ellipse2DLight& rhs) {

    smj_axis = rhs.smj_axis;
    smn_axis = rhs.smn_axis;
    center = rhs.center;
    rot_angle = rhs.rot_angle;
}

void Ellipse2DLight::generate_points_on_the_ellipse(osg::Vec2dArray* data, int num) const {

    double step = TWO_PI/num;
    osg::Vec2d e1(cos(rot_angle), sin(rot_angle));
    osg::Vec2d e2(-sin(rot_angle), cos(rot_angle));
    for(double d = 0.0; d < TWO_PI; d += step)
        data->push_back(center + e1*smj_axis*cos(d) + e2*smn_axis*sin(d));
}

void Ellipse2DLight::generate_points_on_the_ellipse(osg::Vec2dArray* data, int start, int num) const {

    double step = TWO_PI/num;
    osg::Vec2d e1(cos(rot_angle), sin(rot_angle));
    osg::Vec2d e2(-sin(rot_angle), cos(rot_angle));
    for(double d = 0.0; d < TWO_PI; d += step)
        (*data)[start++] = center + e1*smj_axis*cos(d) + e2*smn_axis*sin(d);

}

std::ostream& operator<<(std::ostream& out, const Ellipse2DLight& ellipse) {

    out << "center  : " << ellipse.center.x() << " " << ellipse.center.y() << std::endl;
    out << "smj-axis: " << ellipse.smj_axis << std::endl;
    out << "smn-axis: " << ellipse.smn_axis << std::endl;
    out << "rotation: " << ellipse.rot_angle << std::endl;
    return out;
}
