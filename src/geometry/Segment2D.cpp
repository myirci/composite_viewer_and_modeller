#include "Segment2D.hpp"

Segment2D::Segment2D() : pt1(osg::Vec2d(0,0)), pt2(osg::Vec2d(0,0)) { }

Segment2D::Segment2D(const osg::Vec2d& p1, const osg::Vec2d& p2) : pt1(p1), pt2(p2) { }

double Segment2D::length() const {
    return (pt1 - pt2).length();
}

double Segment2D::half_length() const {
    return ((pt1 - pt2).length()) / 2.0;
}

Segment2D& Segment2D::operator=(const Segment2D& rhs) {

    pt1 = rhs.pt1;
    pt2 = rhs.pt2;
}

osg::Vec2d Segment2D::direction() const {
    osg::Vec2d diff = pt2 - pt1;
    diff.normalize();
    return diff;
}

osg::Vec2d Segment2D::mid_point() const {
    return (pt1 + pt2) / 2.0;
}

void Segment2D::translate(const osg::Vec2d& tvec) {
    pt1 += tvec;
    pt2 += tvec;
}

/*
 * Rotates the segment around its mid point in counter clockwise direction by an angle
 * given in radian.
*/
void Segment2D::rotate(double angle) {

    osg::Vec2d midpt = mid_point();

    translate(-midpt);

    double c = cos(angle);
    double s = sin(angle);

    double x = c * pt1.x() - s * pt1.y();
    double y = s * pt1.x() + c * pt1.y();
    pt1.x() = x;
    pt1.y() = y;

    x = c * pt2.x() - s * pt2.y();
    y = s * pt2.x() + c * pt2.y();
    pt2.x() = x;
    pt2.y() = y;

    translate(midpt);
}


