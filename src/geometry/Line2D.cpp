#include "Line2D.hpp"
#include <iostream>


Line2D::Line2D(double a, double b, double c) : line(osg::Vec3d(a, b, c)) { }

Line2D::Line2D(const osg::Vec3d& coeff) : line (coeff) { }

Line2D::Line2D(const osg::Vec3d& pt1, const osg::Vec3d& pt2) {
    line = pt1 ^ pt2;
}

Line2D::Line2D(const osg::Vec2d& pt1, const osg::Vec2d& pt2) {

    osg::Vec3d pt1_h(pt1.x(), pt1.y(), 1);
    osg::Vec3d pt2_h(pt2.x(), pt2.y(), 1);
    line = pt1_h ^ pt2_h;
}

Line2D::Line2D(const Line2D& other, const osg::Vec2d& pt) {

    line.x() = other.line.y();
    line.y() = -other.line.x();
    line.z() =  other.line.x() * pt.y() - other.line.y() * pt.x();
}

Line2D::Line2D(const Line2D& other) {
    line = other.line;
}

void Line2D::get_normal_vec1(osg::Vec2d& vec1) const {

    vec1.x() = line.x();
    vec1.y() = line.y();
}

void Line2D::get_normal_vec2(osg::Vec2d& vec2) const {

    vec2.x() = -line.x();
    vec2.y() = -line.y();
}

void Line2D::get_tangent_vec1(osg::Vec2d& vec1) const {

    vec1.x() = line.y();
    vec1.y() = -line.x();
}

void Line2D::get_tangent_vec2(osg::Vec2d& vec2) const {

    vec2.x() = -line.y();
    vec2.y() = line.x();
}

bool Line2D::is_parallel(const Line2D& other) const {
    return (line.x() / other.line.x()) == (line.y() / other.line.y());
}

bool Line2D::intersect(const Line2D& other, osg::Vec2d& pt) const {

    if(is_parallel(other)) return false;
    osg::Vec3d intersecting_pt = line ^ other.line;
    pt.x() = intersecting_pt.x() / intersecting_pt.z();
    pt.y() = intersecting_pt.y() / intersecting_pt.z();
}

void Line2D::print() const {

    std::cout << "coeff: " ;
    std::cout << line.x() << " " << line.y() << " " << line.z() << std::endl;
    std::cout << "tangent vectors: (" << line.y() << ", " << -line.x() << ") and (" << -line.y() << ", " << line.x() << ")" << std::endl;
    std::cout << "normal vectors: (" << line.x() << ", " << line.y() << ") and (" << -line.x() << ", " << -line.y() << ")" << std::endl;
}

double Line2D::get_x_at_y(double y) const {
    return (-line.y()/line.x())*y - line.z() / line.x();
}

double Line2D::get_y_at_x(double x) const {
    return (-line.x()/line.y())*x - line.z() / line.y();
}

double Line2D::evaluate(const osg::Vec2d& pt) const {
    return pt.x()*line.x() + pt.y()*line.y() + line.z();
}

double Line2D::distance(const osg::Vec2d& pt) const {
    return std::abs(line.x()*pt.x() + line.y()*pt.y() + line.z()) / std::sqrt(line.x()*line.x() + line.y()*line.y());
}

void Line2D::project(const osg::Vec2d& pt, osg::Vec2d& proj) const {

    double var = evaluate(pt)/(line.x()*line.x() + line.y()*line.y());
    proj.x() = pt.x() - line.x()*var;
    proj.y() = pt.y() - line.y()*var;

}


