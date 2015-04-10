#ifndef LINE2D_HPP
#define LINE2D_HPP

#include <osg/Vec2d>
#include <osg/Vec3d>

class Line2D {
public:

    // construct from coefficients
    Line2D(double a, double b, double c);
    Line2D(const osg::Vec3d& coeff);

    // from two points: l = p1xp2, p1 & p2 are in homemogenous coordinates
    Line2D(const osg::Vec3d& pt1, const osg::Vec3d& pt2);
    Line2D(const osg::Vec2d& pt1, const osg::Vec2d& pt2);

    // creates the line which is perpendicular to the given line at the given point pt
    Line2D(const Line2D& other, const osg::Vec2d& pt);

    // copy constructor
    Line2D(const Line2D& other);

    bool intersect(const Line2D& other, osg::Vec2d& pt) const;
    bool is_parallel(const Line2D& other) const;
    double distance(const osg::Vec2d& pt) const;
    double evaluate(const osg::Vec2d& pt) const;
    void project(const osg::Vec2d& pt, osg::Vec2d& proj) const;
    double get_x_at_y(double y) const;
    double get_y_at_x(double x) const;

    void get_tangent_vec1(osg::Vec2d& vec1) const;
    void get_tangent_vec2(osg::Vec2d& vec2) const;
    void get_normal_vec1(osg::Vec2d& vec1) const;
    void get_normal_vec2(osg::Vec2d& vec2) const;

    const osg::Vec3d& get_coefficients() const { return line; }
    osg::Vec3d& get_coefficients() { return line; }
    virtual void print() const;

protected:
    osg::Vec3d line;
};

#endif // LINE2D_HPP
