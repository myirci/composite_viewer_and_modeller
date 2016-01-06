#ifndef SEGMENT2D_HPP
#define SEGMENT2D_HPP

#include <osg/Vec2d>

class Segment2D {
public:

    osg::Vec2d pt1;
    osg::Vec2d pt2;

    Segment2D();
    Segment2D(const osg::Vec2d& p1, const osg::Vec2d& p2);
    Segment2D& operator=(const Segment2D& rhs);

    osg::Vec2d direction() const;
    osg::Vec2d mid_point() const;
    double length() const;
    double half_length() const;

    void translate(const osg::Vec2d& tvec);
    void rotate(double angle);


};

#endif // SEGMENT2D_HPP
