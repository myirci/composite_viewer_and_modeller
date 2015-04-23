#ifndef ELLIPSE_2D_LIGHT_HPP
#define ELLIPSE_2D_LIGHT_HPP

#include "Primitives.hpp"
#include <vector>
#include <osg/Array>

struct Ellipse2DLight {
public:
    Ellipse2DLight(double smj = 0.0, double smn = 0.0, double rot = 0.0, const osg::Vec2d& center_pt = osg::Vec2d(0,0));
    virtual ~Ellipse2DLight() { }
    Ellipse2DLight& operator=(const Ellipse2DLight& rhs);
    void generate_points_on_the_ellipse(osg::Vec2dArray* data, int num) const;
    void generate_points_on_the_ellipse(osg::Vec2dArray* data, int start, int num) const;
    osg::Vec2d center;
    double rot_angle;
    double smn_axis;
    double smj_axis;
};

std::ostream& operator<<(std::ostream& out, const Ellipse2DLight& ellipse);

#endif // ELLIPSE_2D_LIGHT_HPP
