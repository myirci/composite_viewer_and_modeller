#ifndef RAY3D_HPP
#define RAY3D_HPP

#include <osg/Vec3d>

class Ray3D {
public:

    Ray3D() : start(osg::Vec3d(0,0,0)), direction(osg::Vec3d(0,0,0)) { }
    Ray3D(const osg::Vec3d& s, const osg::Vec3d& d);
    Ray3D(const osg::Vec3d& d);
    void translate(const osg::Vec3d& tvec);
    void move_to(const osg::Vec3d& pt);
    const osg::Vec3d& start_point() const       { return start; }
    osg::Vec3d& start_point()                   { return start; }
    const osg::Vec3d& direction_vector() const  { return direction; }
    osg::Vec3d& direction_vector()              { return direction; }

private:
    osg::Vec3d start;
    osg::Vec3d direction;
};

#endif // RAY3D_HPP
