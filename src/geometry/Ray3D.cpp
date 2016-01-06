#include "Ray3D.hpp"

Ray3D::Ray3D(const osg::Vec3d& s, const osg::Vec3d& d) : start(s), direction(d) { }

Ray3D::Ray3D(const osg::Vec3d& d) : start(osg::Vec3d(0,0,0)), direction(d) { }

void Ray3D::translate(const osg::Vec3d& tvec) {
    start = start + tvec;
}

void Ray3D::move_to(const osg::Vec3d& pt) {
    start = pt;
}








