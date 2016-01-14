#ifndef PLANE3D_HPP
#define PLANE3D_HPP

#include <osg/Vec4d>
#include <osg/Uniform>
#include <Eigen/Dense>

class Ray3D;

class Plane3D {
public:

    Plane3D() { }
    Plane3D(const osg::Vec4d& vec);
    Plane3D(double a, double b, double c, double d);
    Plane3D(const osg::Vec3d& normal, const osg::Vec3d& p);
    Plane3D(const Plane3D& pl);
    Plane3D(const Eigen::Vector3d& normal, const Eigen::Vector3d& pt); 
    Plane3D(const osg::Vec3d& p1, const osg::Vec3d& p2, const osg::Vec3d& p3); // points has to be non - collinear

    bool is_on_plane(const osg::Vec3d& pt) const;
    void get_normal(osg::Vec3d& normal) const;
    void get_unit_normal(osg::Vec3d& normal) const;
    void print() const;
    void rotate(const osg::Matrix3d& mat);
    void translate(const osg::Vec3d& tvec);
    void transform(const osg::Matrixd& mat);
    void hessian_normal_form();
    const osg::Vec4d& get_plane() const { return plane; }
    osg::Vec4d& get_plane()             { return plane; }

private:
    osg::Vec4d plane;
};

#endif // Plane3D_HPP
