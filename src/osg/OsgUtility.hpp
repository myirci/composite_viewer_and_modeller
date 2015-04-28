#ifndef OSGUTILIY_HPP
#define OSGUTILIY_HPP

#include <osg/Camera>
#include <osgGA/CameraManipulator>

class wxSize;
class Circle3D;

osg::Camera* create_background_camera(int left, int right, int bottom, int top);
osg::Geode* create_textured_quad(osg::Image* image, wxSize& size);
osg::Geometry* create_3D_circle(const Circle3D& circle, int approx);
osg::MatrixTransform* display_vector3d(const osg::Vec3d& pt, const osg::Vec3d& vec, const osg::Vec4d& color);
double squared_distance(const osg::Vec3d& pt1, const osg::Vec3d& pt2);
void calculate_transformation_matrix(const Circle3D& circle1, const Circle3D& circle2, osg::Matrixd& mat);
void transpose(const osg::Matrixd& mat, osg::Matrixd& mat_transposed);
void transpose(osg::Matrixd& mat);
void print_projection_matrix(osg::Camera* cam);
void print_modelview_matrix(osg::Camera* cam);
void print_window_matrix(osg::Camera* cam);
void print_camera_manipulator_matrix(osgGA::CameraManipulator* manipulator);
void print_frustrum(osg::Camera* cam);
void print_transposed_osg_matrix(const osg::Matrixd& mat);
void print_osg_matrix(const osg::Matrixd& mat, std::ostream& out);

#endif // OSGUTILIY_HPP
