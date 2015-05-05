#ifndef OSGUTILIY_HPP
#define OSGUTILIY_HPP

#include <iostream>
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
void print_camera_projection_matrix(const osg::Camera* const cam);
void print_camera_modelview_matrix(const osg::Camera* const cam);
void print_camera_viewport_mapping_matrix(const osg::Camera* constcam);
void print_camera_manipulator_matrix(osgGA::CameraManipulator* manipulator);
void print_camera_frustrum(const osg::Camera* const cam);
void print_camera_orientation(const osg::Camera* const cam);
void print_camera_calibration_matrix(const osg::Camera * const cam);
void print_3x4_camera_projection_matrix(const osg::Camera * const cam);

template<typename Matrix>
void print_matrix(const Matrix& mat, int rows, int cols) {

    for(int i = 0; i < rows; ++i) {
        for(int j = 0; j < cols; ++j)
            std::cout << mat(i, j) << " ";
        std::cout << std::endl;
    }
}

template<typename Matrix>
void transpose(Matrix& mat, int size) {

    Matrix tr_mat;
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            tr_mat(i,j) = mat(j,i);
    mat = tr_mat;
}

template<typename Matrix>
void transpose(const Matrix& mat1, Matrix& mat2, int size) {

    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            mat2(i,j) = mat1(j,i);
}

template<typename Matrix1, typename Matrix2>
void transpose(const Matrix1& mat1, int rows, int cols, Matrix2& mat2) {

    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            mat2(i,j) = mat1(j,i);
}

template<typename Matrix1, typename Matrix2, typename Matrix3>
void multiply(const Matrix1& mat1, int r1, int c1, const Matrix2& mat2, int r2, int c2, Matrix3& mul) {

    if(c1 != r2) {
        std::cerr << "matrix size mismatch!" << std::endl;
        return;
    }

    for(int i = 0; i < r1; ++i) {
        for(int j = 0; j < c2; ++j) {
            mul(i,j) = 0;
            for(int k = 0; k < r2; ++k)
                mul(i,j) += mat1(i,k) * mat2(k,j);
        }
    }
}

#endif // OSGUTILIY_HPP
