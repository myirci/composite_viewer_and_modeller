#ifndef PERS_PROJ_PARAM_HPP
#define PERS_PROJ_PARAM_HPP

#include "../utility/Utility.hpp"
#include <osg/Array>
#include <ostream>

struct PersProjParam {

    double fovy;
    double aspect;
    int width;
    int height;
    double far;
    double near;

    PersProjParam(double fy, int w, int h, double n, double f) : fovy(fy), width(w), height(h), near(n), far(f) {
        aspect = static_cast<double>(width) / static_cast<double>(height);
        c1 = near * tan(deg2rad(fovy/2.0));
    }

    void convert_from_image_coordinates_to_opengl_screen_coordinates(const double x_img, const double y_img, double& x_gls, double& y_gls) {
        x_gls = x_img;
        y_gls = -y_img + height - 1;
    }

    void convert_from_opengl_screen_coordinates_to_normalized_device_coordinates(const double x_gls, const double y_gls, double& x_ndc, double& y_ndc) {
        x_ndc = (2.0/static_cast<double>(width))*(x_gls) - 1;
        y_ndc = (2.0/static_cast<double>(height))*(y_gls) - 1;
    }

    void convert_from_normalized_device_coordinates_to_projected_point_coordinates(const double x_ndc, const double y_ndc, double& x_prj, double& y_prj) {
        x_prj = x_ndc * aspect * c1;
        y_prj = y_ndc * c1;
    }

    void convert_from_opengl_screen_coordinates_to_projected_point_coordinates(const double x_gls, const double y_gls, double& x_prj, double& y_prj) {
        double x_ndc = (2.0/static_cast<double>(width))*(x_gls) - 1;
        double y_ndc = (2.0/static_cast<double>(height))*(y_gls) - 1;
        x_prj = x_ndc * aspect * c1;
        y_prj = y_ndc * c1;
    }

    void convert_from_opengl_screen_coordinates_to_projected_point_coordinates(osg::Vec2dArray* vertices) {
        for(auto it = vertices->begin(); it != vertices->end(); ++it) {
            it->x() = ((2.0/static_cast<double>(width))*(it->x()) - 1) * aspect * c1;
            it->y() = ((2.0/static_cast<double>(height))*(it->y()) - 1) * c1;
        }
    }

    void construct_projection_matrix(osg::Matrixd& proj_mat);

    void construct_viewport_mapping_matrix(osg::Matrixd& vp_mat);

private:
    double c1;
};

std::ostream& operator<<(std::ostream& out, const PersProjParam& ppp);

#endif // PERS_PROJ_PARAM_HPP
