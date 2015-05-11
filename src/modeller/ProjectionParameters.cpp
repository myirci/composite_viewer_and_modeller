#include "ProjectionParameters.hpp"
#include "../geometry/Ellipse2D.hpp"

ProjectionParameters::ProjectionParameters(double fy, int w, int h, double n, double f) : fovy(fy), width(w), height(h), near(n), far(f) {

    aspect = static_cast<double>(width) / static_cast<double>(height);
    c1 = (near * tan(deg2rad(fovy/2.0)));
}

void ProjectionParameters::convert_from_image_coordinates_to_logical_device_coordinates(const osg::Vec2d& img_coord, osg::Vec2d& log_coord) {

    log_coord.x() = img_coord.x();
    log_coord.y() = -img_coord.y() + height - 1;
}

void ProjectionParameters::convert_from_logical_device_coordinates_to_image_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& img_coord) {

    img_coord.x() = log_coord.x();
    img_coord.y() = -log_coord.y() + height - 1;
}

void ProjectionParameters::convert_from_viewport_coordinates_to_logical_device_coordinates(const osg::Vec2d& vp_coord, osg::Vec2d& log_coord) {

    log_coord.x() = vp_coord.x() + static_cast<double>(width) / 2.0;
    log_coord.y() = vp_coord.y() + static_cast<double>(height) / 2.0;
}

void ProjectionParameters::convert_from_logical_device_coordinates_to_viewport_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& vp_coord) {

    vp_coord.x() = log_coord.x() - static_cast<double>(width) / 2.0;
    vp_coord.y() = log_coord.y() - static_cast<double>(height) / 2.0;
}

void ProjectionParameters::convert_from_normalized_device_coordinates_to_viewport_coordinates(const osg::Vec2d& ndc_coord, osg::Vec2d& vp_coord) {

    vp_coord.x() = (static_cast<double>(width) / 2.0)*ndc_coord.x();
    vp_coord.y() = (static_cast<double>(height) / 2.0)*ndc_coord.y();
}

void ProjectionParameters::convert_from_viewport_coordinates_to_normalized_device_coordinates(const osg::Vec2d& vp_coord, osg::Vec2d& ndc_coord) {

    ndc_coord.x() = (2.0 / static_cast<double>(width))*vp_coord.x();
    ndc_coord.y() = (2.0 / static_cast<double>(height))*vp_coord.y();
}

void ProjectionParameters::convert_from_projected_coordinates_to_normalized_device_coordinates(const osg::Vec2d& prj_coord, osg::Vec2d& ndc_coord) {

    ndc_coord.x() = prj_coord.x() / (aspect * c1);
    ndc_coord.y() = prj_coord.y() / c1;
}

void ProjectionParameters::convert_from_normalized_device_coordinates_to_projected_coordinates(const osg::Vec2d& ndc_coord, osg::Vec2d& prj_coord) {

    prj_coord.x() = ndc_coord.x() * aspect * c1;
    prj_coord.y() = ndc_coord.y() * c1;
}

void ProjectionParameters::convert_from_logical_device_coordinates_to_projected_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& prj_coord) {

    // logical device coord >> viewport coordinates
    osg::Vec2d vp_coord, ndc_coord;
    convert_from_logical_device_coordinates_to_viewport_coordinates(log_coord, vp_coord);
    // viewport coord >> normalized device coordinates
    convert_from_viewport_coordinates_to_normalized_device_coordinates(vp_coord, ndc_coord);
    // normalized device coord >> projected coordinates
    convert_from_normalized_device_coordinates_to_projected_coordinates(ndc_coord, prj_coord);
}

void ProjectionParameters::convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(const Ellipse2D& elp_dev, Ellipse2D& elp_prj,  bool with_coeffs) {

    // rot angle (angle between the major axis and the positive x-axis) does not change
    elp_prj.rot_angle = elp_dev.rot_angle;

    for(size_t i = 0; i < 4; ++i)
        convert_from_logical_device_coordinates_to_projected_coordinates(elp_dev.points[i], elp_prj.points[i]);

    // calculate the remaining parameters and the coefficients
    elp_prj.center.x() = (elp_prj.points[0].x() + elp_prj.points[1].x())/2.0;
    elp_prj.center.y() = (elp_prj.points[0].y() + elp_prj.points[1].y())/2.0;
    elp_prj.smj_axis = (elp_prj.points[0] - elp_prj.center).length();
    elp_prj.smn_axis = (elp_prj.points[2] - elp_prj.center).length();

    if(with_coeffs)
        elp_prj.calculate_coefficients_from_parameters();
}

void ProjectionParameters::construct_perpective_projection_matrix(osg::Matrixd& proj_mat) {

    proj_mat = osg::Matrixd::identity();
    double half_fov = deg2rad(0.5*fovy);
    proj_mat(0,0) = cos(half_fov)/(sin(half_fov)*aspect);
    proj_mat(1,1) = cos(half_fov)/sin(half_fov);
    proj_mat(2,2) = -(far + near) / (far - near);
    proj_mat(2,3) = (-2*far*near) / (far - near);
    proj_mat(3,2) = -1;
    proj_mat(3,3) = 0;
}

void ProjectionParameters::construct_viewport_mapping_matrix(osg::Matrixd& vp_mat) {

    double fVal = 1.0;
    double nVal = 0.0;
    vp_mat = osg::Matrixd::identity();
    vp_mat(0,0) = width/2.0;
    vp_mat(0,3) = width/2.0;
    vp_mat(1,1) = height/2.0;
    vp_mat(1,3) = height/2.0;
    vp_mat(2,2) = (fVal - nVal) / 2.0;
    vp_mat(2,3) = (fVal + nVal) / 2.0;
}

std::ostream& operator<<(std::ostream& out, const ProjectionParameters& pp) {

    out << "far   : " << pp.far << std::endl;
    out << "near  : " << pp.near << std::endl;
    out << "fovy  : " << pp.fovy << std::endl;
    out << "aspect: " << pp.aspect << std::endl;
    out << "width : " << pp.width << std::endl;
    out << "height: " << pp.height << std::endl;
    return out;
}
