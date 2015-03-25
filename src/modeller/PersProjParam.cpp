#include "PersProjParam.hpp"

void PersProjParam::construct_projection_matrix(osg::Matrixd& proj_mat) {
    proj_mat = osg::Matrixd::identity();
    double half_fov = deg2rad(0.5*fovy);
    proj_mat(0,0) = cos(half_fov)/(sin(half_fov)*aspect);
    proj_mat(1,1) = cos(half_fov)/sin(half_fov);
    proj_mat(2,2) = -(far + near) / (far - near);
    proj_mat(2,3) = (-2*far*near) / (far - near);
    proj_mat(3,2) = -1;
    proj_mat(3,3) = 0;
}

void PersProjParam::construct_viewport_mapping_matrix(osg::Matrixd& vp_mat) {
    vp_mat = osg::Matrixd::identity();
    vp_mat(0,0) = width/2.0;
    vp_mat(0,3) = width/2.0;
    vp_mat(1,1) = height/2.0;
    vp_mat(1,3) = height/2.0;
}

std::ostream& operator<<(std::ostream& out, const PersProjParam& ppp) {
    out << "far   : " << ppp.far << std::endl;
    out << "near  : " << ppp.near << std::endl;
    out << "fovy  : " << ppp.fovy << std::endl;
    out << "aspect: " << ppp.aspect << std::endl;
    out << "width : " << ppp.width << std::endl;
    out << "height: " << ppp.height << std::endl;
    return out;
}
