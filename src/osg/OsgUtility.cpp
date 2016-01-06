#include "OsgUtility.hpp"
#include "../geometry/Circle3D.hpp"

#include <osg/MatrixTransform>
#include <osgViewer/View>
#include <osgDB/ReadFile>
#include <osg/Depth>
#include <osg/Texture2D>
#include <osg/ShapeDrawable>
#include <osg/LineWidth>

#include <wx/gdicmn.h>

#include <Eigen/Dense>

osg::Camera* create_background_camera(int left, int right, int bottom, int top) {

    osg::Camera* camera = new osg::Camera;
    camera->setClearMask(0);
    camera->setCullingActive(false);
    camera->setAllowEventFocus(false);
    camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
    camera->setRenderOrder(osg::Camera::POST_RENDER);
    camera->setProjectionMatrix(osg::Matrix::ortho2D(left, right, bottom, top));
    camera->setViewMatrix(osg::Matrix::identity());
    camera->setViewport(0, 0, right, top);
    osg::StateSet* ss = camera->getOrCreateStateSet();
    ss->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
    ss->setAttributeAndModes(new osg::Depth(osg::Depth::LEQUAL, 1.0, 1.0));
    return camera;
}

osg::Geode* create_textured_quad( osg::Image* image, wxSize& size) {

    size.x = image->s();
    size.y = image->t();

    // Create the geometry
    osg::Vec3 pos_vec = osg::Vec3(0.0f, 0.0f, 0.0f);
    osg::Vec3 width_vec(image->s(), 0.0f, 0.0);
    osg::Vec3 height_vec(0.0, image->t(), 0.0);
    osg::Geometry* quad = osg::createTexturedQuadGeometry(pos_vec, width_vec, height_vec);

    // Add texture to the geometry
    osg::Texture2D* texture = new osg::Texture2D();
    texture->setResizeNonPowerOfTwoHint(false);
    texture->setImage(image);

    osg::Geode* tex_geode = new osg::Geode;
    tex_geode->addDrawable(quad);
    osg::StateSet* stateset = tex_geode->getOrCreateStateSet();
    stateset->setTextureAttributeAndModes(0, texture, osg::StateAttribute::ON);
    return tex_geode;
}

osg::Geometry* create_3D_circle(const Circle3D& circle, int approx) {

    osg::Geometry* circleGeom = new osg::Geometry();
    osg::ref_ptr<osg::Vec3Array> v = new osg::Vec3Array;
    circle.generate_data(v, approx);
    circleGeom->setVertexArray(v.get());
    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(0.0f,1.0f,0.0f,1.0f));
    circleGeom->setColorArray(colors, osg::Array::BIND_OVERALL);
    // set the normal in the same way color.
    osg::Vec3Array* normals = new osg::Vec3Array;
    osg::Vec3d normal(circle.normal[0], circle.normal[1], circle.normal[2]);
    normals->push_back(normal);
    circleGeom->setNormalArray(normals, osg::Array::BIND_OVERALL);
    circleGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, approx));
    return circleGeom;
}

osg::MatrixTransform* display_vector3d(const osg::Vec3d& pt, const osg::Vec3d& vec, const osg::Vec4d& color) {

    // arrow along wframe z-axis
    double size = vec.length();
    float cylinder_radius = size / 20.0;
    float cone_radius     = size / 15.0;
    float cone_height     = size / 5.0;

    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
    osg::ref_ptr<osg::ShapeDrawable> cyl = new osg::ShapeDrawable();
    cyl->setShape(new osg::Cylinder(osg::Vec3d(0.0, 0.0, size/2.0), cylinder_radius, size));
    cyl->setColor(color);
    geode->addDrawable(cyl.get());
    osg::ref_ptr<osg::ShapeDrawable> cone = new osg::ShapeDrawable;
    cone->setShape(new osg::Cone(osg::Vec3d(0.0, 0.0, size), cone_radius, cone_height));
    cone->setColor(color);
    geode->addDrawable(cone.get());

    osg::MatrixTransform* mtrans = new osg::MatrixTransform();
    mtrans->setMatrix(osg::Matrix::rotate(acos(vec.z() / size), osg::Vec3d(0,0,1) ^ vec) * osg::Matrix::translate(pt));
    mtrans->addChild(geode.get());
    return mtrans;
}

osg::Geode* display_lines(osg::Vec3Array* vertices, const osg::Vec4d& color) {

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(color);
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setVertexArray(vertices);
    geom->setColorArray(colors.get());
    geom->setColorBinding(osg::Geometry::BIND_OVERALL);
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertices->size()));
    geode->addDrawable(geom.release());
    osg::LineWidth* linewidth = new osg::LineWidth();
    linewidth->setWidth(0.3);
    geode->getOrCreateStateSet()->setAttributeAndModes(linewidth, osg::StateAttribute::ON);
    geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
    return geode.release();
}

osg::Geode* display_point(const osg::Vec3d& pt, const osg::Vec4d& color) {

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<osg::ShapeDrawable> sphere = new osg::ShapeDrawable;
    sphere->setShape(new osg::Sphere(pt, 0.5));
    sphere->setColor(color);
    geode->addDrawable(sphere.release());
    return geode.release();
}


// given two 3D circles, this function calculates the
// transformation matrix that transforms circle1 to circle2.
void calculate_transformation_matrix(const Circle3D& circle1, const Circle3D& circle2, osg::Matrixd& mat) {

    // 1) translate to the origin
    osg::Matrixd mat_org = osg::Matrixd::identity();
    mat_org(0,3) = -circle1.center[0];
    mat_org(1,3) = -circle1.center[1];
    mat_org(2,3) = -circle1.center[2];

    // 2) scaling
    double scale = circle2.radius / circle1.radius;
    osg::Matrixd mat_scale = osg::Matrixd::identity();
    mat_scale(0,0) = scale;
    mat_scale(1,1) = scale;
    mat_scale(2,2) = scale;

    // 3) rotation
    Eigen::Vector3d rot_axis = circle1.normal.cross(circle2.normal);
    osg::Matrixd mat_rotate = osg::Matrixd::rotate(acos(circle1.normal.dot(circle2.normal)), rot_axis[0], rot_axis[1], rot_axis[2]);
    transpose(mat_rotate, 4); // since osg uses row vector format

    // 4) translate to the position of the circle-2
    osg::Matrixd mat_trans = osg::Matrixd::identity();
    mat_trans(0,3) = circle2.center[0];
    mat_trans(1,3) = circle2.center[1];
    mat_trans(2,3) = circle2.center[2];

    // 5) combined transformation matrix
    mat = mat_trans * mat_rotate * mat_scale * mat_org;
}

void calculate_transformation_matrix_without_scale(const Circle3D& circle1, const Circle3D& circle2, osg::Matrixd& mat) {

    // 1) translation matrix to translate circle-1 to the origin
    osg::Matrixd mat_org = osg::Matrixd::identity();
    mat_org(0,3) = -circle1.center[0];
    mat_org(1,3) = -circle1.center[1];
    mat_org(2,3) = -circle1.center[2];

    // 2) rotation matrix to align circle-1 normal to circle-2 normal
    Eigen::Vector3d rot_axis = circle1.normal.cross(circle2.normal);
    osg::Matrixd mat_rotate = osg::Matrixd::rotate(acos(circle1.normal.dot(circle2.normal)), rot_axis[0], rot_axis[1], rot_axis[2]);
    transpose(mat_rotate, 4);

    // 3) translation matrix to transalte circle-1 to the center of circle-2
    osg::Matrixd mat_trans = osg::Matrixd::identity();
    mat_trans(0,3) = circle2.center[0];
    mat_trans(1,3) = circle2.center[1];
    mat_trans(2,3) = circle2.center[2];

    // 4) combined transformation matrix
    mat = mat_trans * mat_rotate * mat_org;
}

double squared_distance(const osg::Vec3d& pt1, const osg::Vec3d& pt2) {
    return (pt2 - pt1).length2();
}

void print_3x4_camera_projection_matrix(const osg::Camera * const cam) {

    std::cout << "3x4 camera  projection matrix when column vector format is used:" << std::endl;
    osg::Matrixd mat = cam->getProjectionMatrix();
    transpose(mat, 4);
    osg::Matrix3x4d mat_proj;
    for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 4; ++j)
            mat_proj(i,j) = mat(i,j);

    for(size_t j = 0; j < 4; ++j)
        mat_proj(2,j) = mat(3,j);
    print_matrix(mat_proj, 3, 4);
}

void print_camera_projection_matrix(const osg::Camera* const cam) {

    std::cout << "Projection matrix when column vector format is used:" << std::endl;
    osg::Matrixd mat_transposed;
    transpose(cam->getProjectionMatrix(), mat_transposed, 4);
    print_matrix(mat_transposed, 4, 4);
}

void print_camera_modelview_matrix(const osg::Camera* const cam) {

    std::cout << "Model-view matrix when column vector format is used:" << std::endl;
    osg::Matrixd mat_transposed;
    transpose(cam->getViewMatrix(), mat_transposed, 4);
    print_matrix(mat_transposed, 4, 4);
}

void print_camera_viewport_mapping_matrix(const osg::Camera * const cam) {

    std::cout << "viewport mapping matrix when column vector format is used:" << std::endl;
    osg::Matrixd mat_transposed;
    transpose(cam->getViewport()->computeWindowMatrix(), mat_transposed, 4);
    print_matrix(mat_transposed, 4, 4);
}

void print_camera_frustrum(const osg::Camera* const cam) {

    osg::Matrixd mat = cam->getProjectionMatrix();
    double left, right, bottom, top, near, far;
    mat.getFrustum(left, right, bottom, top, near, far);
    std::cout << "Frustrum: " << std::endl;
    std::cout << "left:     " << left   <<  std::endl;
    std::cout << "right:    " << right  <<  std::endl;
    std::cout << "bottom:   " << bottom <<  std::endl;
    std::cout << "top:      " << top    <<  std::endl;
    std::cout << "near:     " << near   <<  std::endl;
    std::cout << "far:      " << far    <<  std::endl;

    double fovy, aspect;
    if(mat.getPerspective(fovy, aspect, near, far)) {
        std::cout << "fovy:     " << fovy   <<  std::endl;
        std::cout << "aspect:   " << aspect <<  std::endl;
    }
}

void print_camera_orientation(const osg::Camera * const cam) {

    std::cout << "Camera orientation: " << std::endl;
    osg::Vec3d eye, center, up;
    cam->getViewMatrixAsLookAt(eye, center, up);
    std::cout << "eye: " << eye.x() << " " << eye.y() << " " << eye.z() << std::endl;
    std::cout << "center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;
    std::cout << "up: " << up.x() << " " << up.y() << " " << up.z() << std::endl;
}

void print_camera_manipulator_matrix(osgGA::CameraManipulator* manipulator) {

    std::cout << "Camera manipulator matrix: " << std::endl;
    osg::Matrixd mat_transposed;
    transpose(manipulator->getMatrix(), mat_transposed, 4);
    print_matrix(mat_transposed, 4, 4);
}

