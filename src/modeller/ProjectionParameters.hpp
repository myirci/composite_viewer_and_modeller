#ifndef PROJECTION_PARAMETERS_HPP
#define PROJECTION_PARAMETERS_HPP

#include "../utility/Utility.hpp"
#include <osg/Geode>
#include <osg/Array>
#include <ostream>

class Ellipse2D;
class Segment2D;

struct ProjectionParameters {

    double fovy;        // vertical field of view
    double aspect;      // aspect ratio vp_width / vp_height
    int width;          // width of the display window = view_port width (for this application)
    int height;         // heigth of the display window = view_port height (for this application)
    double far;         // far clipping plane z = -far, far > 0
    double near;        // near clipping plane z = -near, near > 0

    /* Coordinate systems:
     *
     * 1) Device coordinate system
     *      - Defined w.r.t to the displaying window. Specific to GUI lib.
     *      - For wxWidgets this coordinate system is identical to image coordinate frame.
     *
     * 2) Logical device coordinate system:
     *      - Defined w.r.t to the displaying window. Specific to application, for this application:
     *      - origin >> left bottom corner of the displaying window, x+ >> right, y+ >> up, unit >> pixels, dimension >> 2
     *
     * 3) Image coordinate system
     *      - Defined w.r.t to the displaying window. Image is displayed on entire displaying window.
     *      - origin >> left top corner of the displaying window, x+ >> right, y+ >> down, unit >> pixels, dimension >> 2
     *
     * 4) Viewport coordinate system
     *      - Defined w.r.t to the displaying window.
     *      - origin >> center of the displaying window, x+ >> right, y+ >> up, unit >> pixels, dimension >> 2
     *
     * 5) NDC - normalized device coordinates are defined w.r.t a coordinate frame which is located at the center of a
     *    2x2x2 cube.
     *
     * 6) Projected coordinates - Points in the camera coordinate system are projected onto the near clipping plane. The
     *    coordinates of the projected points in camera coordinate system are called projected coordinates.
    */

    ProjectionParameters(double fy, int w, int h, double n, double f);

    // logical device coordinates <<->> image coordinates
    inline void convert_from_image_coordinates_to_logical_device_coordinates(const osg::Vec2d& img_coord, osg::Vec2d& log_coord);
    inline void convert_from_logical_device_coordinates_to_image_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& img_coord);

    // logical device coordinates <<->> viewport coordinates
    inline void convert_from_viewport_coordinates_to_logical_device_coordinates(const osg::Vec2d& vp_coord, osg::Vec2d& log_coord);
    inline void convert_from_logical_device_coordinates_to_viewport_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& vp_coord);

    // viewport coordinates <<->> normalized device coordinates
    inline void convert_from_viewport_coordinates_to_normalized_device_coordinates(const osg::Vec2d& vp_coord, osg::Vec2d& ndc_coord);
    inline void convert_from_normalized_device_coordinates_to_viewport_coordinates(const osg::Vec2d& ndc_coord, osg::Vec2d& vp_coord);

    // normalized device coordinates <<->> projected coordinates
    inline void convert_from_normalized_device_coordinates_to_projected_coordinates(const osg::Vec2d& ndc_coord, osg::Vec2d& prj_coord);
    inline void convert_from_projected_coordinates_to_normalized_device_coordinates(const osg::Vec2d& prj_coord, osg::Vec2d& ndc_coord);

    // combined transformations
    void convert_from_logical_device_coordinates_to_projected_coordinates(const osg::Vec2d& log_coord, osg::Vec2d& prj_coord);

    // other conversions
    void convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(const Ellipse2D& elp_dev, Ellipse2D& elp_prj, bool with_coeffs = true);
    void convert_segment_from_logical_device_coordinates_to_projected_coordinates(const Segment2D& seg_dev, Segment2D& seg_prj);

    // perspective projection related functions
    void construct_perpective_projection_matrix(osg::Matrixd& proj_mat);
    void construct_viewport_mapping_matrix(osg::Matrixd& vp_mat);

private:
    double c1;

};

std::ostream& operator<<(std::ostream& out, const ProjectionParameters& pp);

#endif // PROJECTION_PARAMETERS_HPP
