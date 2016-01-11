#ifndef IMAGE_MODELLER_HPP
#define IMAGE_MODELLER_HPP

#include <vector>
#include <memory>
#include "../image/algorithms/Algorithms.hpp"
#include "components/GeneralizedCylinder.hpp"
#include <osg/Geometry>
#include <otbImage.h>
#include <wx/gdicmn.h>

class Rectangle2D;
class UIHelper;
class Line2D;
class Segment2D;
class Circle3D;
class Ellipse2D;
class ModelSolver;
class ComponentSolver;
class GeneralizedCylinder;
class ProjectionParameters;
class OsgWxGLCanvas;
class CircleEstimator;

enum class drawing_mode : unsigned char {
    mode_0,     // do nothing
    mode_1,     // semi_major_axis
    mode_2,     // semi_minor_axis
    mode_3      // spine drawing mode
};

enum class component_type : unsigned char {
    generalized_cylinder,
    cuboid,
    sphere
};

enum class spine_constraints : unsigned char {
    none,            // no constraint
    straight_planar, // spine points constitutes a 3D line
    planar           // spine points are planar
};

enum class section_constraints : unsigned char {
    none,
    constant,
    linear_scaling
};

enum class spine_drawing_mode : unsigned char {
    continuous,
    piecewise_linear
};

enum class projection_type : unsigned char {
    perspective,
    orthographic,
    orthogonality_constraint
};

class ImageModeller {
public:
    ImageModeller(const wxString& fpath, const std::shared_ptr<ProjectionParameters>& pp, OsgWxGLCanvas* canvas);
    ~ImageModeller();

    component_type comp_type;           // type of the component being modelled
    spine_constraints sp_constraints;   // constraints on the main axis of the component
    section_constraints sc_constraints; // constraints on the sections of the generalized cylinders and cuboids
    spine_drawing_mode spd_mode;        // either continuous or piecewise

private:

    // otb related data members
    OtbImageType::Pointer m_gimage;                         // for gradient image
    OtbImageType::Pointer m_bimage;                         // for binary image

    // osg related data members
    OsgWxGLCanvas* m_canvas;
    osg::ref_ptr<GeneralizedCylinder> m_gcyl;               // for generalized cylinder modelling
    osg::ref_ptr<osg::Vec2dArray> m_vertices;               // for storing user clicks

    osg::Vec2d m_mouse;                                     // current position of the mouse updated by osgWxGLCanvas
    drawing_mode m_mode;                                    // internal modes of component drawing
    bool m_bimg_exists;                                     // indicates if the binary image exists or not
    bool m_right_click;                                     // set to true in the case of user right click
    bool m_left_click;                                      // set to true in the case of user left click

    rendering_type m_rtype;                                 // rendering type for generalized cylinders

    Circle3D* m_first_circle;
    Circle3D* m_last_circle;
    std::unique_ptr<Ellipse2D> m_first_ellipse;

    std::unique_ptr<Segment2D> m_lsegment; // last segment
    std::unique_ptr<Segment2D> m_dsegment; // dynamic segment

    std::unique_ptr<UIHelper> m_uihelper;
    std::unique_ptr<Rectangle2D> m_rect;
    std::unique_ptr<ModelSolver> m_solver;
    std::unique_ptr<ComponentSolver> m_component_solver;
    std::shared_ptr<ProjectionParameters> m_pp;
    std::unique_ptr<CircleEstimator> m_circle_estimator;

    double m_fixed_depth;
    double m_scale_factor;

    osg::Vec2d m_tvec;

    // for ray cast display
    bool m_display_raycast;
    osg::ref_ptr<osg::Vec2dArray> m_raycast;

public:

    // Public member functions
    void Initialize2DDrawingInterface(osg::Geode* geode);
    void Reset2DDrawingInterface();
    void OnLeftClick(double x, double y);
    void OnRightClick(double x, double y);
    void OnMouseMove(double x, double y);
    void DebugPrint();
    void SaveModel(const std::string& path);
    void DeleteModel();
    void DeleteSelectedComopnents(std::vector<int>& index_vector);
    void SetRenderingType(rendering_type rtype);
    void EnableRayCastDisplay(bool flag);
    void IncrementScaleFactor();
    void DecrementScaleFactor();
    osg::Geode* CreateLocalFramesNode();
    osg::Geode* CreateVertexNormalsNode();
    unsigned int GenerateComponentId();
    ModelSolver* GetModelSolver();

private:

    void model_update();
    void calculate_ellipse();
    void update_dynamic_segment();
    void update_dynamic_segment_with_mirror_point(const osg::Vec2d& pt, bool first);
    void ray_cast_within_binary_image_for_profile_match();     // based on binary images
    void ray_cast_within_gradient_image_for_profile_match();   // based on gradient
    void initialize_spine_drawing_mode(projection_type pt);
    void constrain_mouse_point();

    // estimation of the other circles
    inline void add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
    inline void add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
    inline void add_planar_section_to_the_generalized_cylinder_under_orthogonality_constraint();

    // estimation of the first circle
    inline void estimate_first_circle_under_persective_projection();
    inline void estimate_first_circle_under_orthographic_projection();
    inline void estimate_first_circle_under_orthogonality_constraint();

    // 3D circle estimation functions
    inline int estimate_3d_circles_with_fixed_radius(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_radius);
    inline int estimate_3d_circles_with_fixed_depth(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_depth);
    inline int estimate_unit_3d_circles(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles);
    inline void estimate_3d_circle_under_orthographic_projection(std::unique_ptr<Ellipse2D>& ellipse, Circle3D& circle);


    // selection of the estimated circles under perspective projection
    inline size_t select_first_3d_circle(const Circle3D* const circles);
    inline size_t select_parallel_circle(const Circle3D* const circles);

    // projection functions
    void project_point(const osg::Vec3d& pt3d, osg::Vec2d& pt2d) const;
    void project_points(const osg::Vec3dArray * const pt3darr, osg::Vec2dArray* const pt2darr) const;
    void project_circle(const Circle3D& circle, Ellipse2D& ellipse) const;
    void project_generalized_cylinder(const GeneralizedCylinder& gcyl) const;

    // test
    void test_circle_estimation_from_major_axis(const Segment2D& seg);
    void test_circle_estimation_from_major_axis();

    // projection error calculators
    // inline void constrain_spine_point_in_piecewise_linear_mode();
    // inline void constrain_spine_point_in_continuous_mode();
};

#endif // IMAGE_MODELLER_HPP
