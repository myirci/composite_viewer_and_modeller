#ifndef IMAGE_MODELLER_HPP
#define IMAGE_MODELLER_HPP

#include <vector>
#include <memory>
#include "UIHelper.hpp"
#include "optimization/ModelSolver.hpp"
#include "../geometry/Ellipse2D.hpp"
#include "../geometry/Rectangle2D.hpp"

#include "components/GeneralizedCylinder.hpp"

#include <osg/Geometry>
#include <otbImage.h>
#include <wx/gdicmn.h>

class Circle3D;
class GeneralizedCylinder;
class PersProjParam;
class OsgWxGLCanvas;

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

class ImageModeller {
public:
    ImageModeller(const std::string& fpath, const std::shared_ptr<PersProjParam>& ppp, OsgWxGLCanvas* canvas);
    ~ImageModeller();

    component_type comp_type;           // type of the component being modelled
    spine_constraints sp_constraints;   // constraints on the main axis of the component
    section_constraints sc_constraints; // constraints on the sections of the generalized cylinders and cuboids
    spine_drawing_mode spd_mode;        // either continuous or piecewise

private:

    // otb related data members
    typedef otb::Image<unsigned char, 2> ImageType;
    ImageType::Pointer m_image;                             // for binary image

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

    Circle3D* m_last_circle;
    std::unique_ptr<UIHelper> m_uihelper;
    std::unique_ptr<Ellipse2D> m_ellipse;
    std::unique_ptr<Ellipse2D> m_last_profile;
    std::unique_ptr<Ellipse2D> m_dynamic_profile;
    std::unique_ptr<Rectangle2D> m_rect;
    std::unique_ptr<ModelSolver> m_solver;
    std::shared_ptr<PersProjParam> m_ppp;

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
    osg::Geode* CreateLocalFramesNode();
    osg::Geode* CreateVertexNormalsNode();
    unsigned int GenerateComponentId();
    ModelSolver* GetModelSolver();

private:

    void model_update();
    void calculate_ellipse();
    void generate_dynamic_profile();
    void ray_cast_for_profile_match();
    void initialize_spine_drawing_mode();
    inline void add_planar_section_to_the_generalized_cylinder();

    // 3D circle estimation functions
    void estimate_3d_circles_with_fixed_radius(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_radius);
    void estimate_3d_circles_with_fixed_depth(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_depth);
    void estimate_unit_3d_circles(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles);

    size_t select_first_3d_circle(const Circle3D* const circles);
    size_t select_parallel_circle(const Circle3D* const circles);
    size_t select_planar_circle(const Circle3D* const circles);

    // void constrain_point(osg::Vec2d& pt);
    // inline void constrain_spine_point_in_piecewise_linear_mode();
    // inline void constrain_spine_point_in_continuous_mode();
};

#endif // IMAGE_MODELLER_HPP
