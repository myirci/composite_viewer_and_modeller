#include "ImageModeller.hpp"
#include "optimization/CircleEstimator.hpp"
#include "optimization/ComponentSolver.hpp"
#include "optimization/ModelSolver.hpp"
#include "ProjectionParameters.hpp"
#include "UIHelper.hpp"
#include "../geometry/Circle3D.hpp"
#include "../geometry/Rectangle2D.hpp"
#include "../geometry/Ellipse2D.hpp"
#include "../osg/OsgWxGLCanvas.hpp"
#include "../osg/OsgUtility.hpp"
#include "../utility/Utility.hpp"
#include "../wx/WxUtility.hpp"
#include "../image/algorithms/Algorithms.hpp"

#include <osgDB/WriteFile>
#include <otbImageFileReader.h>

ImageModeller::ImageModeller(const wxString& fpath, const std::shared_ptr<ProjectionParameters>& pp, OsgWxGLCanvas* canvas) :
    m_pp(pp),
    m_canvas(canvas),
    m_gcyl(nullptr),
    m_vertices(nullptr),
    m_rect(nullptr),
    m_rtype(rendering_type::triangle_strip),
    m_solver(new ModelSolver()),
    m_component_solver(new ComponentSolver(-pp->near)),
    m_mode(drawing_mode::mode_0),
    m_left_click(false),
    m_right_click(false),
    m_bimg_exists(false),
    spd_mode(spine_drawing_mode::piecewise_linear),
    sp_constraints(spine_constraints::none),
    sc_constraints(section_constraints::constant),
    comp_type(component_type::generalized_cylinder),
    m_ellipse(new Ellipse2D),
    m_last_profile(new Ellipse2D),
    m_dynamic_profile(new Ellipse2D),
    m_last_circle(new Circle3D),
    m_first_circle(new Circle3D),
    m_uihelper(nullptr),
    m_tilt_angle(0.0),
    m_circle_estimator(new CircleEstimator),
    m_display_raycast(false),
    m_raycast(nullptr) {

    std::string binary_img_path = utilityInsertAfter(fpath, wxT('.'), wxT("_binary"));
    std::ifstream infile(binary_img_path);
    if(infile.good()) {
        std::cout << "INFO: Binary image is loaded" << std::endl;
        otb::ImageFileReader<OtbImageType>::Pointer reader = otb::ImageFileReader<OtbImageType>::New();
        reader->SetFileName(binary_img_path);
        reader->Update();
        m_bimage = reader->GetOutput();
        m_bimg_exists = true;
        // binary image is ready for processing
        OtbImageType::SizeType size = m_bimage->GetLargestPossibleRegion().GetSize();
        m_rect = std::unique_ptr<Rectangle2D>(new Rectangle2D(0, 0, size[0] - 1, size[1] - 1));
    }
    else {

        std::cout << "INFO: No associated binary image file: " << std::endl;
        std::string grad_img_path = utilityInsertAfter(fpath, wxT('.'), wxT("_grad"));
        std::ifstream grad_file(grad_img_path);
        if(grad_file.good()) {
            m_gimage = LoadImage<OtbImageType>(grad_img_path);
            std::cout << "INFO: Gradient image is loaded" << std::endl;
        }
        else {
            OtbFloatVectorImageType::Pointer img = LoadImage<OtbFloatVectorImageType>(fpath.ToStdString());
            m_gimage = GradientMagnitudeImage(img);
            SaveImage<OtbImageType>(m_gimage, grad_img_path);
            std::cout << "INFO: No gradient image! Gradient image is calculated and saved to the path: " << grad_img_path << std::endl;
        }

        OtbImageType::SizeType size = m_gimage->GetLargestPossibleRegion().GetSize();
        m_rect = std::unique_ptr<Rectangle2D>(new Rectangle2D(0, 0, size[0] - 1, size[1] - 1));
    }

    m_fixed_depth = -(m_pp->near + m_pp->far)/2.0;
}

ImageModeller::~ImageModeller() {
    delete m_last_circle;
}

unsigned int ImageModeller::GenerateComponentId() {

    static unsigned int component_id_source = 0;
    return ++component_id_source;
}

ModelSolver* ImageModeller::GetModelSolver() {
    return m_solver.get();
}

void ImageModeller::Initialize2DDrawingInterface(osg::Geode* geode) {

    m_uihelper = std::unique_ptr<UIHelper>(new UIHelper(geode));
    m_vertices = new osg::Vec2dArray;
}

void ImageModeller::SaveModel(const std::string& path) {

    if(!m_gcyl.valid()) {
        std::cout << "ERROR: Generalized cylinder is not valid" << std::endl;
        return;
    }
    osgDB::writeNodeFile(*m_gcyl, path);
}

void ImageModeller::DeleteModel() {
    m_solver->DeleteAllComponents();
}

void ImageModeller::DeleteSelectedComopnents(std::vector<int>& index_vector) {
    m_solver->DeleteSelectedComponents(index_vector);
}

void ImageModeller::SetRenderingType(rendering_type rtype) {
    m_rtype = rtype;
}

void ImageModeller::EnableRayCastDisplay(bool flag) {

    if(m_display_raycast == flag) return;
    m_display_raycast = flag;
    if(m_display_raycast) m_raycast = new osg::Vec2dArray(8);
    else                  m_raycast.release();
}

osg::Geode* ImageModeller::CreateLocalFramesNode() {

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    osg::Vec3dArray* vertices = new osg::Vec3dArray;
    // fill the vertices here
    vertices->push_back(osg::Vec3d(0,0,-50));
    vertices->push_back(osg::Vec3d(100,100,-50));
    geom->setVertexArray(vertices);
    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    geom->setColorArray(colors, osg::Array::BIND_OVERALL);
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, 2));
    geode->addDrawable(geom.get());
    return geode.release();
}

osg::Geode* ImageModeller::CreateVertexNormalsNode() {

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    osg::Vec3dArray* vertices = new osg::Vec3dArray;
    // fill the vertices here
    vertices->push_back(osg::Vec3d(0,0,-50));
    vertices->push_back(osg::Vec3d(100,100,-50));
    geom->setVertexArray(vertices);
    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    geom->setColorArray(colors, osg::Array::BIND_OVERALL);
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, 2));
    geode->addDrawable(geom.get());
    return geode.release();
}

void ImageModeller::DebugPrint() {
    m_solver->Print();
}

void ImageModeller::Reset2DDrawingInterface() {

    m_mode = drawing_mode::mode_0;
    m_left_click = false;
    m_right_click = false;
    m_vertices->clear();
    m_uihelper->Reset();
}

void ImageModeller::OnLeftClick(double x, double y) {

    m_left_click = true;
    m_mouse.set(x,y);
    m_vertices->push_back(m_mouse);
    model_update();
}

void ImageModeller::OnRightClick(double x, double y) {

    if(m_mode == drawing_mode::mode_3) {
        m_right_click = true;
        m_mouse.set(x,y);
        m_vertices->push_back(m_mouse);
        model_update();
    }
}

void ImageModeller::OnMouseMove(double x, double y) {

    if(m_mode != drawing_mode::mode_0) {
        m_mouse.set(x,y);
        model_update();
    }
}

// Execution of the modelling process is done within this function.
void ImageModeller::model_update() {

    if(comp_type == component_type::generalized_cylinder) {

        if(m_mode == drawing_mode::mode_0) {
            if(m_left_click) {
                // first click
                m_left_click = false;
                m_uihelper->InitializeMajorAxisDrawing(m_mouse);
                m_mode = drawing_mode::mode_1;
            }
        }
        else if(m_mode == drawing_mode::mode_1) {
            if(m_left_click) {
                // second click: major axis has been determined.
                m_left_click = false;
                m_ellipse->update_major_axis(m_vertices->at(0), m_vertices->at(1));
                m_uihelper->InitializeMinorAxisDrawing(m_mouse);
                m_mode = drawing_mode::mode_2;
            }
            else {
                // Here we are executing the major axis drawing mode. In this mode we only update the end point of
                // the major axis with the mouse position. The operator has not decided the major axis yet.
                m_uihelper->Updatep1(m_mouse);
            }
        }
        else if(m_mode == drawing_mode::mode_2) {

            // calculate the possible ellpise based on the current mouse position
            calculate_ellipse();

            if(m_left_click) {
                // third click: base ellipse (m_ellipse) has been determined.
                m_left_click = false;
                initialize_spine_drawing_mode();
                m_uihelper->InitializeSpineDrawing(m_ellipse);
                m_mode = drawing_mode::mode_3;
            }
        }

        else if(m_mode == drawing_mode::mode_3) {

            if(spd_mode == spine_drawing_mode::continuous) {

                if(m_left_click) {
                    m_left_click = false;
                    // End the modelling process for the current generalized cylinder with the 4th click. The last
                    // clicked point is accepted as the last sample point.
                    generate_dynamic_profile();
                    *m_last_profile = *m_dynamic_profile;
                    add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                    // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();

                    m_solver->AddComponent(m_gcyl.get());
                    Reset2DDrawingInterface();
                }
                else {
                    // Here we are executing the mode_3 (spine drawing mode) for continuos spine drawing. Generate the
                    // 2D profile along the path of the spine as the spine is being drawn

                    // vector from last validated spine point to the current mouse position
                    osg::Vec2d vec = m_mouse - m_last_profile->points[2];

                    // If the distance between the last validated spine point and the candidate spine point (mouse point)
                    // is bigger than a threshold, then the current spine point is accepted as a sample point.
                    if(vec.length2() > 100) {
                        generate_dynamic_profile();
                        *m_last_profile = *m_dynamic_profile;
                        // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                        add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                    }
                }
            }
            else if(spd_mode == spine_drawing_mode::piecewise_linear) {

                // estimate the elliptic dynamic profile
                generate_dynamic_profile();

                if(m_right_click) {
                    // right click ends the modelling of the current component being modelled.
                    m_right_click = false;
                    m_uihelper->AddSpinePoint(m_mouse);
                    *m_last_profile = *m_dynamic_profile;
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                    // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                    add_planar_section_to_the_generalized_cylinder_constrained();
                    m_solver->AddComponent(m_gcyl.get());
                    Reset2DDrawingInterface();

                    project_generalized_cylinder(*m_gcyl);

                }
                else if(m_left_click) {

                    // left click is a new spine point
                    m_left_click = false;
                    m_uihelper->AddSpinePoint(m_mouse);
                    *m_last_profile = *m_dynamic_profile;
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                    // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                    add_planar_section_to_the_generalized_cylinder_constrained();
                }
                else {
                    m_uihelper->SpinePointCandidate(m_mouse);
                    m_uihelper->UpdateSweepCurve(m_dynamic_profile);
                }
            }
        }
    }
}

void ImageModeller::estimate_first_circle_under_persective_projection() {

    // 1) Estimate the first 3D circle under perspective projection from the ellipse drawn by the operator
    // (m_ellipse): projection of two circles will match with the ellipse.
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_ellipse, circles, m_fixed_depth);
    if(count == 2) {
        // 1.1) Select one of the two estimated circles based on how the user drew the ellipse
        *m_last_circle = circles[select_first_3d_circle(circles)];
    }
    else if (count == 1) {
        *m_last_circle = circles[0];
    }
    else {
        std::cout << "ERROR: Perspective 3D circle estimation error "  << std::endl;
    }

    // 2) Calculate the tilt-angle
    m_tilt_angle = acos(m_last_circle->normal[2]);
    if(m_tilt_angle > HALF_PI)
        m_tilt_angle = PI - m_tilt_angle;

}

void ImageModeller::estimate_first_circle_under_persective_projection__() {

    // 1) Estimate the 3D circles under perspective projection for the initial user drawn ellipse (m_ellipse): projection of two circles
    //    will match with the ellipse.
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_ellipse, circles, m_fixed_depth);
    if(count == 2) {
        // 1.1) Select one of the two estimated circles based on how the user drew the ellipse
        *m_last_circle = circles[select_first_3d_circle(circles)];
    }
    else if (count == 1) {
        *m_last_circle = circles[0];
    }
    else {
        std::cout << "ERROR: Perspective 3D circle estimation error "  << std::endl;
    }

    // 2) Calculate the tilt-angle
    m_tilt_angle = acos(m_last_circle->normal[2]);
    if(m_tilt_angle > HALF_PI)
        m_tilt_angle -= HALF_PI;
}

void ImageModeller::estimate_first_circle_under_orthographic_projection() {

    // 1) Estimate the 3D circles under orthographic projection
    // estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(m_ellipse, *m_last_circle, m_fixed_depth);
    estimate_3d_circles_under_orthographic_projection(m_ellipse, *m_first_circle);
    // estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(m_ellipse, *m_last_circle, m_fixed_depth);

    // 2) Calculate the tilt angle
    m_tilt_angle = acos(m_ellipse->smn_axis / m_ellipse->smj_axis);
}
/*
void ImageModeller::estimate_first_circle_under_orthographic_projection() {

    // 1) Estimate the 3D circles under orthographic projection
    // estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(m_ellipse, *m_last_circle, m_fixed_depth);
    estimate_3d_circles_under_orthographic_projection(m_ellipse, *m_last_circle);
    // estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(m_ellipse, *m_last_circle, m_fixed_depth);

    // 2) Calculate the tilt angle
    m_tilt_angle = acos(m_ellipse->smn_axis / m_ellipse->smj_axis);
}
*/

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection() {

    // 1) Estimate the 3D circles for the current profile and add it to the generatlized cylinder.
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_last_profile, circles, m_fixed_depth);
    if(count == 2) {
        // 1.1) Select one of the two estimated circles based on the angle between the normals
        *m_last_circle = circles[select_parallel_circle(circles)];
    }
    else if (count == 1) {
        *m_last_circle = circles[0];
    }
    else {
        std::cout << "ERROR: Perspective 3D circle estimation error "  << std::endl;
    }

    // 2) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_orthographic_projection() {

    // 1) Estimate the 3D circles under orthographic projection
    // estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(m_last_profile, *m_last_circle, m_fixed_depth);
    estimate_3d_circles_under_orthographic_projection(m_last_profile, *m_last_circle);
    // estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(m_last_profile, *m_last_circle, m_fixed_depth);

    // 2) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_constrained() {

    // transform m_last_profile to projected coordinates
    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*m_last_profile, elp_prj, false);

    osg::Vec2dArray* arr = m_component_solver->GetLocalCoordinateFrameProjections();
    arr->at(3) = elp_prj.center;

    // 1) set the radius : length of the semi-major axis
    m_last_circle->radius = elp_prj.smj_axis;

    // 2) set the center
    m_last_circle->center[0] = elp_prj.center.x();
    m_last_circle->center[1] = elp_prj.center.y();
    m_last_circle->center[2] = -m_pp->near;

    // 3) set the normal : tilt angle (constant for a generalized cylinder) & bend angle (rot angle for the ellipse)
    m_last_circle->normal[0] = sin(m_tilt_angle) * sin(m_last_profile->rot_angle);
    m_last_circle->normal[1] = -sin(m_tilt_angle) * cos(m_last_profile->rot_angle);

    osg::Vec2d smj_vec = elp_prj.points[1] - elp_prj.center;
    osg::Vec2d smn_vec = elp_prj.points[2] - elp_prj.center;
    if(smj_vec.x() * smn_vec.y() - smj_vec.y() * smn_vec.x() > 0)
        m_last_circle->normal[2] = cos(m_tilt_angle);
    else
        m_last_circle->normal[2] = -cos(m_tilt_angle);

    double calculated_depth = -m_pp->near;
    double C1[3];
    double C2[3];

    C1[0] = m_first_circle->center[0];
    C1[1] = m_first_circle->center[1];
    C1[2] = m_first_circle->center[2];

    C2[0] = m_last_circle->center[0];
    C2[1] = m_last_circle->center[1];
    C2[2] = m_last_circle->center[2];

    double radius_1 = m_component_solver->Solve(calculated_depth, C1, C2);

    m_first_circle->center[0] = C1[0];
    m_first_circle->center[1] = C1[1];
    m_first_circle->center[2] = C1[2];

    m_last_circle->center[0] = C2[0];
    m_last_circle->center[1] = C2[1];
    m_last_circle->center[2] = C2[2];

    osg::Vec3d nrm(C1[0] - C2[0], C1[1] - C2[1], C1[1] - C2[1]);
    nrm.normalize();

    m_last_circle->normal[0] = nrm.x();
    m_last_circle->normal[1] = nrm.y();
    m_last_circle->normal[2] = nrm.z();

    m_first_circle->normal[0] = nrm.x();
    m_first_circle->normal[1] = nrm.y();
    m_first_circle->normal[2] = nrm.z();

    m_first_circle->radius = radius_1;

    m_gcyl->AddPlanarSection(*m_first_circle);
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

/*
void ImageModeller::add_planar_section_to_the_generalized_cylinder_constrained() {

    // 1) set the radius : proportional to the length of the semi-major axis
    // m_last_circle->radius = m_last_profile->smj_axis * (m_fixed_depth / -m_pp->near);

    // 2) set the center : should be scaled with respect to the fixed depth
    osg::Vec2d ctr;
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->center, ctr);
    m_last_circle->center[0] = ctr.x();
    m_last_circle->center[1] = ctr.y();
    m_last_circle->center[2] = -m_pp->near;
    m_last_circle->center *= (m_fixed_depth / -m_pp->near);

    // 3) set the normal : tilt angle (constant for a generalized cylinder) & bend angle (rot angle for the ellipse)
    m_last_circle->normal[0] = sin(m_tilt_angle) * sin(m_last_profile->rot_angle);
    m_last_circle->normal[1] = -sin(m_tilt_angle) * cos(m_last_profile->rot_angle);

    osg::Vec2d smj_vec = m_last_profile->points[1] - m_last_profile->center;
    osg::Vec2d smn_vec = m_last_profile->points[2] - m_last_profile->center;
    if(smj_vec.x() * smn_vec.y() - smj_vec.y() * smn_vec.x() > 0)
        m_last_circle->normal[2] = cos(m_tilt_angle);
    else
        m_last_circle->normal[2] = -cos(m_tilt_angle);

    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}


void ImageModeller::initialize_spine_drawing_mode() {

    // modify the user clicked point with its projection on the minor-axis guide line
    m_vertices->at(2) = m_ellipse->points[2];

    // initialize the generalized cylinder as a new node in the scene graph.
    if(m_gcyl.valid()) m_gcyl = nullptr;

    // estimate the first circle
    // estimate_first_circle_under_persective_projection();
    estimate_first_circle_under_orthographic_projection();

    Ellipse2D elp;
    project_circle(*m_last_circle, elp);

    // m_gcyl = new GeneralizedCylinder(GenerateComponentId(), *m_last_circle);
    m_gcyl = new GeneralizedCylinder(GenerateComponentId(), *m_last_circle);
    m_canvas->UsrAddSelectableNodeToDisplay(m_gcyl.get(), m_gcyl->GetComponentId());

    // initialize the constraint line if spine contraint is straight_planar
    if(sp_constraints == spine_constraints::straight_planar) {

//        // let's define the constraint line in the opengl screen coordinate which is also
//        double a = m_last_circle->center[1] * m_last_circle->normal[2] - m_last_circle->center[2] * m_last_circle->normal[1];
//        double b = m_last_circle->center[2] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[2];
//        double c = (m_last_circle->center[1] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[1]) * m_pp->height / (2*tan(m_pp->fovy/2.0)) - a*(m_pp->width/2.0) - b*(m_pp->height/2.0);
//        m_constraint_line = std::unique_ptr<Line2D>(new Line2D(a,b,c));

//        // display the constraint line
//        std::vector<osg::Vec2d> intersection_pts;
//        double var = m_constraint_line->get_x_at_y(0);
//        if(var >= 0 && var < m_pp->width)
//            intersection_pts.push_back(osg::Vec2d(var, 0));

//        var = m_constraint_line->get_x_at_y(m_pp->height - 1);
//        if(var >= 0 && var < m_pp->width)
//            intersection_pts.push_back(osg::Vec2d(var, m_pp->height - 1));

//        var = m_constraint_line->get_y_at_x(0);
//        if(var >= 0 && var < m_pp->height)
//            intersection_pts.push_back(osg::Vec2d(0, var));

//        var = m_constraint_line->get_y_at_x(m_pp->width - 1);
//        if(var >= 0 && var < m_pp->height)
//            intersection_pts.push_back(osg::Vec2d(m_pp->width - 1, var));

//        if(intersection_pts.size() == 2)
//            m_uihelper->DisplayConstraintLine(intersection_pts);
//        else
//            std::cout << "ERROR: Number of intersections must be two!" << std::endl;

//        std::cout << "intersection points: \n";
//        for(auto it = intersection_pts.begin(); it != intersection_pts.end(); ++it)
//            std::cout << it->x() << " " << it->y() << std::endl;

    }

    // Copy the base ellipse into the m_last_profile.
    *m_last_profile = *m_ellipse;
}

*/

void ImageModeller::initialize_spine_drawing_mode() {

    // modify the user clicked point with its projection on the minor-axis guide line
    m_vertices->at(2) = m_ellipse->points[2];

    // initialize the generalized cylinder as a new node in the scene graph.
    if(m_gcyl.valid()) m_gcyl = nullptr;

    // estimate the first circle
    // estimate_first_circle_under_persective_projection();
    estimate_first_circle_under_orthographic_projection();

    Ellipse2D elp;
    project_circle(*m_first_circle, elp);

    m_gcyl = new GeneralizedCylinder(GenerateComponentId());
    m_canvas->UsrAddSelectableNodeToDisplay(m_gcyl.get(), m_gcyl->GetComponentId());

    // initialize the constraint line if spine contraint is straight_planar
    if(sp_constraints == spine_constraints::straight_planar) {

//        // let's define the constraint line in the opengl screen coordinate which is also
//        double a = m_last_circle->center[1] * m_last_circle->normal[2] - m_last_circle->center[2] * m_last_circle->normal[1];
//        double b = m_last_circle->center[2] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[2];
//        double c = (m_last_circle->center[1] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[1]) * m_pp->height / (2*tan(m_pp->fovy/2.0)) - a*(m_pp->width/2.0) - b*(m_pp->height/2.0);
//        m_constraint_line = std::unique_ptr<Line2D>(new Line2D(a,b,c));

//        // display the constraint line
//        std::vector<osg::Vec2d> intersection_pts;
//        double var = m_constraint_line->get_x_at_y(0);
//        if(var >= 0 && var < m_pp->width)
//            intersection_pts.push_back(osg::Vec2d(var, 0));

//        var = m_constraint_line->get_x_at_y(m_pp->height - 1);
//        if(var >= 0 && var < m_pp->width)
//            intersection_pts.push_back(osg::Vec2d(var, m_pp->height - 1));

//        var = m_constraint_line->get_y_at_x(0);
//        if(var >= 0 && var < m_pp->height)
//            intersection_pts.push_back(osg::Vec2d(0, var));

//        var = m_constraint_line->get_y_at_x(m_pp->width - 1);
//        if(var >= 0 && var < m_pp->height)
//            intersection_pts.push_back(osg::Vec2d(m_pp->width - 1, var));

//        if(intersection_pts.size() == 2)
//            m_uihelper->DisplayConstraintLine(intersection_pts);
//        else
//            std::cout << "ERROR: Number of intersections must be two!" << std::endl;

//        std::cout << "intersection points: \n";
//        for(auto it = intersection_pts.begin(); it != intersection_pts.end(); ++it)
//            std::cout << it->x() << " " << it->y() << std::endl;

    }

    // Copy the base ellipse into the m_last_profile.
    *m_last_profile = *m_ellipse;
}

void ImageModeller::calculate_ellipse() {

    // calculate the end points of the minor_axis guide line
    osg::Vec2d vec_mj = m_ellipse->points[1] - m_ellipse->points[0];
    osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
    vec_mn.normalize();
    osg::Vec2d pt_0 = m_ellipse->center - vec_mn * m_ellipse->smj_axis;
    osg::Vec2d pt_1 = m_ellipse->center + vec_mn * m_ellipse->smj_axis;

    osg::Vec2d vec1 = m_mouse - pt_0;
    osg::Vec2d vec2 = pt_1 - pt_0;

    double ratio = (vec1 * vec2) / (vec2 * vec2) ;
    if(ratio >= 0.0 && ratio <= 1.0) {

        // find the projection point
        vec2.normalize();
        osg::Vec2d proj_point = (vec2 * 2 * ratio * m_ellipse->smj_axis) + pt_0;

        // update the m_ellipse
        m_ellipse->update_minor_axis(proj_point);

        // display the ellipse and a small circle on the projection point
        m_uihelper->UpdateBaseEllipse(m_ellipse);
    }
}

void ImageModeller::generate_dynamic_profile() {

    // copy the last profile into the dynamic profile
    *m_dynamic_profile = *m_last_profile;

    // vector from last validated spine point to the current mouse position
    osg::Vec2d vec1 = m_mouse - m_last_profile->points[2];

    if(sp_constraints == spine_constraints::straight_planar) {
        // translate the profile to the current spine point
        m_dynamic_profile->translate(vec1);
    }
    else {
        // rotate the dynamic ellipse according to the bend of the spine curve
        osg::Vec2d vec2 = m_last_profile->points[3] - m_last_profile->points[2];
        vec2.normalize();
        double angle = acos((vec1 * vec2) / vec1.length());

        // determine the orientation of the rotation: CW (negative angles) or CCW(positive angles)
        if(angle > HALF_PI) angle -= PI;
        if(vec2.x()*vec1.y() - vec2.y()*vec1.x() > 0) m_dynamic_profile->rotate(angle);
        else m_dynamic_profile->rotate(-angle);

        // translate the rotated ellipse such that p2 concides back to the mouse point
        // m_dynamic_profile->translate(m_mouse - m_dynamic_profile->points[2]);
        m_dynamic_profile->translate(m_mouse - m_dynamic_profile->center);
    }

    if(m_bimg_exists) ray_cast_within_binary_image_for_profile_match();
    else              ray_cast_within_gradient_image_for_profile_match();

}

int ImageModeller::estimate_3d_circles_with_fixed_radius(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_radius) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    return m_circle_estimator->estimate_3d_circles_with_fixed_radius(elp_prj, circles, m_pp.get(), desired_radius);
}

int ImageModeller::estimate_3d_circles_with_fixed_depth(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_depth) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    return m_circle_estimator->estimate_3d_circles_with_fixed_depth(elp_prj, circles, m_pp.get(), desired_depth);
}

int ImageModeller::estimate_3d_circles_with_fixed_depth__(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_depth) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    return m_circle_estimator->estimate_3d_circles_with_fixed_depth(elp_prj, circles, m_pp.get(), desired_depth);
}

int ImageModeller::estimate_unit_3d_circles(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    return m_circle_estimator->estimate_unit_3d_circles(elp_prj, circles, m_pp.get());
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection(std::unique_ptr<Ellipse2D>& ellipse, Circle3D& circle) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);

    osg::Vec2dArray* arr = m_component_solver->GetLocalCoordinateFrameProjections();
    arr->at(0) = elp_prj.points[2]; // third clicked point
    arr->at(1) = elp_prj.points[0]; // first clicked point
    arr->at(2) = elp_prj.points[1]; // third clicked point

    m_circle_estimator->estimate_3d_circles_under_orthographic_projection(elp_prj, circle, m_pp->near);
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(std::unique_ptr<Ellipse2D>& ellipse, Circle3D &circle, double desired_depth) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    m_circle_estimator->estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(elp_prj, circle, m_pp->near, desired_depth);
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(std::unique_ptr<Ellipse2D>& ellipse, Circle3D& circle, double desired_depth) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    m_circle_estimator->estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(elp_prj, circle, m_pp->near, desired_depth);
    std::cout << circle.normal << std::endl;
}

size_t ImageModeller::select_first_3d_circle(const Circle3D* const circles) {

    // get projection and viewport mapping matrix
    const osg::Camera* const cam = m_canvas->UsrGetMainCamera();
    osg::Matrix projectionMatrix = cam->getProjectionMatrix();
    osg::Matrix windowMatrix = cam->getViewport()->computeWindowMatrix();

    // two point on the normal of the circle
    osg::Vec3 ctr(circles[0].center[0], circles[0].center[1], circles[0].center[2]);
    Eigen::Vector3d normal_tip = circles[0].center + circles[0].normal;
    osg::Vec3 ntip(normal_tip[0], normal_tip[1], normal_tip[2]);

    // project two points
    ctr = ctr * projectionMatrix * windowMatrix;
    ntip = ntip * projectionMatrix * windowMatrix;

    // construct 2D vectors:
    Vector2D<double> vec1(ntip.x() - ctr.x(), ntip.y() - ctr.y());
    Vector2D<double> vec2(ctr.x() - m_ellipse->points[2].x(), ctr.y() - m_ellipse->points[2].y());

    return (vec1.dot(vec2) < 0) ? 0 : 1;
}

size_t ImageModeller::select_parallel_circle(const Circle3D * const circles) {

    double a = (m_last_circle->normal).dot(circles[0].normal);
    double b = (m_last_circle->normal).dot(circles[1].normal);
    return std::abs(a) > std::abs(b) ? 0 : 1;
}

void ImageModeller::project_point(const osg::Vec3d& pt3d, osg::Vec2d& pt2d) const {

    osg::Vec4d pt4d(pt3d,1);
    const osg::Matrixd& Mproj = m_canvas->UsrGetMainCamera()->getProjectionMatrix();
    osg::Matrixd Mvp = m_canvas->UsrGetMainCamera()->getViewport()->computeWindowMatrix();
    osg::Vec4d pt_vp = pt4d * Mproj * Mvp;
    pt2d.x() = pt_vp.x() / pt_vp.w();
    pt2d.y() = pt_vp.y() / pt_vp.w();
}

void ImageModeller::project_points(const osg::Vec3dArray * const pt3darr, osg::Vec2dArray * const pt2darr) const {

    osg::Matrixd M = m_canvas->UsrGetMainCamera()->getProjectionMatrix() *
                     m_canvas->UsrGetMainCamera()->getViewport()->computeWindowMatrix();
    osg::Vec4d pt4d;
    for(auto it = pt3darr->begin(); it != pt3darr->end(); ++it) {
        pt4d = osg::Vec4d(*it, 1) * M;
        pt2darr->push_back(osg::Vec2d(pt4d.x()/pt4d.w(), pt4d.y()/pt4d.w()));
    }
}

void ImageModeller::project_circle(const Circle3D& circle, Ellipse2D& ellipse) const {

    osg::Matrixd M = m_canvas->UsrGetMainCamera()->getProjectionMatrix();
    transpose(M,4);
    Eigen::Matrix<double, 3, 4> Mprj;
    for(size_t i = 0; i < 4; ++i) {
        Mprj(0,i) = M(0,i);
        Mprj(1,i) = M(1,i);
        Mprj(2,i) = M(3,i);
    }

    Eigen::Matrix4d Qs;
    circle.get_matrix_representation(Qs);
    Eigen::Matrix3d dConic = Mprj * Qs * Mprj.transpose();

    if(dConic.determinant() != 0) {
        Eigen::Matrix3d conic = dConic.inverse();
        conic /= conic(0,0);
        ellipse.coeff[0] = conic(0,0);
        ellipse.coeff[1] = 2*conic(0,1);
        ellipse.coeff[2] = conic(1,1);
        ellipse.coeff[3] = 2*conic(0,2);
        ellipse.coeff[4] = 2*conic(1,2);
        ellipse.coeff[5] = conic(2,2);
        ellipse.calculate_parameters_from_coeffients();
    }
    else {
        std::cout << "camera::project_camera_circle3d: dual conic is degenerate" << std::endl;
    }

    // Ellipse is in projected coordinatres. We need to convert it to the logical device coordinates
    osg::Matrixd Mvp = m_canvas->UsrGetMainCamera()->getViewport()->computeWindowMatrix();

    osg::Vec4d pt4d;
    for(size_t i = 0; i < 4; ++i) {
        pt4d = osg::Vec4d(ellipse.points[i].x(), ellipse.points[i].y(), -m_pp->near, 1) * Mvp;
        ellipse.points[i].x() = pt4d.x() / pt4d.w();
        ellipse.points[i].y() = pt4d.y() / pt4d.w();
    }

    // calculate the remaining parameters and the coefficients
    ellipse.center.x() = (ellipse.points[0].x() + ellipse.points[1].x())/2.0;
    ellipse.center.y() = (ellipse.points[0].y() + ellipse.points[1].y())/2.0;
    ellipse.smj_axis = (ellipse.points[0] - ellipse.center).length();
    ellipse.smn_axis = (ellipse.points[2] - ellipse.center).length();
    ellipse.calculate_coefficients_from_parameters();
}

void ImageModeller::project_generalized_cylinder(const GeneralizedCylinder& gcyl)  const {

    std::vector<Ellipse2D> ellipses;
    const std::vector<Circle3D>& circles =  gcyl.GetGeometry()->GetSections();
    osg::ref_ptr<osg::Vec3dArray> main_axis = new osg::Vec3dArray;
    Ellipse2D elp;
    for(size_t i = 0; i < circles.size(); ++i) {
        main_axis->push_back(osg::Vec3d(circles[i].center[0], circles[i].center[1], circles[i].center[2]));
        project_circle(circles[i], elp);
        ellipses.push_back(elp);
    }

    osg::ref_ptr<osg::Vec2dArray> main_axis_prj = new osg::Vec2dArray;
    project_points(main_axis.get(), main_axis_prj.get());

    osg::Vec2d left, right, dir;
    std::vector<osg::Vec2d> left_silhouette;
    std::vector<osg::Vec2d> right_silhouette;
    for(int i = 0; i < main_axis_prj->size()-1; ++i) {
        dir = main_axis_prj->at(i+1) - main_axis_prj->at(i);
        ellipses[i].get_tangent_points(dir, left, right);
        left_silhouette.push_back(left);
        right_silhouette.push_back(right);
    }
    ellipses.back().get_tangent_points(dir,left, right);
    left_silhouette.push_back(left);
    right_silhouette.push_back(right);


    // displat the silhouettes
    m_uihelper->DisplayLineStrip(left_silhouette, osg::Vec4(0,1,0,1));
    m_uihelper->DisplayLineStrip(right_silhouette, osg::Vec4(0,0,1,1));

    // display feature curves
    std::vector<osg::Vec2d> elp_pts;
    for(int i = 0; i < ellipses.size(); ++ i) {
        ellipses[i].generate_points_on_the_ellipse(elp_pts, 40);
        m_uihelper->DisplayLineLoop(elp_pts, osg::Vec4(0.1,0.1,0.1,1));
    }

}

void ImageModeller::ray_cast_within_binary_image_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p0(static_cast<int>(m_dynamic_profile->points[0].x()), static_cast<int>(m_dynamic_profile->points[0].y()));
    Point2D<int> p1(static_cast<int>(m_dynamic_profile->points[1].x()), static_cast<int>(m_dynamic_profile->points[1].y()));
    m_canvas->UsrDeviceToLogical(p0);
    m_canvas->UsrDeviceToLogical(p1);

    // 2) calculate the casting direction vector,the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p1.x - p0.x) * 0.2), static_cast<int>((p1.y - p0.y) * 0.2));

    // 3) perform ray casts and display shot rays variables for ray casting
    bool hit_result[4] = { false, false, false, false };
    Point2D<int> hit[4];

    // 3.1) ray cast from p0-center direction
    Point2D<int> end = p0 + dir_vec;
    if(m_rect->intersect(p0, end))
        hit_result[0] = BinaryImageRayCast(m_bimage, p0, end, hit[0]);

    // 3.2) ray cast from p0-outside direction
    end = p0 - dir_vec;
    if(m_rect->intersect(p0, end))
        hit_result[1] = BinaryImageRayCast(m_bimage, p0, end, hit[1]);

    // 3.3) ray cast from p1-center direction
    end = p1 - dir_vec;
    if(m_rect->intersect(p1, end))
        hit_result[2] = BinaryImageRayCast(m_bimage, p1, end, hit[2]);

    // 3.4) ray cast from p1-outside direction
    end = p1 + dir_vec;
    if(m_rect->intersect(p1, end))
        hit_result[3] = BinaryImageRayCast(m_bimage, p1, end, hit[3]);

    // 4) analyze the result of the ray casts
    Point2D<int> p0_hit = p0;
    if(hit_result[0] && hit_result[1])  p0_hit = p0;
    else if(hit_result[0])              p0_hit = hit[0];
    else if(hit_result[1])              p0_hit = hit[1];
    m_canvas->UsrDeviceToLogical(p0_hit);

    Point2D<int> p1_hit;
    if(hit_result[2] && hit_result[3])  p1_hit = p1;
    else if(hit_result[2])              p1_hit = hit[2];
    else if(hit_result[3])              p1_hit = hit[3];
    m_canvas->UsrDeviceToLogical(p1_hit);

    if((hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {          // update p0, p1
        m_dynamic_profile->points[0] = osg::Vec2d(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = osg::Vec2d(p1_hit.x, p1_hit.y);
    }
    else if((hit_result[0] || hit_result[1]) && !(hit_result[2] || hit_result[3])) {    // update p0, p1 is the mirror of p0
        osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = m_dynamic_profile->points[1] + (m_dynamic_profile->points[0] - new_p0);
        m_dynamic_profile->points[0] = new_p0;
    }
    else if(!(hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {    // update p1, p0 is the mirror of p1
        osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
        m_dynamic_profile->points[0] = m_dynamic_profile->points[0] + (m_dynamic_profile->points[1] - new_p1);
        m_dynamic_profile->points[1] = new_p1;
    }
    else {
        return;
    }

    m_dynamic_profile->center = (m_dynamic_profile->points[0] + m_dynamic_profile->points[1]) / 2.0;    // update center
    double smj_new = (m_dynamic_profile->points[0] - m_dynamic_profile->center).length();               // new smj
    m_dynamic_profile->smn_axis *= (smj_new / m_dynamic_profile->smj_axis);                             // update_smn
    m_dynamic_profile->smj_axis = smj_new;                                                              // update_smj
    osg::Vec2d vec_mn(m_dynamic_profile->points[0].y() - m_dynamic_profile->points[1].y(),
                      m_dynamic_profile->points[1].x() - m_dynamic_profile->points[0].x());
    vec_mn.normalize();                                                                                 // normalized minor axis vector
    m_dynamic_profile->points[2] = m_dynamic_profile->center - vec_mn * m_dynamic_profile->smn_axis;    // update p2
    m_dynamic_profile->points[3] = m_dynamic_profile->center + vec_mn * m_dynamic_profile->smn_axis;    // update p3
}

void ImageModeller::ray_cast_within_gradient_image_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p0(static_cast<int>(m_dynamic_profile->points[0].x()), static_cast<int>(m_dynamic_profile->points[0].y()));
    Point2D<int> p1(static_cast<int>(m_dynamic_profile->points[1].x()), static_cast<int>(m_dynamic_profile->points[1].y()));
    m_canvas->UsrDeviceToLogical(p0);
    m_canvas->UsrDeviceToLogical(p1);

    // 2) calculate the casting direction vector, the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p1.x - p0.x) * 0.35), static_cast<int>((p1.y - p0.y) * 0.35));

    // 3) perform ray casts and display shot rays variables for ray casting
    OtbImageType::PixelType hit_val[4];
    Point2D<int> hit_idx[4];

    // 3.1) ray cast from p0-center direction
    Point2D<int> end = p0 + dir_vec;
    if(m_rect->intersect(p0, end))
        hit_val[0] = GradientImageRayCast(m_gimage, p0, end, hit_idx[0]);

    if(m_display_raycast) {
        m_raycast->at(0).x() = end.x;
        m_raycast->at(0).y() = end.y;
    }

    // 3.2) ray cast from p0-outside direction
    end = p0 - dir_vec;
    if(m_rect->intersect(p0, end))
        hit_val[1] = GradientImageRayCast(m_gimage, p0, end, hit_idx[1]);

    if(m_display_raycast) {
        m_raycast->at(1).x() = end.x;
        m_raycast->at(1).y() = end.y;
    }

    // 3.3) ray cast from p1-center direction
    end = p1 - dir_vec;
    if(m_rect->intersect(p1, end))
        hit_val[2] = GradientImageRayCast(m_gimage, p1, end, hit_idx[2]);

    if(m_display_raycast) {
        m_raycast->at(2).x() = end.x;
        m_raycast->at(2).y() = end.y;
    }

    // 3.4) ray cast from p1-outside direction
    end = p1 + dir_vec;
    if(m_rect->intersect(p1, end))
        hit_val[3] = GradientImageRayCast(m_gimage, p1, end, hit_idx[3]);

    if(m_display_raycast) {
        m_raycast->at(3).x() = end.x;
        m_raycast->at(3).y() = end.y;

        for(size_t i = 0; i < 4; ++i) {
            m_raycast->at(i+4).x() = hit_idx[i].x;
            m_raycast->at(i+4).y() = hit_idx[i].y;
        }

        for(size_t i = 0; i < 8; ++i)
            m_canvas->UsrDeviceToLogical(m_raycast->at(i));

        m_uihelper->DisplayRayCast(m_raycast);
    }

    // 4) analyze the result of the ray casts
    OtbImageType::IndexType p0Idx, p1Idx;
    p0Idx[0] = p0.x;
    p0Idx[1] = p0.y;
    OtbImageType::PixelType p0val = m_gimage->GetPixel(p0Idx);
    p1Idx[0] = p1.x;
    p1Idx[1] = p1.y;
    OtbImageType::PixelType p1val = m_gimage->GetPixel(p1Idx);

    bool p0_updated = false;
    Point2D<int> p0_hit = p0;
    if(hit_val[0] > p0val && hit_val[1] > p0val) {
        if(hit_val[0] > hit_val[1]) p0_hit = hit_idx[0];
        else                        p0_hit = hit_idx[1];
        p0_updated = true;
    }
    else if(hit_val[0] > p0val) {
        p0_hit = hit_idx[0];
        p0_updated = true;
    }
    else if(hit_val[1] > p0val) {
        p0_hit = hit_idx[1];
        p0_updated = true;
    }
    m_canvas->UsrDeviceToLogical(p0_hit);

    bool p1_updated = false;
    Point2D<int> p1_hit = p1;
    if(hit_val[2] > p1val && hit_val[3] > p1val) {
        if(hit_val[2] > hit_val[3]) p1_hit = hit_idx[2];
        else                        p1_hit = hit_idx[3];
        p1_updated = true;
    }
    else if(hit_val[2] > p1val) {
        p1_hit = hit_idx[2];
        p1_updated = true;
    }
    else if(hit_val[3] > p1val) {
        p1_hit = hit_idx[3];
        p1_updated = true;
    }
    m_canvas->UsrDeviceToLogical(p1_hit);

    if(p0_updated && p1_updated) { // update p0, p1
        m_dynamic_profile->points[0] = osg::Vec2d(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = osg::Vec2d(p1_hit.x, p1_hit.y);
    }
    else if(p0_updated) { // update p0, p1 is the mirror of p0
        osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = m_dynamic_profile->points[1] + (m_dynamic_profile->points[0] - new_p0);
        m_dynamic_profile->points[0] = new_p0;
    }
    else if(p1_updated) { // update p1, p0 is the mirror of p1
        osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
        m_dynamic_profile->points[0] = m_dynamic_profile->points[0] + (m_dynamic_profile->points[1] - new_p1);
        m_dynamic_profile->points[1] = new_p1;
    }
    // not hits
    else {
        return;
    }

    m_dynamic_profile->center = (m_dynamic_profile->points[0] + m_dynamic_profile->points[1]) / 2.0;    // update center
    double smj_new = (m_dynamic_profile->points[0] - m_dynamic_profile->center).length();               // new smj
    m_dynamic_profile->smn_axis *= (smj_new / m_dynamic_profile->smj_axis);                             // update_smn
    m_dynamic_profile->smj_axis = smj_new;                                                              // update_smj

    osg::Vec2d vec_mn(m_dynamic_profile->points[0].y() - m_dynamic_profile->points[1].y(),
                      m_dynamic_profile->points[1].x() - m_dynamic_profile->points[0].x());
    vec_mn.normalize();                                                                                 // normalized minor axis vector
    m_dynamic_profile->points[2] = m_dynamic_profile->center - vec_mn * m_dynamic_profile->smn_axis;    // update p2
    m_dynamic_profile->points[3] = m_dynamic_profile->center + vec_mn * m_dynamic_profile->smn_axis;    // update p3
}

void ImageModeller::constrain_mouse_point() {

    switch(sp_constraints) {
    case spine_constraints::none:
        break;
    case spine_constraints::straight_planar:
        // mouse point should be on the projection of the 3D line
        break;
    case spine_constraints::planar:
        break;
    default:
        break;
    }
}

/*
// The spine points may be restricted by constraints. This function checks
// the mouse position along with the possible constraints and updates the
// spine point. If there are not any constraints, it does update the spine
// point with the mouse point. The current spine point is stored at the back
// of the m_vertices buffer.
void ImageModeller::constrain_spine_point_in_continuous_mode() {

    switch(sp_constraints) {
    case spine_constraints::none:
        m_vertices->back() = m_mouse;
        break;
    case spine_constraints::straight_planar:
        osg::Vec2d vec_mj = m_ellipse->points[1] - m_ellipse->points[0];
        osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
        vec_mn.normalize();
        osg::Vec2d mouse_vec = m_mouse - m_ellipse->center;
        osg::Vec2d proj_vec = vec_mn * (mouse_vec * vec_mn);
        osg::Vec2d proj_pt = m_ellipse->center + proj_vec;
        m_vertices->back() = proj_pt;
        break;
    }
}

// The spine point may be restricted by a constraint. This function checks the last clicked point or current mouse point
// along with the constraints and updates the spine point.
void ImageModeller::constrain_spine_point_in_piecewise_linear_mode() {

    // In piecewise linear mode, the last clicked point or the mouse point is at the back of the m_vertices.
    switch(sp_constraints) {
    case spine_constraints::none: // do nothing
        break;
    case spine_constraints::straight_planar: // project the point to the minor_axis

        // major axis vector of the base ellipse
        osg::Vec2d vec_mj = m_ellipse->points[1] - m_ellipse->points[0];
        // minor axis vector of the base ellipse
        osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
        vec_mn.normalize();
        osg::Vec2d vec = m_vertices->back() - m_ellipse->center;
        osg::Vec2d proj_vec = vec_mn * (vec * vec_mn);
        osg::Vec2d proj_pt = m_ellipse->center + proj_vec;
        m_vertices->back() = proj_pt;
        break;
    }
}
*/
