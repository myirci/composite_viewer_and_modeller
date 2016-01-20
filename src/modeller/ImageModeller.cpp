#include "ImageModeller.hpp"
#include "optimization/CircleEstimator.hpp"
#include "optimization/ComponentSolver.hpp"
#include "optimization/ModelSolver.hpp"
#include "ProjectionParameters.hpp"
#include "UIHelper.hpp"
#include "../geometry/Circle3D.hpp"
#include "../geometry/Rectangle2D.hpp"
#include "../geometry/Ellipse2D.hpp"
#include "../geometry/Segment2D.hpp"
#include "../geometry/Plane3D.hpp"
#include "../geometry/Ray3D.hpp"
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
    m_gcyl_dmode(gcyl_drawing_mode::mode_0),
    m_left_click(false),
    m_right_click(false),
    m_bimg_exists(false),
    spd_mode(spine_drawing_mode::piecewise_linear),
    sp_constraints(spine_constraints::planar),
    sc_constraints(section_constraints::none),
    comp_type(component_type::generalized_cylinder),
    m_first_ellipse(new Ellipse2D),
    m_final_ellipse(new Ellipse2D),
    m_lsegment(new Segment2D),
    m_dsegment(new Segment2D),
    m_last_circle(new Circle3D),
    m_first_circle(new Circle3D),
    m_uihelper(nullptr),
    m_circle_estimator(new CircleEstimator),
    m_display_raycast(false),
    m_raycast(nullptr),
    m_scale_factor(0.35),
    m_num_right_click(0) {

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
    delete m_first_circle;
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

    m_gcyl_dmode = gcyl_drawing_mode::mode_0;
    m_left_click = false;
    m_right_click = false;
    m_vertices->clear();
    m_num_right_click = 0;
    m_segments.clear();
    m_uihelper->Reset();
}

void ImageModeller::DeleteLastSection() {

    if(m_gcyl.valid()) {
        m_gcyl->DeleteLastSection();
        m_uihelper->DeleteLastSpinePoint();
        m_segments.pop_back();
        // think a method to update the m_lsegmet, etc..
    }
}

void ImageModeller::OnLeftClick(double x, double y) {

    m_left_click = true;
    m_mouse.set(x,y);
    m_vertices->push_back(m_mouse);
    model_update();
}

void ImageModeller::OnRightClick(double x, double y) {

    if(m_gcyl_dmode == gcyl_drawing_mode::mode_3) {
        m_right_click = true;
        m_mouse.set(x,y);
        m_vertices->push_back(m_mouse);
        model_update();
    }
}

void ImageModeller::OnMouseMove(double x, double y) {

    if(m_gcyl_dmode != gcyl_drawing_mode::mode_0) {
        m_mouse.set(x,y);
        model_update();
    }
}

void ImageModeller::IncrementScaleFactor() {

    if(m_scale_factor < 1)
        m_scale_factor += 0.05;
    std::cout << "Current scale factor: " << m_scale_factor << std::endl;
}

void ImageModeller::DecrementScaleFactor() {

    if(m_scale_factor > 0.1)
        m_scale_factor -= 0.05;
    std::cout << "Current scale factor: " << m_scale_factor << std::endl;
}

// Execution of the modelling process is done within this function.
void ImageModeller::model_update() {

    if(comp_type == component_type::generalized_cylinder) model_generalized_cylinder();
    else                                                  std::cout << "Component type is unknown" << std::endl;
}

void ImageModeller::model_generalized_cylinder() {

    if(m_gcyl_dmode == gcyl_drawing_mode::mode_0) {

        if(m_left_click) {
            // first click
            m_left_click = false;
            m_uihelper->InitializeMajorAxisDrawing(m_mouse);
            m_gcyl_dmode = gcyl_drawing_mode::mode_1;
        }
    }
    else if(m_gcyl_dmode == gcyl_drawing_mode::mode_1) {

        if(m_left_click) {
            // second click: major axis has been determined.
            m_left_click = false;
            m_first_ellipse->update_major_axis(m_vertices->at(0), m_vertices->at(1));

            Segment2D seg;
            m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(Segment2D(m_vertices->at(0), m_vertices->at(1)), seg);
            m_segments.push_back(seg);

            m_uihelper->InitializeMinorAxisDrawing(m_mouse);
            m_gcyl_dmode = gcyl_drawing_mode::mode_2;
        }
        else {
            // Here we are executing the major axis drawing mode. In this mode we only update the end point of
            // the major axis with the mouse position. The user has not decided the major axis yet.
            m_uihelper->Updatep1(m_mouse);
        }
    }
    else if(m_gcyl_dmode == gcyl_drawing_mode::mode_2) {

        // calculate the possible ellpise based on the current mouse position
        calculate_ellipse(m_first_ellipse);
        m_uihelper->UpdateBaseEllipse(m_first_ellipse);

        if(m_left_click) {
            // third click: base ellipse (m_first_ellipse) has been determined.
            m_left_click = false;
            initialize_spine_drawing_mode(projection_type::perspective);
            m_uihelper->InitializeSpineDrawing(m_first_ellipse);
            m_gcyl_dmode = gcyl_drawing_mode::mode_3;

            /* to enable the test
             * uncomment the following line and comment out the m_gcyl_dmode = gcyl_drawing_mode::mode_3
            */
            // test_circle_estimation_from_major_axis();
        }
    }
    else if(m_gcyl_dmode == gcyl_drawing_mode::mode_3) {

        if(spd_mode == spine_drawing_mode::piecewise_linear) {

            update_dynamic_segment();

            if(m_left_click) {

                // left click is a new spine point
                m_left_click = false;
                *m_lsegment = *m_dsegment;
                m_uihelper->AddSpinePoint(m_dsegment->mid_point());
                // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_straight_cylinder();
                add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                // add_planar_section_to_the_generalized_cylinder_under_orthogonality_constraint();
            }

            else if(m_right_click) {

                m_right_click = false;
                *m_lsegment = *m_dsegment;

                if(sp_constraints == spine_constraints::constant_depth) {

                    // right click ends the modelling of the current generalized cylinder
                    m_uihelper->AddSpinePoint(m_dsegment->mid_point());
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_straight_cylinder();
                    add_planar_section_to_the_generalized_cylinder_under_perspective_projection();
                    // add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                    // add_planar_section_to_the_generalized_cylinder_under_orthogonality_constraint();
                    Reset2DDrawingInterface();
                    // project_generalized_cylinder(*m_gcyl);
                }
                else if(sp_constraints == spine_constraints::planar || sp_constraints == spine_constraints::straight_planar) {

                    ++m_num_right_click;
                    if(m_num_right_click == 1) {
                        m_uihelper->InitializeFinalEllipseDisplay(m_dsegment);
                        m_final_ellipse->update_major_axis(m_lsegment->pt1, m_lsegment->pt2);
                        Segment2D seg;
                        m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
                        m_segments.push_back(seg);
                    }
                    else if(m_num_right_click == 2) {
                        calculate_ellipse(m_final_ellipse);
                        m_uihelper->UpdateFinalEllipse(m_final_ellipse);
                        // compute_generalized_cylinder();
                        recompute_generalized_cylinder();
                        Reset2DDrawingInterface();
                    }
                }
            }
            else {
                if(sp_constraints == spine_constraints::constant_depth ) {
                    m_uihelper->UpdateSweepCurve(m_dsegment);
                    m_uihelper->SpinePointCandidate(m_dsegment->mid_point());
                }
                else if(sp_constraints == spine_constraints::planar || sp_constraints == spine_constraints::straight_planar) {
                    if(m_num_right_click == 0) {
                        m_uihelper->UpdateSweepCurve(m_dsegment);
                        m_uihelper->SpinePointCandidate(m_dsegment->mid_point());
                    }
                    else {
                        calculate_ellipse(m_final_ellipse);
                        m_uihelper->UpdateFinalEllipse(m_final_ellipse);
                    }
                }                
            }
        }
        else if(spd_mode == spine_drawing_mode::continuous) {
            std::cout << "Continuous mode is not implemented yet!" << std::endl;
        }
    }
}

void ImageModeller::estimate_first_circle_under_persective_projection() {

    // 1) Estimate the first 3D circle under perspective projection from the user drawn ellipse (m_first_ellipse)
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_first_ellipse, circles, m_fixed_depth);

    // 2) Select one of the two estimated circles based on how the user drew the ellipse
    if(count == 2)       *m_first_circle = circles[select_first_3d_circle(circles, m_first_ellipse)];
    else if (count == 1) *m_first_circle = circles[0];
    else                  std::cout << "ERROR: Perspective 3D circle estimation error " << std::endl;

    // 3) copy the first circle to the last circle
    *m_last_circle = *m_first_circle;
}

void ImageModeller::estimate_first_circle_under_orthographic_projection() {

    // 1) Estimate the first 3D circle
    estimate_3d_circle_under_orthographic_projection(m_first_ellipse, *m_first_circle);

    // 2) copy the first circle to the last circle
    *m_last_circle = *m_first_circle;
}

void ImageModeller::estimate_first_circle_under_orthogonality_constraint() {

    Circle3D circles[2];
    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*m_first_ellipse, elp_prj);
    m_circle_estimator->estimate_3d_circles_using_orthogonality_constraint(elp_prj, -m_pp->near, circles, false);

    *m_first_circle = circles[0]; // for now! Select the one that fits better!
    double ratio = (m_fixed_depth / m_first_circle->center[2]);

    m_first_circle->radius *= ratio;
    m_first_circle->center *= ratio;

    // copy the first circle to the last circle
    *m_last_circle = *m_first_circle;
}

void ImageModeller::compute_generalized_cylinder() {

    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_final_ellipse, circles, m_fixed_depth);

    Circle3D final_circle;
    // 2) Select one of the two estimated circles based on how the user drew the ellipse
    if(count == 2)       final_circle = circles[select_first_3d_circle(circles, m_final_ellipse)];
    else if (count == 1) final_circle = circles[0];
    else                 std::cout << "ERROR: Perspective 3D circle estimation error " << std::endl;

    if(sp_constraints == spine_constraints::planar) {

        // find the plane of the generalized cylinder's axis
        Eigen::Vector3d n = final_circle.normal.cross(m_first_circle->normal);
        n.normalize();
        Plane3D main_axis_plane(n, m_first_circle->center);

        // update the sections
        std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
        for(int i = 1; i < sections.size(); ++i) {
            double factor =  (- main_axis_plane.get_plane().w()) / (n.dot(sections[i].center));
            sections[i].center *= factor;
            sections[i].radius *= factor;

            sections[i].normal[2] = (-n[0]*sections[i].normal[0] - n[1]*sections[i].normal[1]) / n[2];
            sections[i].normal.normalize();
        }
        m_gcyl->Recalculate();

        // scale the last circle so that it's center lies on the computed plane
        double scale = m_first_circle->center.dot(n) / final_circle.center.dot(n);
        final_circle.radius *= scale;
        final_circle.center *= scale;
        if(final_circle.normal.dot(sections.back().normal) < 0)
            final_circle.normal *= -1;
    }
    else if(sp_constraints == spine_constraints::straight_planar) {

        m_component_solver->SolveDepth(*m_first_circle, final_circle);

        Eigen::Vector3d vec = m_first_circle->center - final_circle.center;
        Eigen::MatrixXd A = Eigen::MatrixXd(3, 2);
        Eigen::Vector2d sol;
        A(0,1) = vec[0];
        A(1,1) = vec[1];
        A(2,1) = vec[2];

        vec.normalize();
        if(vec.dot(m_first_circle->normal) < 0) vec *= -1;

        std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
        for(int i = 1; i < sections.size(); ++i) {

            A(0,0) = sections[i].center[0];
            A(1,0) = sections[i].center[1];
            A(2,0) = sections[i].center[2];
            sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(m_first_circle->center);
            std::cout << "sol[0]: " << sol[0] << std::endl;
            sections[i].center *= sol[0];
            sections[i].radius *= sol[0];
            // sections[i].normal = m_first_circle->normal;
            sections[i].normal = vec;
        }
        m_gcyl->Recalculate();

        final_circle.normal = vec;
    }

    m_gcyl->AddPlanarSection(final_circle);
    m_gcyl->Update();
}

void ImageModeller::recompute_generalized_cylinder() {

    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_final_ellipse, circles, m_fixed_depth);

    Circle3D final_circle;
    // 2) Select one of the two estimated circles based on how the user drew the ellipse
    if(count == 2)       final_circle = circles[select_first_3d_circle(circles, m_final_ellipse)];
    else                 final_circle = circles[0];

    if(sp_constraints == spine_constraints::planar) {

        // find the plane of the generalized cylinder's axis
        Eigen::Vector3d n = final_circle.normal.cross(m_first_circle->normal);
        n.normalize();
        Plane3D main_axis_plane(n, m_first_circle->center);

        // update the sections
        std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
        for(int i = 1; i < sections.size(); ++i) {
            // update the normal
            sections[i].normal[2] = (-n[0]*sections[i].normal[0] - n[1]*sections[i].normal[1]) / n[2];
            sections[i].normal.normalize();
            m_circle_estimator->estimate_unit_3d_circle_from_major_axis(m_segments[i], -m_pp->near, sections[i]);
            double factor =  (- main_axis_plane.get_plane().w()) / (n.dot(sections[i].center));
            sections[i].center *= factor;
            sections[i].radius *= factor;
        }
        m_gcyl->Recalculate();

        // scale the last circle so that it's center lies on the computed plane
        double scale = m_first_circle->center.dot(n) / final_circle.center.dot(n);
        final_circle.radius *= scale;
        final_circle.center *= scale;
        if(final_circle.normal.dot(sections.back().normal) < 0)
            final_circle.normal *= -1;
    }
    else if(sp_constraints == spine_constraints::straight_planar) {

        m_component_solver->SolveDepth(*m_first_circle, final_circle);

        Eigen::Vector3d vec = m_first_circle->center - final_circle.center;
        Eigen::MatrixXd A = Eigen::MatrixXd(3, 2);
        Eigen::Vector2d sol;
        A(0,1) = vec[0];
        A(1,1) = vec[1];
        A(2,1) = vec[2];

        vec.normalize();
        if(vec.dot(m_first_circle->normal) < 0) vec *= -1;

        std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
        for(int i = 1; i < sections.size(); ++i) {

            A(0,0) = sections[i].center[0];
            A(1,0) = sections[i].center[1];
            A(2,0) = sections[i].center[2];
            sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(m_first_circle->center);
            std::cout << "sol[0]: " << sol[0] << std::endl;
            sections[i].center *= sol[0];
            sections[i].radius *= sol[0];
            // sections[i].normal = m_first_circle->normal;
            sections[i].normal = vec;
        }
        m_gcyl->Recalculate();

        final_circle.normal = vec;
    }

    m_gcyl->AddPlanarSection(final_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection() {

    // 1) Estimate the normal of the circle
    m_tvec.normalize();
    m_last_circle->normal[0] = m_tvec.x();
    m_last_circle->normal[1] = m_tvec.y();
    m_last_circle->normal[2] = 0;
    const std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    if(m_last_circle->normal.dot(sections.back().normal) < 0) m_last_circle->normal *= -1;

    // 2) Set the depth of the last circle
    m_last_circle->center[2] = m_fixed_depth;

    // 3) Based on normal and depth estimation, estimate the 3D circle
    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
    m_segments.push_back(seg);
    m_circle_estimator->estimate_3d_circle_from_major_axis_when_circle_depth_is_fixed(seg, -m_pp->near, *m_last_circle);

    // 4) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection_straight_cylinder() {

    m_last_circle->normal = m_first_circle->normal;
    m_last_circle->center[2] = m_fixed_depth;
    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
    m_segments.push_back(seg);
    m_circle_estimator->estimate_3d_circle_from_major_axis_when_circle_depth_is_fixed(seg, -m_pp->near, *m_last_circle);
    m_component_solver->SolveDepth(*m_first_circle, *m_last_circle);
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_orthographic_projection() {

    // set the radius: proportional to the length of the semi-major axis
    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
    m_segments.push_back(seg);
    m_last_circle->radius = seg.half_length() * (m_fixed_depth / -m_pp->near);

    // set the center: should be scaled with respect to the fixed depth
    osg::Vec2d ctr = seg.mid_point();
    m_last_circle->center[0] = ctr.x();
    m_last_circle->center[1] = ctr.y();
    m_last_circle->center[2] = -m_pp->near;
    m_last_circle->center *= (m_fixed_depth / -m_pp->near);

    // update the circle normal
    m_tvec.normalize();
    m_last_circle->normal[0] = m_tvec.x();
    m_last_circle->normal[1] = m_tvec.y();
    m_last_circle->normal[2] = 0;
    std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    if(m_last_circle->normal.dot(sections.back().normal) < 0) m_last_circle->normal *= -1;

    if(sections.size() == 1)
        m_last_circle->normal = m_first_circle->normal;

    /*
    if(sections.size() == 1) {
        sections.back().normal = m_last_circle->normal;
        m_gcyl->Recalculate();
    }
    */


    // update the circle normal alternative but does not work when the mouse pointer does nor coincide
    // the center of the ellipse.
    // const std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    // m_last_circle->normal = (m_last_circle->center - sections.back().center).normalized();

    // add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_orthogonality_constraint() {

    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
    m_segments.push_back(seg);

    std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    if(sections.size() == 1) {

        Circle3D circles[2];
        Ellipse2D elp_prj;
        m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*m_first_ellipse, elp_prj);
        elp_prj.points[3] = seg.mid_point();
        m_circle_estimator->estimate_3d_circles_using_orthogonality_constraint(elp_prj, -m_pp->near, circles, true);

        *m_first_circle = circles[0]; // for now! Select the one that fits better!
        double ratio = (m_fixed_depth / m_first_circle->center[2]);

        m_first_circle->radius *= ratio;
        m_first_circle->center *= ratio;
        sections[0] = *m_first_circle;
        m_gcyl->Recalculate();
    }

    // set the radius: proportional to the length of the semi-major axis
    m_last_circle->radius = seg.half_length() * (m_fixed_depth / -m_pp->near);

    // set the center: should be scaled with respect to the fixed depth
    osg::Vec2d ctr = seg.mid_point();
    m_last_circle->center[0] = ctr.x();
    m_last_circle->center[1] = ctr.y();
    m_last_circle->center[2] = -m_pp->near;
    m_last_circle->center *= (m_fixed_depth / -m_pp->near);

    // update the circle normal
    m_tvec.normalize();
    m_last_circle->normal[0] = m_tvec.x();
    m_last_circle->normal[1] = m_tvec.y();
    m_last_circle->normal[2] = 0;
    if(m_last_circle->normal.dot(sections.back().normal) < 0) m_last_circle->normal *= -1;

    // if(sections.size() == 1) m_last_circle->normal = m_first_circle->normal;

    // add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

void ImageModeller::initialize_spine_drawing_mode(projection_type pt) {

    // modify the user clicked point with its projection on the minor-axis guide line
    m_vertices->at(2) = m_first_ellipse->points[2];

    // copy the base ellipse major axis into the last segment
    m_lsegment->pt1 = m_first_ellipse->points[0];
    m_lsegment->pt2 = m_first_ellipse->points[1];

    // initialize the generalized cylinder as a new node in the scene graph.
    if(m_gcyl.valid()) m_gcyl = nullptr;

    // estimate the first circle
    if(pt == projection_type::perspective) {
        estimate_first_circle_under_persective_projection();
    }
    else if(pt == projection_type::orthographic) {
        estimate_first_circle_under_orthographic_projection();
    }
    if(pt == projection_type::orthogonality_constraint) {
        estimate_first_circle_under_orthogonality_constraint();
    }

    m_gcyl = new GeneralizedCylinder(GenerateComponentId(), *m_first_circle, m_rtype);
    m_canvas->UsrAddSelectableNodeToDisplay(m_gcyl.get(), m_gcyl->GetComponentId());
}

void ImageModeller::calculate_ellipse(std::unique_ptr<Ellipse2D>& ellipse) {

    // calculate the end points of the minor_axis guide line
    osg::Vec2d vec_mj = ellipse->points[1] - ellipse->points[0];
    osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
    vec_mn.normalize();
    osg::Vec2d pt_0 = ellipse->center - vec_mn * ellipse->smj_axis;
    osg::Vec2d pt_1 = ellipse->center + vec_mn * ellipse->smj_axis;

    osg::Vec2d vec1 = m_mouse - pt_0;
    osg::Vec2d vec2 = pt_1 - pt_0;

    double ratio = (vec1 * vec2) / (vec2 * vec2) ;
    if(ratio >= 0.0 && ratio <= 1.0) {

        // find the projection point
        vec2.normalize();
        osg::Vec2d proj_point = (vec2 * 2 * ratio * ellipse->smj_axis) + pt_0;

        // update the m_first_ellipse
        ellipse->update_minor_axis(proj_point);
    }
}

void ImageModeller::update_dynamic_segment() {

    // copy the last segment into the dynamic segment
    *m_dsegment = *m_lsegment;

    // vector from last validated spine point to the current mouse position
    m_tvec = m_mouse - m_lsegment->mid_point();

    if(sp_constraints == spine_constraints::straight_planar) {
        // translate the dynamic segment to the current mouse point
        // under the constraint of straight line

        osg::Vec2d vec = m_first_ellipse->points[2] - m_first_ellipse->points[3];
        vec.normalize();
        m_tvec = vec * (vec * m_tvec);
        m_dsegment->translate(m_tvec);
    }
    else {

        // rotate the dynamic segment according to the bend of the spine curve
        osg::Vec2d dir = m_lsegment->direction();
        double angle = acos((m_tvec * dir) / m_tvec.length());

        if(dir.x()*m_tvec.y() - dir.y()*m_tvec.x() > 0) m_dsegment->rotate(angle - HALF_PI);
        else                                            m_dsegment->rotate(HALF_PI - angle);

        // translate the dynamic segment to the current mouse point
        m_dsegment->translate(m_tvec);
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

int ImageModeller::estimate_unit_3d_circles(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    return m_circle_estimator->estimate_unit_3d_circles(elp_prj, circles, m_pp.get());
}

void ImageModeller::estimate_3d_circle_under_orthographic_projection(std::unique_ptr<Ellipse2D>& ellipse, Circle3D& circle) {

    Ellipse2D elp_prj;
    m_pp->convert_ellipse_from_logical_device_coordinates_to_projected_coordinates(*ellipse, elp_prj);
    m_circle_estimator->estimate_3d_circles_under_orthographic_projection(elp_prj, circle, -m_pp->near);

    // perspective scaling:
    // -----------------------------------------------------------------------------------------------------
    double ratio = (m_fixed_depth / -m_pp->near);
    circle.radius *= ratio;
    circle.center *= ratio;
}

size_t ImageModeller::select_first_3d_circle(const Circle3D* const circles, std::unique_ptr<Ellipse2D>& ellipse) {

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
    Vector2D<double> vec2(ctr.x() - ellipse->points[2].x(), ctr.y() - ellipse->points[2].y());

    return (vec1.dot(vec2) < 0) ? 0 : 1;
}

size_t ImageModeller::select_parallel_circle(const Circle3D * const circles) {

    double a = (m_last_circle->normal).dot(circles[0].normal);
    double b = (m_last_circle->normal).dot(circles[1].normal);
    return std::abs(a) > std::abs(b) ? 0 : 1;
}

void ImageModeller::update_dynamic_segment_with_mirror_point(const osg::Vec2d& pt, bool first) {

    if(first) {
        m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - pt);
        m_dsegment->pt1 = pt;
    }
    else {
        m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - pt);
        m_dsegment->pt2 = pt;
    }
}

void ImageModeller::ray_cast_within_binary_image_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p0(static_cast<int>(m_dsegment->pt1.x()), static_cast<int>(m_dsegment->pt1.y()));
    Point2D<int> p1(static_cast<int>(m_dsegment->pt2.x()), static_cast<int>(m_dsegment->pt2.y()));
    m_canvas->UsrDeviceToLogical(p0);
    m_canvas->UsrDeviceToLogical(p1);

    // 2) calculate the casting direction vector,the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p1.x - p0.x) * 0.35), static_cast<int>((p1.y - p0.y) * 0.35));

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

    osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
    osg::Vec2d new_p1(p1_hit.x, p1_hit.y);

    if((hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {          // update p0, p1

        if((m_dsegment->pt1 - new_p0).length2() < (m_dsegment->pt2 - new_p1).length2()) {
            m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - new_p1);
            m_dsegment->pt2 = new_p1;
        }
        else {
            m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - new_p0);
            m_dsegment->pt1 = new_p0;
        }
    }
    else if((hit_result[0] || hit_result[1]) && !(hit_result[2] || hit_result[3])) {    // update p0, p1 is the mirror of p0
        m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - new_p0);
        m_dsegment->pt1 = new_p0;
    }
    else if(!(hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {    // update p1, p0 is the mirror of p1
        osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
        m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - new_p1);
    }
    else {
        return;
    }
}

void ImageModeller::ray_cast_within_gradient_image_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p1(static_cast<int>(m_dsegment->pt1.x()), static_cast<int>(m_dsegment->pt1.y()));
    m_canvas->UsrDeviceToLogical(p1);
    Point2D<int> p2(static_cast<int>(m_dsegment->pt2.x()), static_cast<int>(m_dsegment->pt2.y()));
    m_canvas->UsrDeviceToLogical(p2);

    // 2) calculate the vectors for ray casting in both directions.
    //    The length of the inward vector should not be larger than the half of the current segment length
    //    No limit for the outward vector length (only limited by the scale factor)
    Vector2D<int> outward_vec(static_cast<int>((p2.x - p1.x) * m_scale_factor),
                              static_cast<int>((p2.y - p1.y) * m_scale_factor));
    Vector2D<int> inward_vec = outward_vec;
    if(m_scale_factor > 0.5) {
        inward_vec *= (0.5 / m_scale_factor);
    }

    // 3) perform ray casts and display shot rays variables for ray casting
    OtbImageType::PixelType hit_val[4];
    Point2D<int> hit_idx[4];

    // 3.1) ray cast from p1-center direction
    Point2D<int> end = p1 + inward_vec;
    if(m_rect->intersect(p1, end))
        hit_val[0] = GradientImageRayCast(m_gimage, p1, end, hit_idx[0]);

    if(m_display_raycast) {
        m_raycast->at(0).x() = end.x;
        m_raycast->at(0).y() = end.y;
    }

    // 3.2) ray cast from p1-outside direction
    end = p1 - outward_vec;
    if(m_rect->intersect(p1, end))
        hit_val[1] = GradientImageRayCast(m_gimage, p1, end, hit_idx[1]);

    if(m_display_raycast) {
        m_raycast->at(1).x() = end.x;
        m_raycast->at(1).y() = end.y;
    }

    // 3.3) ray cast from p2-center direction
    end = p2 - inward_vec;
    if(m_rect->intersect(p2, end))
        hit_val[2] = GradientImageRayCast(m_gimage, p2, end, hit_idx[2]);

    if(m_display_raycast) {
        m_raycast->at(2).x() = end.x;
        m_raycast->at(2).y() = end.y;
    }

    // 3.4) ray cast from p2-outside direction
    end = p2 + outward_vec;
    if(m_rect->intersect(p2, end))
        hit_val[3] = GradientImageRayCast(m_gimage, p2, end, hit_idx[3]);

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
    OtbImageType::IndexType p1Idx, p2Idx;
    p1Idx[0] = p1.x; p1Idx[1] = p1.y;
    OtbImageType::PixelType p1val = m_gimage->GetPixel(p1Idx);
    p2Idx[0] = p2.x; p2Idx[1] = p2.y;
    OtbImageType::PixelType p2val = m_gimage->GetPixel(p2Idx);

    bool p1_updated = false;
    Point2D<int> p1_hit = p1;
    if(hit_val[0] > p1val && hit_val[1] > p1val) {              // both side hit
        if(hit_val[0] > hit_val[1]) p1_hit = hit_idx[0];        // select based on the value of the hit
        else                        p1_hit = hit_idx[1];
        p1_updated = true;
    }
    else if(hit_val[0] > p1val) {                               // only one side hit
        p1_hit = hit_idx[0];
        p1_updated = true;
    }
    else if(hit_val[1] > p1val) {                               // only one side hit
        p1_hit = hit_idx[1];
        p1_updated = true;
    }
    m_canvas->UsrDeviceToLogical(p1_hit);                       // convert hit point to logical coordinates

    bool p2_updated = false;
    Point2D<int> p2_hit = p2;
    if(hit_val[2] > p2val && hit_val[3] > p2val) {              // both side hit
        if(hit_val[2] > hit_val[3]) p2_hit = hit_idx[2];        // select based on the value of the hit
        else                        p2_hit = hit_idx[3];
        p2_updated = true;
    }
    else if(hit_val[2] > p2val) {                               // only one side hit
        p2_hit = hit_idx[2];
        p2_updated = true;
    }
    else if(hit_val[3] > p2val) {                               // only one side hit
        p2_hit = hit_idx[3];
        p2_updated = true;
    }
    m_canvas->UsrDeviceToLogical(p2_hit);                       // convert hit point to logical coordinates

    osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
    osg::Vec2d new_p2(p2_hit.x, p2_hit.y);
    osg::Vec2d mid_p = m_dsegment->mid_point();

    if(p1_updated && p2_updated) { // both points are updated

        m_dsegment->pt1 = new_p1;
        m_dsegment->pt2 = new_p2;
    }
    else if(p1_updated) { // update p1, p2 is the mirror of p1
        m_dsegment->pt1 = new_p1;
        m_dsegment->pt2 = (mid_p * 2) - new_p1;
    }
    else if(p2_updated) { // update p2, p1 is the mirror of p2
        m_dsegment->pt1 = (mid_p * 2) - new_p2;
        m_dsegment->pt2 = new_p2;
    }
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

void ImageModeller::test_circle_estimation_from_major_axis() {

    std::cout.precision(16);
    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);

    Circle3D circles[2];
    estimate_3d_circles_with_fixed_depth(m_first_ellipse, circles, m_fixed_depth);

    Circle3D c1;
    c1.normal = circles[0].normal;
    c1.center[2] = circles[0].center[2];
    Circle3D c2;
    c2.normal = circles[1].normal;
    c2.radius = circles[1].radius;

    m_circle_estimator->estimate_3d_circle_from_major_axis_when_circle_depth_is_fixed(seg, -m_pp->near, c1);
    m_circle_estimator->estimate_3d_circle_from_major_axis_when_circle_radius_is_fixed(seg, -m_pp->near, c2);

    std::cout << "----------------------" << std::endl;
    std::cout << circles[0] << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << c1 << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << circles[1] << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << c2 << std::endl;
    std::cout << "----------------------" << std::endl;


    /*
    Circle3D circles2[2];
    estimate_unit_3d_circles(m_first_ellipse, circles2);

    Circle3D c3;
    c3.normal = circles2[0].normal;
    c3.radius = 1;

    m_circle_estimator->estimate_unit_3d_circle_from_major_axis(seg, -m_pp->near, c3);

    std::cout << "----------------------" << std::endl;
    std::cout << circles2[0] << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << c3 << std::endl;
    */
}
