#include "ImageModeller.hpp"
#include "optimization/CircleEstimator.hpp"
#include "optimization/ModelSolver.hpp"
#include "PersProjParam.hpp"
#include "UIHelper.hpp"
#include "../geometry/Circle3D.hpp"
#include "../geometry/Line2D.hpp"
#include "../geometry/Rectangle2D.hpp"
#include "../geometry/Ellipse2D.hpp"
#include "../osg/OsgWxGLCanvas.hpp"
#include "../osg/OsgUtility.hpp"
#include "../utility/Utility.hpp"
#include "../image/algorithms/Algorithms.hpp"
#include <osgDB/WriteFile>
#include <otbImageFileReader.h>

ImageModeller::ImageModeller(const std::string& fpath, const std::shared_ptr<PersProjParam>& ppp, OsgWxGLCanvas* canvas) :
    m_ppp(ppp),
    m_canvas(canvas),
    m_gcyl(nullptr),
    m_vertices(nullptr),
    m_rect(nullptr),
    m_rtype(rendering_type::triangle_strip),
    m_solver(new ModelSolver()),
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
    m_uihelper(nullptr),
    m_constraint_line(nullptr) {

    std::ifstream infile(fpath);
    if(!infile.good()) {
        std::cout << "INFO: No associated binary image file: " << std::endl;
    }
    else {
        otb::ImageFileReader<ImageType>::Pointer reader = otb::ImageFileReader<ImageType>::New();
        reader->SetFileName(fpath);
        reader->Update();
        m_image = reader->GetOutput();
        m_bimg_exists = true;
        // binary image is ready for processing
        ImageType::SizeType size = m_image->GetLargestPossibleRegion().GetSize();
        m_rect = std::unique_ptr<Rectangle2D>(new Rectangle2D(0, 0, size[0] - 1, size[1] - 1));
    }
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
                    add_planar_section_to_the_generalized_cylinder();

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
                        add_planar_section_to_the_generalized_cylinder();
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
                    add_planar_section_to_the_generalized_cylinder();
                    // update_component_local_frame();
                    m_solver->AddComponent(m_gcyl.get());
                    Reset2DDrawingInterface();
                }
                else if(m_left_click) {
                    // left click is a new spine point
                    m_left_click = false;
                    m_uihelper->AddSpinePoint(m_mouse);
                    *m_last_profile = *m_dynamic_profile;
                    add_planar_section_to_the_generalized_cylinder();
                }
                else {
                    m_uihelper->SpinePointCandidate(m_mouse);
                    m_uihelper->UpdateDynamicProfile(m_dynamic_profile);
                }
            }
        }
    }
}

void ImageModeller::estimate_first_circle() {



}

void ImageModeller::initialize_spine_drawing_mode() {

    // modify the user clicked point with its projection on the minor-axis guide line
    m_vertices->at(2) = m_ellipse->points[2];

    // estimate the 3D cicrles.
    Circle3D circles[2];
    // estimate_unit_3d_circles(m_ellipse, circles);
    estimate_3d_circles_with_fixed_depth(m_ellipse, circles, -(m_ppp->near + m_ppp->far)/2.0);
    // estimate_3d_circles_under_orthographic_projection(m_ellipse, circles[0]);
    // estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(m_ellipse, circles[0], -(m_ppp->near + m_ppp->far)/2.0);
    // estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(m_ellipse, circles[0], -(m_ppp->near + m_ppp->far)/2.0);

    // initialize the generalized cylinder as a new node in the scene graph.
    if(m_gcyl.valid()) m_gcyl = nullptr;
    *m_last_circle = circles[select_first_3d_circle(circles)];
    // *m_last_circle = circles[0];
    m_gcyl = new GeneralizedCylinder(GenerateComponentId(), *m_last_circle);
    m_canvas->UsrAddSelectableNodeToDisplay(m_gcyl.get(), m_gcyl->GetComponentId());

    // initialize the constraint line if spine contraint is straight_planar
    if(sp_constraints == spine_constraints::straight_planar) {

        // let's define the constraint line in the opengl screen coordinate which is also
        double a = m_last_circle->center[1] * m_last_circle->normal[2] - m_last_circle->center[2] * m_last_circle->normal[1];
        double b = m_last_circle->center[2] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[2];
        double c = (m_last_circle->center[1] * m_last_circle->normal[0] - m_last_circle->center[0] * m_last_circle->normal[1]) * m_ppp->height / (2*tan(m_ppp->fovy/2.0)) - a*(m_ppp->width/2.0) - b*(m_ppp->height/2.0);
        m_constraint_line = std::unique_ptr<Line2D>(new Line2D(a,b,c));

        std::cout << "Constraint line: " << std::endl;
        m_constraint_line->print();

        std::cout << "Perspective projection parameters " << std::endl;
        std::cout << *m_ppp << std::endl;

        // display the constraint line
        std::vector<osg::Vec2d> intersection_pts;
        double var = m_constraint_line->get_x_at_y(0);
        if(var >= 0 && var < m_ppp->width)
            intersection_pts.push_back(osg::Vec2d(var, 0));

        var = m_constraint_line->get_x_at_y(m_ppp->height - 1);
        if(var >= 0 && var < m_ppp->width)
            intersection_pts.push_back(osg::Vec2d(var, m_ppp->height - 1));

        var = m_constraint_line->get_y_at_x(0);
        if(var >= 0 && var < m_ppp->height)
            intersection_pts.push_back(osg::Vec2d(0, var));

        var = m_constraint_line->get_y_at_x(m_ppp->width - 1);
        if(var >= 0 && var < m_ppp->height)
            intersection_pts.push_back(osg::Vec2d(m_ppp->width - 1, var));

        if(intersection_pts.size() == 2)
            m_uihelper->DisplayConstraintLine(intersection_pts);
        else
            std::cout << "ERROR: Number of intersections must be two!" << std::endl;

        std::cout << "intersection points: \n";
        for(auto it = intersection_pts.begin(); it != intersection_pts.end(); ++it)
            std::cout << it->x() << " " << it->y() << std::endl;
    }

    // Copy the base ellipse into the m_last_profile.
    *m_last_profile = *m_ellipse;
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder() {

    // estimate the 3D circle for the current profile and add it to the generatlized cylinder.
    Circle3D circles[2];
    // estimate_unit_3d_circles(m_last_profile, circles);
    estimate_3d_circles_with_fixed_depth(m_last_profile, circles, -(m_ppp->near + m_ppp->far)/2.0);
    // estimate_3d_circles_under_orthographic_projection(m_last_profile, circles[0]);
    // estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(m_last_profile, circles[0], -(m_ppp->near + m_ppp->far)/2.0);
    // estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(m_last_profile, circles[0], -(m_ppp->near + m_ppp->far)/2.0);

    *m_last_circle = circles[select_parallel_circle(circles)];
    // *m_last_circle = circles[0];
    // *m_last_circle = circles[select_first_3d_circle(circles)];
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
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
        m_dynamic_profile->translate(m_mouse - m_dynamic_profile->points[2]);
        // m_dynamic_profile->translate(m_mouse - m_dynamic_profile->center);
    }

    if(m_bimg_exists) ray_cast_for_profile_match();
}

void ImageModeller::estimate_3d_circles_with_fixed_radius(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_radius) {

    CircleEstimator estimator;
    ellipse->calculate_algebraic_equation_in_projected_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0));
    estimator.estimate_3d_circles_with_fixed_radius(*ellipse, circles, m_ppp.get(), desired_radius);
}

void ImageModeller::estimate_3d_circles_with_fixed_depth(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles, double desired_depth) {

    CircleEstimator estimator;
    ellipse->calculate_algebraic_equation_in_projected_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0));
    estimator.estimate_3d_circles_with_fixed_depth(*ellipse, circles, m_ppp.get(), desired_depth);
}

void ImageModeller::estimate_unit_3d_circles(std::unique_ptr<Ellipse2D>& ellipse, Circle3D* circles) {

    CircleEstimator estimator;
    ellipse->calculate_algebraic_equation_in_projected_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0));
    estimator.estimate_unit_3d_circles(*ellipse, circles, m_ppp.get());
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection(std::unique_ptr<Ellipse2D>& ellipse, Circle3D& circle) {

    CircleEstimator estimator;
    Ellipse2D elp_prj;
    convert_ellipse_from_logical_device_coordinates_to_projected_point_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0), *ellipse, elp_prj);
    estimator.estimate_3d_circles_under_orthographic_projection(elp_prj, circle, m_ppp->near);
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(std::unique_ptr<Ellipse2D>& ellipse, Circle3D &circle, double desired_depth) {

    CircleEstimator estimator;
    Ellipse2D elp_prj;
    convert_ellipse_from_logical_device_coordinates_to_projected_point_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0), *ellipse, elp_prj);
    estimator.estimate_3d_circles_under_orthographic_projection_and_scale_orthographically(elp_prj, circle, m_ppp->near, desired_depth);
}

void ImageModeller::estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(std::unique_ptr<Ellipse2D>& ellipse, Circle3D &circle, double desired_depth) {

    CircleEstimator estimator;
    Ellipse2D elp_prj;
    convert_ellipse_from_logical_device_coordinates_to_projected_point_coordinates(m_ppp->width, m_ppp->height, m_ppp->near, deg2rad(m_ppp->fovy/2.0), *ellipse, elp_prj);
    estimator.estimate_3d_circles_under_orthographic_projection_and_scale_perspectively(elp_prj, circle, m_ppp->near, desired_depth);
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
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    return std::abs(a) > std::abs(b) ? 0 : 1;
}

size_t ImageModeller::select_planar_circle(const Circle3D* const circles) { }

void ImageModeller::ray_cast_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p0(static_cast<int>(m_dynamic_profile->points[0].x()), static_cast<int>(m_dynamic_profile->points[0].y()));
    Point2D<int> p1(static_cast<int>(m_dynamic_profile->points[1].x()), static_cast<int>(m_dynamic_profile->points[1].y()));
    m_canvas->UsrTransformCoordinates(p0);
    m_canvas->UsrTransformCoordinates(p1);

    // 2) calculate the ray_cast direction vector,the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p1.x - p0.x) * 0.35), static_cast<int>((p1.y - p0.y) * 0.35));

    // 3) perform ray casts and display shot rays variables for ray casting
    bool hit_result[4] = { false, false, false, false };
    Point2D<int> hit[4];

    // 3.1) ray cast from p0-center direction
    Point2D<int> end = p0 + dir_vec;
    if(m_rect->intersect(p0, end)) hit_result[0] = ray_cast(m_image, p0, end, hit[0]);

    // 3.2) ray cast from p0-outside direction
    end = p0 - dir_vec;
    if(m_rect->intersect(p0, end)) hit_result[1] = ray_cast(m_image, p0, end, hit[1]);

    // 3.3) ray cast from p1-center direction
    end = p1 - dir_vec;
    if(m_rect->intersect(p1, end)) hit_result[2] = ray_cast(m_image, p1, end, hit[2]);

    // 3.4) ray cast from p1-outside direction
    end = p1 + dir_vec;
    if(m_rect->intersect(p1, end)) hit_result[3] = ray_cast(m_image, p1, end, hit[3]);

    // 5) analyze the result of the ray casts
    Point2D<int> p0_hit = p0;
    if(hit_result[0] && hit_result[1])  p0_hit = p0;
    else if(hit_result[0])              p0_hit = hit[0];
    else if(hit_result[1])              p0_hit = hit[1];
    m_canvas->UsrTransformCoordinates(p0_hit);

    Point2D<int> p1_hit;
    if(hit_result[2] && hit_result[3])  p1_hit = p1;
    else if(hit_result[2])              p1_hit = hit[2];
    else if(hit_result[3])              p1_hit = hit[3];
    m_canvas->UsrTransformCoordinates(p1_hit);

    // if p0 and p1 hits
    if((hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {
        // update p0, p1 and center, semi-major and semi-minor
        m_dynamic_profile->points[0] = osg::Vec2d(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = osg::Vec2d(p1_hit.x, p1_hit.y);
        m_dynamic_profile->center = (m_dynamic_profile->points[0] + m_dynamic_profile->points[1]) / 2.0;

        // osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
        // m_dynamic_profile->points[1] = m_dynamic_profile->points[1] + (m_dynamic_profile->points[0] - new_p0);
        // m_dynamic_profile->points[0] = new_p0;
    }
    // only p0_hits
    else if((hit_result[0] || hit_result[1]) && !(hit_result[2] || hit_result[3])) {
        // the second point is the mirror of the hit_point
        osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
        m_dynamic_profile->points[1] = m_dynamic_profile->points[1] + (m_dynamic_profile->points[0] - new_p0);
        m_dynamic_profile->points[0] = new_p0;
    }
    // only p1_hits
    else if(!(hit_result[0] || hit_result[1]) && (hit_result[2] || hit_result[3])) {
        // the second point is the mirror of the hit_point
        osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
        m_dynamic_profile->points[0] = m_dynamic_profile->points[0] + (m_dynamic_profile->points[1] - new_p1);
        m_dynamic_profile->points[1] = new_p1;
    }
    else {
        return;
    }

    osg::Vec2d vec_mn(p0_hit.y - p1_hit.y, p1_hit.x - p0_hit.x);
    vec_mn.normalize();
    double smj_new = (m_dynamic_profile->points[0] - m_dynamic_profile->center).length();
    m_dynamic_profile->smn_axis *= (smj_new / m_dynamic_profile->smj_axis);
    m_dynamic_profile->points[2] = m_dynamic_profile->center - vec_mn * m_dynamic_profile->smn_axis;
    m_dynamic_profile->points[3] = m_dynamic_profile->center + vec_mn * m_dynamic_profile->smn_axis;
    m_dynamic_profile->smj_axis = smj_new;
}


void ImageModeller::constrain_mouse_point() {

    switch(sp_constraints) {
    case spine_constraints::none:
        break;
    case spine_constraints::straight_planar:
    {
        // mouse point should be on the projection of the 3D line
        break;
    }
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
