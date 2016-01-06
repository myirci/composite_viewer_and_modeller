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
    m_mode(drawing_mode::mode_0),
    m_left_click(false),
    m_right_click(false),
    m_bimg_exists(false),
    spd_mode(spine_drawing_mode::piecewise_linear),
    sp_constraints(spine_constraints::none),
    sc_constraints(section_constraints::constant),
    comp_type(component_type::generalized_cylinder),
    m_first_ellipse(new Ellipse2D),
    m_lsegment(new Segment2D),
    m_dsegment(new Segment2D),
    m_last_circle(new Circle3D),
    m_first_circle(new Circle3D),
    m_uihelper(nullptr),
    m_circle_estimator(new CircleEstimator),
    m_display_raycast(false),
    m_raycast(nullptr),
    m_scale_factor(0.35) {

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
                m_first_ellipse->update_major_axis(m_vertices->at(0), m_vertices->at(1));
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
                // third click: base ellipse (m_first_ellipse) has been determined.
                m_left_click = false;
                initialize_spine_drawing_mode(projection_type::orthographic);
                m_uihelper->InitializeSpineDrawing(m_first_ellipse);
                m_mode = drawing_mode::mode_3;
            }
        }
        else if(m_mode == drawing_mode::mode_3) {

            if(spd_mode == spine_drawing_mode::continuous) {

                /*
                if(m_left_click) {
                    m_left_click = false;
                    // End the modelling process for the current generalized cylinder with the 4th click. The last
                    // clicked point is accepted as the last sample point.
                    update_dynamic_profile();
                    // *m_last_profile = *m_dynamic_profile;
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_1();
                    add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
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
                        update_dynamic_profile();
                        *m_last_profile = *m_dynamic_profile;
                        add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                        // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_1();
                    }
                }
                */
            }
            else if(spd_mode == spine_drawing_mode::piecewise_linear) {

                update_dynamic_segment();

                if(m_right_click) {
                    // right click ends the modelling of the current component being modelled.
                    m_right_click = false;
                    m_uihelper->AddSpinePoint(m_mouse);
                    *m_lsegment = *m_dsegment;
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_1();
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_2();
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_3();
                    add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                    Reset2DDrawingInterface();
                    // m_component_solver->SolveGeneralizedCylinder(m_gcyl.get());
                    // m_solver->AddComponent(m_gcyl.get());
                    // project_generalized_cylinder(*m_gcyl);
                }
                else if(m_left_click) {

                    // left click is a new spine point
                    m_left_click = false;
                    m_uihelper->AddSpinePoint(m_mouse);
                    *m_lsegment = *m_dsegment;
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_1();
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_2();
                    // add_planar_section_to_the_generalized_cylinder_under_perspective_projection_3();
                    add_planar_section_to_the_generalized_cylinder_under_orthographic_projection();
                }
                else {
                    m_uihelper->UpdateSweepCurve(m_dsegment);
                    m_uihelper->SpinePointCandidate(m_mouse);
                }
            }
        }
    }
}

void ImageModeller::estimate_first_circle_under_persective_projection() {

    // 1) Estimate the first 3D circle under perspective projection from the user drawn ellipse (m_first_ellipse)
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_first_ellipse, circles, m_fixed_depth);

    // 2) Select one of the two estimated circles based on how the user drew the ellipse
    if(count == 2)       *m_first_circle = circles[select_first_3d_circle(circles)];
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

/*
 * The minor axis of the proceeding 2D elliptic profiles are estimated depending on the update of the
 * major axis of the profile. For instance, if the major axis is extended 10%, then the minor axis is
 * also extended 10%. The 3D section is estimated based on this 2D elliptic profile.
*/
void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection_1() {

    /*
    // 1) Estimate the 3D circles for the current profile and add it to the generalized cylinder.
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_last_profile, circles, m_fixed_depth);

    // 2) Select one of the two estimated circles based on the angle between the normals
    if(count == 2)       *m_last_circle = circles[select_parallel_circle(circles)];
    else if (count == 1) *m_last_circle = circles[0];
    else                  std::cout << "ERROR: Perspective 3D circle estimation error "  << std::endl;

    // 3) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
    */
}

/*
 * The minor axis of the proceeding 2D elliptic profiles are estimated depending on the update of the
 * major axis of the profile. For instance, if the major axis is extended 10%, then the minor axis is
 * also extended 10%. The 3D section is estimated based on this 2D elliptic profile. The normal of the
 * estimated circle is changed such that it becomes parallel to the main axis of the generalized cylinder
*/
void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection_2() {
/*
    // 1) Estimate the 3D circles for the current profile and add it to the generalized cylinder.
    Circle3D circles[2];
    int count = estimate_3d_circles_with_fixed_depth(m_last_profile, circles, m_fixed_depth);

    // 2) Select one of the two estimated circles based on the angle between the normals
    if(count == 2)       *m_last_circle = circles[select_parallel_circle(circles)];
    else if (count == 1) *m_last_circle = circles[0];
    else                  std::cout << "ERROR: Perspective 3D circle estimation error "  << std::endl;

    const std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    m_last_circle->normal = (sections.back().center - m_last_circle->center).normalized();
    if(m_first_circle->normal.dot(m_last_circle->normal))
        m_last_circle->normal *= -1;

    // 3) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
    */
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_perspective_projection_3() {
/*
    // 1) radius : proportional to the length of the semi-major axis
    m_last_circle->radius = m_first_circle->radius * (m_last_profile->smj_axis / m_first_ellipse->smj_axis);

    m_last_circle->normal = m_first_circle->normal;

    // 3) calculate the center

    // find the points on the projection plane
    osg::Vec2d p1, p2, ctr;
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->points[0], p1);
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->points[1], p2);
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->center, ctr);

    // calculate the tangent vectors
    osg::Vec3d p1vec(p1, -m_pp->near);
    osg::Vec3d p2vec(p2, -m_pp->near);
    osg::Vec3d direction_vec(p2vec.y() - p1vec.y(), p1vec.x() - p2vec.x(), 0);
    osg::Vec3d W1_ = direction_vec ^ p1vec;
    osg::Vec3d W2_ = direction_vec ^ p2vec;
    Eigen::Vector3d W1(W1_.x(), W1_.y(), W1_.z());
    Eigen::Vector3d W2(W2_.x(), W2_.y(), W2_.z());

    // calculate the offset planes
    Ray3D ray1(osg::Vec3d(0,0,0), p1vec);
    Ray3D ray2(osg::Vec3d(0,0,0), p2vec);
    Eigen::Vector3d N;
    m_last_circle->get_scaled_normal(N);
    Eigen::Matrix3d Nx3;
    Nx3 << 0, -N[2], N[1], N[2], 0, -N[0], -N[1], N[0], 0;
    Eigen::Matrix3d Nx3_sq = Nx3 * Nx3;
    Plane3D offset_pl1(W1[0], W1[1], W1[2], std::sqrt(-W1.transpose() * Nx3_sq * W1));
    if(!is_intersecting(ray2, offset_pl1))
        offset_pl1.get_plane().w() *= -1;

    Plane3D offset_pl2(W2[0], W2[1], W2[2], std::sqrt(-W2.transpose() * Nx3_sq * W2));
    if(!is_intersecting(ray1, offset_pl2))
        offset_pl2.get_plane().w() *= -1;

    // intersect two offset planes:
    osg::Vec3d start_pt;
    osg::Vec3d dir_vec;
    intersect_two_planes(offset_pl1, offset_pl2, start_pt, dir_vec);

    // calculate the center
    double coeff1[18];
    coeff1[0] = dir_vec.x() * dir_vec.x();
    coeff1[1] = 2 * dir_vec.x() * start_pt.x();
    coeff1[2] = start_pt.x() * start_pt.x() - N[1]*N[1] - N[2]*N[2];
    coeff1[3] = dir_vec.x() * dir_vec.y();
    coeff1[4] = dir_vec.x() * start_pt.y() + dir_vec.y() * start_pt.x();
    coeff1[5] = N[0]*N[1] + start_pt.x() * start_pt.y();
    coeff1[6] = dir_vec.x() * dir_vec.z();
    coeff1[7] = dir_vec.x() * start_pt.z() + dir_vec.z() * start_pt.x();
    coeff1[8] = N[0]*N[2] + start_pt.x() * start_pt.z();
    coeff1[9] = dir_vec.y() * dir_vec.y();
    coeff1[10] = 2 * dir_vec.y() * start_pt.y() ;
    coeff1[11] = start_pt.y() * start_pt.y() - N[0]*N[0] - N[2]*N[2];
    coeff1[12] = dir_vec.y() * dir_vec.z();
    coeff1[13] = dir_vec.y() * start_pt.z() + dir_vec.z() * start_pt.y();
    coeff1[14] = N[1]*N[2] + start_pt.y() * start_pt.z();
    coeff1[15] = dir_vec.z() * dir_vec.z();
    coeff1[16] = 2 * dir_vec.z() * start_pt.z();
    coeff1[17] = start_pt.z() * start_pt.z() - N[0]*N[0] - N[1]*N[1];

    double coeff2[17];
    coeff2[0] = coeff1[9]*coeff1[17] + coeff1[10]*coeff1[16] + coeff1[11]*coeff1[15] - 2*coeff1[12]*coeff1[14] - coeff1[13]*coeff1[13];
    coeff2[1] = coeff1[10]*coeff1[17] + coeff1[11]*coeff1[16] - 2*coeff1[13]*coeff1[14];
    coeff2[2] = coeff1[11]*coeff1[17] - coeff1[14]*coeff1[14];
    coeff2[3] = coeff1[6]*coeff1[14] + coeff1[7]*coeff1[13] + coeff1[8]*coeff1[12] - coeff1[3]*coeff1[17] - coeff1[4]*coeff1[16] - coeff1[5]*coeff1[15];
    coeff2[4] = coeff1[7]*coeff1[14] + coeff1[8]*coeff1[13] - coeff1[4]*coeff1[17] - coeff1[5]*coeff1[16];
    coeff2[5] = coeff1[8]*coeff1[14] - coeff1[5]*coeff1[17];
    coeff2[6] = coeff1[3]*coeff1[14] + coeff1[4]*coeff1[13] + coeff1[5]*coeff1[12] - coeff1[6]*coeff1[11] - coeff1[7]*coeff1[10] - coeff1[8]*coeff1[9];
    coeff2[7] = coeff1[4]*coeff1[14] + coeff1[5]*coeff1[13] - coeff1[7]*coeff1[11] - coeff1[8]*coeff1[10];
    coeff2[8] = coeff1[5]*coeff1[14] - coeff1[8]*coeff1[11];
    coeff2[9] = coeff1[0]*coeff1[17] + coeff1[1]*coeff1[16] + coeff1[2]*coeff1[15] - 2*coeff1[6]*coeff1[8] - coeff1[7]*coeff1[7];
    coeff2[10] = coeff1[1]*coeff1[17] + coeff1[2]*coeff1[16] - 2*coeff1[7]*coeff1[8];
    coeff2[11] = coeff1[2]*coeff1[17] - coeff1[8]*coeff1[8];
    coeff2[12] = coeff1[3]*coeff1[8] + coeff1[4]*coeff1[7] + coeff1[5]*coeff1[6] - coeff1[0]*coeff1[14] - coeff1[1]*coeff1[13] - coeff1[2]*coeff1[12];
    coeff2[13] = coeff1[4]*coeff1[8] + coeff1[5]*coeff1[7] - coeff1[1]*coeff1[14] - coeff1[2]*coeff1[13];
    coeff2[14] = coeff1[5]*coeff1[8] - coeff1[2]*coeff1[14];
    coeff2[15] = coeff1[0]*coeff1[11] + coeff1[1]*coeff1[10] + coeff1[2]*coeff1[9] - 2*coeff1[3]*coeff1[5] - coeff1[4]*coeff1[4];
    coeff2[16] = coeff1[1]*coeff1[11] + coeff1[2]*coeff1[10] - 2*coeff1[4]*coeff1[5];
    coeff2[17] = coeff1[2]*coeff1[11] - coeff1[5]*coeff1[5];

    Eigen::Matrix3d A;
    A << 1/m_pp->near, 0, 0, 0, 1/m_pp->near, 0, 0, 0, -1;

    Eigen::Vector3d vec(ctr.x(), ctr.y(), 1);
    double k0 = A.row(0).dot(vec);
    double k1 = A.row(1).dot(vec);
    double k2 = A.row(2).dot(vec);

    double coeff3[9];
    coeff3[0] = coeff2[0] *  k0*A(0,0) + coeff2[3] * (k1*A(0,0) + k0*A(1,0)) + coeff2[6] * (k2*A(0,0) + k0*A(2,0)) + coeff2[9] * k1*A(1,0) + coeff2[12] *(k2*A(1,0) + k1*A(2,0)) + coeff2[15] * k2*A(2,0);
    coeff3[1] = coeff2[1] *  k0*A(0,0) + coeff2[4] * (k1*A(0,0) + k0*A(1,0)) + coeff2[7] * (k2*A(0,0) + k0*A(2,0)) + coeff2[10] * k1*A(1,0) + coeff2[13] *(k2*A(1,0) + k1*A(2,0)) + coeff2[16] * k2*A(2,0);
    coeff3[2] = coeff2[2] *  k0*A(0,0) + coeff2[5] * (k1*A(0,0) + k0*A(1,0)) + coeff2[8] * (k2*A(0,0) + k0*A(2,0)) + coeff2[11] * k1*A(1,0) + coeff2[14] *(k2*A(1,0) + k1*A(2,0)) + coeff2[17] * k2*A(2,0);
    coeff3[3] = coeff2[0] *  k0*A(0,1) + coeff2[3] * (k1*A(0,1) + k0*A(1,1)) + coeff2[6] * (k2*A(0,1) + k0*A(2,1)) + coeff2[9] *  k1*A(1,1) + coeff2[12] *(k2*A(1,1) + k1*A(2,1)) + coeff2[15] * k2*A(2,1);
    coeff3[4] = coeff2[1] *  k0*A(0,1) + coeff2[4] * (k1*A(0,1) + k0*A(1,1)) + coeff2[7] * (k2*A(0,1) + k0*A(2,1)) + coeff2[10] * k1*A(1,1) + coeff2[13] *(k2*A(1,1) + k1*A(2,1)) + coeff2[16] * k2*A(2,1);
    coeff3[5] = coeff2[2] *  k0*A(0,1) + coeff2[5] * (k1*A(0,1) + k0*A(1,1)) + coeff2[8] * (k2*A(0,1) + k0*A(2,1)) + coeff2[11] * k1*A(1,1) + coeff2[14] *(k2*A(1,1) + k1*A(2,1)) + coeff2[17] * k2*A(2,1);
    coeff3[6] = coeff2[0] *  k0*A(0,2) + coeff2[3] * (k1*A(0,2) + k0*A(1,2)) + coeff2[6] * (k2*A(0,2) + k0*A(2,2)) + coeff2[9] *  k1*A(1,2) + coeff2[12] *(k2*A(1,2) + k1*A(2,2)) + coeff2[15] * k2*A(2,2);
    coeff3[7] = coeff2[1] *  k0*A(0,2) + coeff2[4] * (k1*A(0,2) + k0*A(1,2)) + coeff2[7] * (k2*A(0,2) + k0*A(2,2)) + coeff2[10] * k1*A(1,2) + coeff2[13] *(k2*A(1,2) + k1*A(2,2)) + coeff2[16] * k2*A(2,2);
    coeff3[8] = coeff2[2] *  k0*A(0,2) + coeff2[5] * (k1*A(0,2) + k0*A(1,2)) + coeff2[8] * (k2*A(0,2) + k0*A(2,2)) + coeff2[11] * k1*A(1,2) + coeff2[14] *(k2*A(1,2) + k1*A(2,2)) + coeff2[17] * k2*A(2,2);
    double t = (coeff3[0]*coeff3[5] - coeff3[2]*coeff3[3]) / (coeff3[1]*coeff3[3] - coeff3[0]*coeff3[4]);

    osg::Vec3d calculated_ctr = start_pt + dir_vec * t;
    m_last_circle->center[0] = calculated_ctr.x();
    m_last_circle->center[1] = calculated_ctr.y();
    m_last_circle->center[2] = calculated_ctr.z();

    std::cout << "Final Last Circle: " << std::endl;
    std::cout << *m_last_circle << std::endl;
    std::cout << "------------------------" << std::endl;

    // 4) Add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
    */
}

void ImageModeller::add_planar_section_to_the_generalized_cylinder_under_orthographic_projection() {

    // set the radius: proportional to the length of the semi-major axis
    Segment2D seg;
    m_pp->convert_segment_from_logical_device_coordinates_to_projected_coordinates(*m_lsegment, seg);
    m_last_circle->radius = seg.half_length() * (m_fixed_depth / -m_pp->near);

    // set the center: should be scaled with respect to the fixed depth
    osg::Vec2d ctr = seg.mid_point();
    m_last_circle->center[0] = ctr.x();
    m_last_circle->center[1] = ctr.y();
    m_last_circle->center[2] = -m_pp->near;
    m_last_circle->center *= (m_fixed_depth / -m_pp->near);

    // update the circle normal
    const std::vector<Circle3D>& sections = m_gcyl->GetGeometry()->GetSections();
    m_last_circle->normal = (m_last_circle->center - sections.back().center).normalized();

    // add estimated 3D circle to the generalized cylinder
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

/*
void ImageModeller::add_planar_section_to_the_generalized_cylinder_constrained() {

    // 1) set the radius : proportional to the length of the semi-major axis
    m_last_circle->radius = m_first_circle->radius * (m_last_profile->smj_axis / m_first_ellipse->smj_axis);

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
    if(smj_vec.x() * smn_vec.y() - smj_vec.y() * smn_vec.x() > 0) m_last_circle->normal[2] = cos(m_tilt_angle);
    else                                                          m_last_circle->normal[2] = -cos(m_tilt_angle);

    osg::ref_ptr<osg::Vec2dArray> projections = new osg::Vec2dArray(2);
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->points[0], projections->at(0));
    m_pp->convert_from_logical_device_coordinates_to_projected_coordinates(m_last_profile->points[1], projections->at(1));

    m_component_solver->SolveForSingleCircle(projections.get(), *m_last_circle);
    m_gcyl->AddPlanarSection(*m_last_circle);
    m_gcyl->Update();
}

*/

void ImageModeller::initialize_spine_drawing_mode(projection_type pt) {

    // modify the user clicked point with its projection on the minor-axis guide line
    m_vertices->at(2) = m_first_ellipse->points[2];

    // initialize the generalized cylinder as a new node in the scene graph.
    if(m_gcyl.valid()) m_gcyl = nullptr;

    // estimate the first circle
    if(pt == projection_type::perspective)
        estimate_first_circle_under_persective_projection();
    else if(pt == projection_type::orthographic)
        estimate_first_circle_under_orthographic_projection();

    m_gcyl = new GeneralizedCylinder(GenerateComponentId(), *m_first_circle, m_rtype);
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

    // copy the base ellipse major axis into the last segment
    m_lsegment->pt1 = m_first_ellipse->points[0];
    m_lsegment->pt2 = m_first_ellipse->points[1];
}

void ImageModeller::calculate_ellipse() {

    // calculate the end points of the minor_axis guide line
    osg::Vec2d vec_mj = m_first_ellipse->points[1] - m_first_ellipse->points[0];
    osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
    vec_mn.normalize();
    osg::Vec2d pt_0 = m_first_ellipse->center - vec_mn * m_first_ellipse->smj_axis;
    osg::Vec2d pt_1 = m_first_ellipse->center + vec_mn * m_first_ellipse->smj_axis;

    osg::Vec2d vec1 = m_mouse - pt_0;
    osg::Vec2d vec2 = pt_1 - pt_0;

    double ratio = (vec1 * vec2) / (vec2 * vec2) ;
    if(ratio >= 0.0 && ratio <= 1.0) {

        // find the projection point
        vec2.normalize();
        osg::Vec2d proj_point = (vec2 * 2 * ratio * m_first_ellipse->smj_axis) + pt_0;

        // update the m_first_ellipse
        m_first_ellipse->update_minor_axis(proj_point);

        // display the ellipse and a small circle on the projection point
        m_uihelper->UpdateBaseEllipse(m_first_ellipse);
    }
}

void ImageModeller::update_dynamic_segment() {

    // copy the last segment into the dynamic segment
    *m_dsegment = *m_lsegment;

    // vector from last validated spine point to the current mouse position
    osg::Vec2d tvec = m_mouse - m_lsegment->mid_point();

    if(sp_constraints == spine_constraints::straight_planar) {
        // translate the dynamic segment to the current mouse point
        m_dsegment->translate(tvec);
    }
    else {

        // rotate the dynamic segment according to the bend of the spine curve
        osg::Vec2d dir = m_lsegment->direction();
        double angle = acos((tvec * dir) / tvec.length());

        if(dir.x()*tvec.y() - dir.y()*tvec.x() > 0) m_dsegment->rotate(angle - HALF_PI);
        else                                        m_dsegment->rotate(HALF_PI - angle);

        // translate the dynamic segment to the current mouse point
        m_dsegment->translate(tvec);
    }

    if(m_bimg_exists)
        ray_cast_within_binary_image_for_profile_match();
    else
        ray_cast_within_gradient_image_for_profile_match();
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
    m_circle_estimator->estimate_3d_circles_under_orthographic_projection(elp_prj, circle, m_pp->near);

    // perspective scaling:
    // -----------------------------------------------------------------------------------------------------

    /* radius of the circle is equal to the length of the smj_vec_2d = smj_vec_3d = tilted_smn_vec_3d = ellipse.smj_axis on the
     * near plane (image plane). On the near plane, depth = 0 w.r.t computer graphics and depth = -near w.r.t mathematical notion.
     * Here we need to use the mathematical notion and take the depth = -near on the image plane. At the desired depth, the radius
     * will be calculated as below from the similar triangles.
     *
     * radius_at_desired_depth = (desired_depth / depth_at_near_plane) * radius_on_near_plane
     *                         = (desired_depth / -near) * length(smj_vec_2d)
     *
     */
    circle.radius *= (m_fixed_depth / -m_pp->near);

    /* center of the ellipse on the image plane is (ellipse.center.x(), ellipse.center.y(), -near)
     * however circle is moved to the desired depth. We need to multiply the center with the ratio (new_radius / old_radius) or
     * (new_depth / old_depth)
    */
    circle.center *= (m_fixed_depth / -m_pp->near);

    // othographic scaling:
    // -----------------------------------------------------------------------------------------------------
    /*
    // radius and x, y coordinates of the centre do not change. Only the z coordinate of the circle needs
    // to be updated with the desired depth.
    circle.center[2] = m_fixed_depth;
    */
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
    Vector2D<double> vec2(ctr.x() - m_first_ellipse->points[2].x(), ctr.y() - m_first_ellipse->points[2].y());

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
    Point2D<int> p1(static_cast<int>(m_dsegment->pt1.x()), static_cast<int>(m_dsegment->pt1.y()));
    Point2D<int> p2(static_cast<int>(m_dsegment->pt2.x()), static_cast<int>(m_dsegment->pt2.y()));
    m_canvas->UsrDeviceToLogical(p1);
    m_canvas->UsrDeviceToLogical(p2);

    // 2) calculate the casting direction vector,the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p2.x - p1.x) * m_scale_factor), static_cast<int>((p2.y - p1.y) * m_scale_factor));

    // 3) perform ray casts and display shot rays variables for ray casting
    bool hit_result[4] = { false, false, false, false };
    Point2D<int> hit[4];

    // 3.1) ray cast from p0-center direction
    Point2D<int> end = p1 + dir_vec;
    if(m_rect->intersect(p1, end))
        hit_result[0] = BinaryImageRayCast(m_bimage, p1, end, hit[0]);

    // 3.2) ray cast from p0-outside direction
    end = p1 - dir_vec;
    if(m_rect->intersect(p1, end))
        hit_result[1] = BinaryImageRayCast(m_bimage, p1, end, hit[1]);

    // 3.3) ray cast from p1-center direction
    end = p2 - dir_vec;
    if(m_rect->intersect(p2, end))
        hit_result[2] = BinaryImageRayCast(m_bimage, p2, end, hit[2]);

    // 3.4) ray cast from p1-outside direction
    end = p2 + dir_vec;
    if(m_rect->intersect(p2, end))
        hit_result[3] = BinaryImageRayCast(m_bimage, p2, end, hit[3]);

    // 4) analyze the result of the ray casts
    bool flagp1 = hit_result[0] || hit_result[1];
    bool flagp2 = hit_result[2] || hit_result[3];

    if(!flagp1 && !flagp2) { // no_hit
        return;
    }
    else if(flagp1 && !flagp2) { // only p1 hit

        Point2D<int> p1_hit;
        if(hit_result[0]) p1_hit = hit[0];
        else p1_hit = hit[1];
        m_canvas->UsrDeviceToLogical(p1_hit);
        update_dynamic_segment_with_mirror_point(osg::Vec2d(p1_hit.x, p1_hit.y), true);
    }
    else if(!flagp1 && flagp2) { // only p2 hit

        Point2D<int> p2_hit;
        if(hit_result[0]) p2_hit = hit[2];
        else p2_hit = hit[3];
        m_canvas->UsrDeviceToLogical(p2_hit);
        update_dynamic_segment_with_mirror_point(osg::Vec2d(p2_hit.x, p2_hit.y), false);
    }
    else { // both hit

        Point2D<int> p1_hit;
        if(hit_result[0]) p1_hit = hit[0];
        else p1_hit = hit[1];
        m_canvas->UsrDeviceToLogical(p1_hit);
        update_dynamic_segment_with_mirror_point(osg::Vec2d(p1_hit.x, p1_hit.y), true);

        /*
        if(hit_result[0] && hit_result[1]) { // hit both ways

            Point2D<int> p1_hit0 = hit[0];
            m_canvas->UsrDeviceToLogical(p1_hit0);
            osg::Vec2d vec1(p1_hit0.x, p1_hit0.y);

            Point2D<int> p1_hit1 = hit[1];
            m_canvas->UsrDeviceToLogical(p1_hit1);
            osg::Vec2d vec2(p1_hit1.x, p1_hit1.y);

            double l1 = (m_dsegment->pt1 - vec1).length2();
            double l2 = (m_dsegment->pt1 - vec2).length2();

            if(l1 < l2) { // apply same with the shrink

            }
            else if(l2 > l1) { // apply same with the enlarge

            }
            else {

            }
        }

        else if(!hit_result[0] && hit_result[1]) { // enlarge

            if(hit_result[2] && hit_result[3]) { // hit both ways
                std::cout << "Hit both sides: noth handled!" << std::endl;
            }
            else if(!hit_result[2] && hit_result[3]) { // enlarge
                // take the biggest enlargement
                Point2D<int> p1_hit = hit[1];
                m_canvas->UsrDeviceToLogical(p1_hit);
                osg::Vec2d p1h(p1_hit.x, p1_hit.y);

                Point2D<int> p2_hit = hit[3];
                m_canvas->UsrDeviceToLogical(p2_hit);
                osg::Vec2d p2h(p2_hit.x, p2_hit.y);

                double l1 = (p1h - m_dsegment->pt1).length2();
                double l2 = (p2h - m_dsegment->pt2).length2();
                if(l1 < l2) {
                    m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - p2h);
                    m_dsegment->pt2 = p2h;
                }
                else if(l1 > l2) {
                    m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - p1h);
                    m_dsegment->pt1 = p1h;
                }
                else {
                    m_dsegment->pt1 = osg::Vec2d(p1_hit.x, p1_hit.y);
                    m_dsegment->pt2 = osg::Vec2d(p2_hit.x, p2_hit.y);
                }
            }
            else { // shrink
                std::cout << "Opposite directions: noth handled!" << std::endl;
            }
        }
        else {  // shrink
            if(hit_result[2] && hit_result[3]) { // hit both ways
                std::cout << "Hit both sides: noth handled!" << std::endl;
            }
            else if(!hit_result[2] && hit_result[3]) { // enlarge
                std::cout << "Opposite directions: noth handled!" << std::endl;
            }
            else { // shrink
                // take the smallest enlargement
                Point2D<int> p1_hit = hit[0];
                m_canvas->UsrDeviceToLogical(p1_hit);
                osg::Vec2d p1h(p1_hit.x, p1_hit.y);

                Point2D<int> p2_hit = hit[2];
                m_canvas->UsrDeviceToLogical(p2_hit);
                osg::Vec2d p2h(p2_hit.x, p2_hit.y);

                double l1 = (p1h - m_dsegment->pt1).length2();
                double l2 = (p2h - m_dsegment->pt2).length2();
                if(l1 < l2) {
                    m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - p1h);
                    m_dsegment->pt1 = p1h;
                }
                else if (l2 < l1) {
                    osg::Vec2d p2h(p2_hit.x, p2_hit.y);
                    m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - p2h);
                    m_dsegment->pt2 = p2h;
                }
                else {
                    m_dsegment->pt1 = osg::Vec2d(p1_hit.x, p1_hit.y);
                    m_dsegment->pt2 = osg::Vec2d(p2_hit.x, p2_hit.y);
                }
            }
        }
        */
    }
}

void ImageModeller::ray_cast_within_gradient_image_for_profile_match() {

    // 1) transform the point coordinates to pixel coordinates
    Point2D<int> p0(static_cast<int>(m_dsegment->pt1.x()), static_cast<int>(m_dsegment->pt1.y()));
    Point2D<int> p1(static_cast<int>(m_dsegment->pt2.x()), static_cast<int>(m_dsegment->pt2.y()));
    m_canvas->UsrDeviceToLogical(p0);
    m_canvas->UsrDeviceToLogical(p1);

    // 2) calculate the casting direction vector, the change in major-axis length is limited by a factor
    Vector2D<int> dir_vec(static_cast<int>((p1.x - p0.x) * m_scale_factor), static_cast<int>((p1.y - p0.y) * m_scale_factor));

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
        m_dsegment->pt1 = osg::Vec2d(p0_hit.x, p0_hit.y);
        m_dsegment->pt2 = osg::Vec2d(p1_hit.x, p1_hit.y);
    }
    else if(p0_updated) { // update p0, p1 is the mirror of p0
        osg::Vec2d new_p0(p0_hit.x, p0_hit.y);
        m_dsegment->pt2 = m_dsegment->pt2 + (m_dsegment->pt1 - new_p0);
        m_dsegment->pt1 = new_p0;
    }
    else if(p1_updated) { // update p1, p0 is the mirror of p1
        osg::Vec2d new_p1(p1_hit.x, p1_hit.y);
        m_dsegment->pt1 = m_dsegment->pt1 + (m_dsegment->pt2 - new_p1);
        m_dsegment->pt2 = new_p1;
    }
    // not hits
    else {
        return;
    }
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
        osg::Vec2d vec_mj = m_first_ellipse->points[1] - m_first_ellipse->points[0];
        osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
        vec_mn.normalize();
        osg::Vec2d mouse_vec = m_mouse - m_first_ellipse->center;
        osg::Vec2d proj_vec = vec_mn * (mouse_vec * vec_mn);
        osg::Vec2d proj_pt = m_first_ellipse->center + proj_vec;
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
        osg::Vec2d vec_mj = m_first_ellipse->points[1] - m_first_ellipse->points[0];
        // minor axis vector of the base ellipse
        osg::Vec2d vec_mn(-vec_mj.y(), vec_mj.x());
        vec_mn.normalize();
        osg::Vec2d vec = m_vertices->back() - m_first_ellipse->center;
        osg::Vec2d proj_vec = vec_mn * (vec * vec_mn);
        osg::Vec2d proj_pt = m_first_ellipse->center + proj_vec;
        m_vertices->back() = proj_pt;
        break;
    }
}
*/
