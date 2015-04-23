#include "UIHelper.hpp"

#include "../geometry/Ellipse2D.hpp"

#include <memory>
#include <osg/Geometry>

UIHelper::UIHelper(osg::Geode* geode) :
    m_sweepline_vertices(nullptr),
    m_base_elp_vertices(nullptr),
    m_spine_vertices(nullptr),
    m_constraint_vertices(nullptr) {

    geode->addDrawable(initialize_sweepline_display());
    geode->addDrawable(initialize_base_ellipse_display());
    geode->addDrawable(initialize_spine_display());
    geode->addDrawable(initialize_constraint_display());

    // for smooth lines
    geode->getOrCreateStateSet()->setMode(GL_LINE_SMOOTH, osg::StateAttribute::ON);
}

osg::Geometry* UIHelper::initialize_sweepline_display() {

    // for dynamic ellipse display
    // 2 for sweepline                  Line     : m_sweepline_arrays[0], Color: Red
    // 10 for center                    Line Loop: m_sweepline_arrays[1], Color: Gray
    // 10 for p0                        Line Loop: m_sweepline_arrays[2], Color: Light Blue
    // 10 for p1                        Line Loop: m_sweepline_arrays[3], Color: Royal Blue

    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setUseDisplayList(false);
    geom->setUseVertexBufferObjects(true);
    geom->setDataVariance(osg::Object::DYNAMIC);

    m_sweepline_vertices = new osg::Vec2dArray(32);
    geom->setVertexArray(m_sweepline_vertices);

    m_sweepline_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINES));
    geom->addPrimitiveSet(m_sweepline_arrays.back());
    for(int i = 0; i < 3; ++i) {
        m_sweepline_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP));
        geom->addPrimitiveSet(m_sweepline_arrays.back());
    }

    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));                   // Red
    colors->push_back(osg::Vec4(0.2901961f, 0.2901961f, 0.2901961f, 1.0f)); // Gray
    colors->push_back(osg::Vec4(0.749019608f, 0.937254902f, 1.0f, 1.0f));   // Light Blue
    colors->push_back(osg::Vec4(0.2f, 0.2f, 1.0f, 1.0f));                   // Royal Blue
    geom->setColorArray(colors, osg::Array::BIND_PER_PRIMITIVE_SET);
    return geom.release();
}

osg::Geometry* UIHelper::initialize_base_ellipse_display() {

    // for base ellipse display
    // 2 for major axis             Lines: m_base_elp_arrays[0], Color: Red
    // 2 for minor axis guideline   Lines: m_base_elp_arrays[1], Color: Red
    // 40 for base ellipse          Line Loop: m_base_elp_arrays[2], Color: Cyan
    // 10 for p0                    Line Loop: m_base_elp_arrays[3], Color: Lawn Green
    // 10 for p1                    Line Loop: m_base_elp_arrays[4], Color: Lawn Green
    // 10 for center                Line Loop: m_base_elp_arrays[5], Color: Lawn Green
    // 10 for p2                    Line Loop: m_base_elp_arrays[6], Color: Lawn Green

    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setUseDisplayList(false);
    geom->setUseVertexBufferObjects(true);
    geom->setDataVariance(osg::Object::DYNAMIC);

    m_base_elp_vertices = new osg::Vec2dArray(84);
    geom->setVertexArray(m_base_elp_vertices);

    for(int i = 0; i < 2; ++i) {
        m_base_elp_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINES));
        geom->addPrimitiveSet(m_base_elp_arrays.back());
    }

    for(int i = 0; i < 5; ++i) {
        m_base_elp_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP));
        geom->addPrimitiveSet(m_base_elp_arrays.back());
    }

    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));                  // Red
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));                  // Red
    colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f));                  // Cyan
    colors->push_back(osg::Vec4(0.48627451f, 0.988235294f, 0.0f, 1.0f));   // Lawn Green
    colors->push_back(osg::Vec4(0.48627451f, 0.988235294f, 0.0f, 1.0f));   // Lawn Green
    colors->push_back(osg::Vec4(0.48627451f, 0.988235294f, 0.0f, 1.0f));   // Lawn Green
    colors->push_back(osg::Vec4(0.48627451f, 0.988235294f, 0.0f, 1.0f));   // Lawn Green

    geom->setColorArray(colors, osg::Array::BIND_PER_PRIMITIVE_SET);
    return geom.release();
}

osg::Geometry* UIHelper::initialize_spine_display() {

    // for spine display : extending with each spine point
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setUseDisplayList(false);
    geom->setUseVertexBufferObjects(true);
    geom->setDataVariance(osg::Object::DYNAMIC);

    m_spine_vertices = new osg::Vec2dArray;
    geom->setVertexArray(m_spine_vertices);
    m_spine_array = new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP);
    geom->addPrimitiveSet(m_spine_array);

    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(0.3373f, 0.102f, 0.5137f, 1.0f));
    geom->setColorArray(colors, osg::Array::BIND_PER_PRIMITIVE_SET);
    return geom.release();
}

osg::Geometry* UIHelper::initialize_constraint_display() {

    // for dynamic ellipse display
    // 2 for constraint line         : lines,     m_constraint_arrays[0], Color: Gray
    // 10 for mouse point projection : line_loop, m_constraint_arrays[1], Color: Gray

    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setUseDisplayList(false);
    geom->setUseVertexBufferObjects(true);
    geom->setDataVariance(osg::Object::DYNAMIC);

    m_constraint_vertices = new osg::Vec2dArray(12);
    geom->setVertexArray(m_constraint_vertices);

    m_constraint_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINES));
    geom->addPrimitiveSet(m_constraint_arrays.back());
    m_constraint_arrays.push_back(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP));
    geom->addPrimitiveSet(m_constraint_arrays.back());

    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f)); // Gray
    colors->push_back(osg::Vec4(0.65f, 0.65f, 0.65f, 1.0f)); // Gray
    geom->setColorArray(colors, osg::Array::BIND_PER_PRIMITIVE_SET);
    return geom.release();
}

void UIHelper::Reset() {

    for(auto it = m_sweepline_arrays.begin(); it != m_sweepline_arrays.end(); ++it)
        (*it)->setCount(0);
    for(auto it = m_base_elp_arrays.begin(); it != m_base_elp_arrays.end(); ++it)
        (*it)->setCount(0);
    for(auto it = m_constraint_arrays.begin(); it != m_constraint_arrays.end(); ++it)
        (*it)->setCount(0);

    m_spine_array->setCount(0);
    m_spine_vertices->clear();

    m_spine_vertices->dirty();
    m_sweepline_vertices->dirty();
    m_base_elp_vertices->dirty();
    m_constraint_vertices->dirty();
}

// pt = p0 (first click)
void UIHelper::InitializeMajorAxisDrawing(const osg::Vec2d& pt) {

    // store the first clicked point (p0) into the first and second slots
    m_base_elp_vertices->at(0) = pt;
    m_base_elp_vertices->at(1) = pt;

    // enable the display of the major axis
    m_base_elp_arrays[0]->setFirst(0);
    m_base_elp_arrays[0]->setCount(2);

    // display a small circle at the first click(p0)
    Ellipse2DLight tmp(3, 3, 0, pt);
    tmp.generate_points_on_the_ellipse(m_base_elp_vertices, 44, 10);
    m_base_elp_arrays[3]->setFirst(44);
    m_base_elp_arrays[3]->setCount(10);

    m_base_elp_vertices->dirty();
}

void UIHelper::Updatep1(const osg::Vec2d& pt) {

    m_base_elp_vertices->at(1) = pt;
    m_base_elp_vertices->dirty();
}

// pt = p1 (second click)
void UIHelper::InitializeMinorAxisDrawing(const osg::Vec2d& pt) {

    // store the second clicked point (p1) into the second slot
    m_base_elp_vertices->at(1) = pt;

    // display a small circle at the second clicked point (p1)
    Ellipse2DLight tmp(3, 3, 0, pt);
    tmp.generate_points_on_the_ellipse(m_base_elp_vertices, 54, 10);
    m_base_elp_arrays[4]->setFirst(54);
    m_base_elp_arrays[4]->setCount(10);

    // display a small cirle at the center of the ellipse
    tmp.center = (m_base_elp_vertices->at(0) + m_base_elp_vertices->at(1)) / 2.0;
    tmp.generate_points_on_the_ellipse(m_base_elp_vertices, 64, 10);
    m_base_elp_arrays[5]->setFirst(64);
    m_base_elp_arrays[5]->setCount(10);

    //display the minor axis guideline
    osg::Vec2d vec1 = m_base_elp_vertices->at(1) - tmp.center;
    osg::Vec2d vec2(-vec1.y(), vec1.x());
    m_base_elp_vertices->at(2) = tmp.center - vec2;
    m_base_elp_vertices->at(3) = tmp.center + vec2;
    m_base_elp_arrays[1]->setFirst(2);
    m_base_elp_arrays[1]->setCount(2);

    // display of the projection of mouse pointer onto the minor axis guideline
    // currently the projection is the same with the center of the ellipse
    tmp.generate_points_on_the_ellipse(m_base_elp_vertices, 74, 10);
    m_base_elp_arrays[6]->setFirst(74);
    m_base_elp_arrays[6]->setCount(10);

    m_base_elp_vertices->dirty();
}

void UIHelper::UpdateBaseEllipse(const std::unique_ptr<Ellipse2D>& elp) {

    elp->generate_points_on_the_ellipse(m_base_elp_vertices, 4, 40);
    Ellipse2DLight tmp(3, 3, 0, elp->points[2]);
    tmp.generate_points_on_the_ellipse(m_base_elp_vertices, 74, 10);
    m_base_elp_arrays[2]->setFirst(4);
    m_base_elp_arrays[2]->setCount(40);
    m_base_elp_vertices->dirty();
}

void UIHelper::InitializeSpineDrawing(const std::unique_ptr<Ellipse2D>& ellipse) {

    // initialize the spine drawing
    m_spine_vertices->push_back(ellipse->points[2]);
    m_spine_vertices->push_back(ellipse->points[2]);
    m_spine_array->setFirst(0);
    m_spine_array->setCount(2);
    m_spine_vertices->dirty();

    // initialize the major and minor axes display
    m_sweepline_vertices->at(0) = ellipse->points[0];
    m_sweepline_vertices->at(1) = ellipse->points[1];
    m_sweepline_vertices->at(2) = ellipse->points[2];
    m_sweepline_vertices->at(3) = ellipse->points[3];
    m_sweepline_arrays[0]->setFirst(0);
    m_sweepline_arrays[0]->setCount(2);

    // display the center
    Ellipse2DLight tmp(3, 3, 0, ellipse->center);
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 2, 10);
    m_sweepline_arrays[1]->setFirst(2);
    m_sweepline_arrays[1]->setCount(10);

    // display p0
    tmp.center = ellipse->points[0];
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 12, 10);
    m_sweepline_arrays[2]->setFirst(12);
    m_sweepline_arrays[2]->setCount(10);

    // display p1
    tmp.center = ellipse->points[1];
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 22, 10);
    m_sweepline_arrays[3]->setFirst(22);
    m_sweepline_arrays[3]->setCount(10);

    m_sweepline_vertices->dirty();
}

void UIHelper::UpdateSweepline(const std::unique_ptr<Ellipse2D>& ellipse) {

    m_sweepline_vertices->at(0) = ellipse->points[0];
    m_sweepline_vertices->at(1) = ellipse->points[1];

    Ellipse2DLight tmp(3, 3, 0, ellipse->center);
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 2, 10);
    tmp.center = ellipse->points[0];
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 12, 10);
    tmp.center = ellipse->points[1];
    tmp.generate_points_on_the_ellipse(m_sweepline_vertices, 22, 10);

    m_sweepline_vertices->dirty();
}

void UIHelper::SpinePointCandidate(const osg::Vec2d& pt) {

    m_spine_vertices->back() = pt;
    m_spine_vertices->dirty();
}

void UIHelper::AddSpinePoint(const osg::Vec2d& pt) {

    m_spine_vertices->back() = pt;
    m_spine_vertices->push_back(pt);
    m_spine_array->setFirst(0);
    m_spine_array->setCount(m_spine_vertices->size());
    m_spine_vertices->dirty();
}

void UIHelper::DisplayConstraintLine(const std::vector<osg::Vec2d>& pts) {

    m_constraint_vertices->at(0) = pts[0];
    m_constraint_vertices->at(1) = pts[1];
    m_constraint_arrays[0]->setCount(2);
    m_constraint_vertices->dirty();
}
























