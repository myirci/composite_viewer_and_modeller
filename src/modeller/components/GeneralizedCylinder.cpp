
#include "GeneralizedCylinder.hpp"
#include "../../osg/OsgUtility.hpp"
#include "../../geometry/Circle3D.hpp"

#include <osg/Geode>
#include <osg/Switch>
#include <osg/MatrixTransform>

GeneralizedCylinder::GeneralizedCylinder(unsigned int component_id, rendering_type rtype, unsigned int numpoints_per_section, const osg::Vec4& color) :
    ComponentBase(component_id),
    m_geometry(new GeneralizedCylinderGeometry(numpoints_per_section, color, rtype)),
    m_snormals_color(osg::Vec4(1,1,0,1)),
    m_vnormals_color(osg::Vec4(1,1,1,1)),
    m_snormals(new osg::Switch),
    m_vnormals(new osg::Switch),
    m_display_section_normals(false),
    m_display_vertex_normals(false) {

    // add the geometry
    osg::Geode* geode = new osg::Geode;
    geode->addDrawable(m_geometry.get());
    addChild(geode);
    addChild(m_snormals.get());
    addChild(m_vnormals.get());
}

GeneralizedCylinder::GeneralizedCylinder(unsigned int component_id, const Circle3D& base_circle, rendering_type rtype, unsigned int numpoints_per_section, const osg::Vec4& color) :
    ComponentBase(component_id),
    m_geometry(new GeneralizedCylinderGeometry(base_circle, numpoints_per_section, color, rtype)),
    m_snormals_color(osg::Vec4(1,1,0,1)),
    m_vnormals_color(osg::Vec4(1,1,1,1)),
    m_snormals(new osg::Switch),
    m_vnormals(new osg::Switch),
    m_display_section_normals(false),
    m_display_vertex_normals(false) {

    // add the geometry
    osg::Geode* geode = new osg::Geode;
    geode->addDrawable(m_geometry.get());
    addChild(geode);

    // add the section normal
    add_to_section_normals(base_circle, 2);
    addChild(m_snormals.get());

    // add vertex normals
    add_to_vertex_normals(0);
    addChild(m_vnormals.get());
}

void GeneralizedCylinder::SetSectionNormalsColor(const osg::Vec4& color) {
    m_snormals_color = color;
}

void GeneralizedCylinder::SetVertexNormalsColor(const osg::Vec4& color) {
    m_vnormals_color = color;
}

void GeneralizedCylinder::AddPlanarSection(const Circle3D& circle) {

    m_geometry->AddPlanarSection(circle);                         // add the section
    add_to_section_normals(circle, 2);                            // add the section normal
    add_to_vertex_normals(m_geometry->GetNumberOfSections() - 1); // add vertex normals
}

void GeneralizedCylinder::Recalculate() {

    // 1) Clear the existing geometry and section normals and vertex normals
    Clear(false);

    // 2) Recalculate the geometry
    m_geometry->Recalculate();

    // 3) Recalculate the vertex and section normals
    for(int i = 0; i < m_geometry->GetSections().size(); ++i) {
        add_to_section_normals(m_geometry->GetSections()[i], 2);
        add_to_vertex_normals(i);
    }
}

void GeneralizedCylinder::DeleteLastSection() {

    // need a far better implementation but for now all the generalized cylinder is
    // recalculated.
    m_geometry->GetSections().pop_back();
    Recalculate();
}

void GeneralizedCylinder::MakeTransparent() {

    unsigned int numch = getNumChildren();
    for(int i = 0; i < numch; ++i) {
        osg::Node* child = getChild(i);
        osg::StateSet* pStateSet = child->getOrCreateStateSet();
        pStateSet->setMode(GL_BLEND, osg::StateAttribute::ON);
        pStateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
    }
}

void GeneralizedCylinder::Clear(bool update_flag) {

    m_snormals->removeChildren(0, m_snormals->getNumChildren());
    m_vnormals->removeChildren(0, m_vnormals->getNumChildren());
    m_geometry->Clear(update_flag);
}

void GeneralizedCylinder::Update() {
    m_geometry->Update();
}

void GeneralizedCylinder::DisplaySectionNormals(bool flag) {

    for(int i = 0;  i < m_snormals->getNumChildren(); ++i)
        m_snormals->setValue(i, flag);
    m_display_section_normals = flag;
}

void GeneralizedCylinder::DisplayVertexNormals(bool flag) {

    for(int i = 0;  i < m_snormals->getNumChildren(); ++i)
        m_vnormals->setValue(i, flag);
    m_display_vertex_normals = flag;
}

void GeneralizedCylinder::ChangeRenderingType(rendering_type rtype) {

    m_geometry->ChangeRenderingType(rtype);
    m_geometry->Update();
}

void GeneralizedCylinder::add_to_section_normals(const Circle3D& circle, unsigned int scale) {

    osg::Vec3d ctr(circle.center[0], circle.center[1], circle.center[2]);
    osg::Vec3d nrm(circle.normal[0], circle.normal[1], circle.normal[2]);
    nrm *= scale;
    m_snormals->addChild(display_vector3d(ctr, nrm, m_snormals_color), m_display_section_normals);
}

void GeneralizedCylinder::add_to_vertex_normals(size_t section_index) {

    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    m_geometry->GetVertexNormals(section_index, vertices.get());
    m_vnormals->addChild(display_lines(vertices.get(), m_vnormals_color), m_display_vertex_normals);
}
