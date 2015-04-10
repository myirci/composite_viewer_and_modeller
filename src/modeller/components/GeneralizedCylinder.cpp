
#include "GeneralizedCylinder.hpp"
#include "../../osg/OsgUtility.hpp"
#include "../../geometry/Circle3D.hpp"

#include <osg/Geode>
#include <osg/Switch>
#include <osg/MatrixTransform>

GeneralizedCylinder::GeneralizedCylinder(unsigned int component_id, unsigned int numpoints_per_section, const osg::Vec4& color) :
    ComponentBase(component_id),
    m_geometry(new GeneralizedCylinderGeometry(numpoints_per_section, color)),
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

GeneralizedCylinder::GeneralizedCylinder(unsigned int component_id, const Circle3D& base_circle, unsigned int numpoints_per_section, const osg::Vec4& color) :
    ComponentBase(component_id),
    m_geometry(new GeneralizedCylinderGeometry(base_circle, numpoints_per_section, color)),
    m_snormals_color(osg::Vec4(1,1,0,1)),
    m_vnormals_color(osg::Vec4(1,1,1,1)),
    m_snormals(new osg::Switch),
    m_vnormals(new osg::Switch),
    m_display_section_normals(false),
    m_display_vertex_normals(false) {

    // add the geometry
    osg::Geode* geode1 = new osg::Geode;
    geode1->addDrawable(m_geometry.get());
    addChild(geode1);

    // add the section normal
    add_to_section_normals(base_circle);
    addChild(m_snormals.get());

    // add vertex normals
    osg::Geode* geode2 = new osg::Geode;
    m_geometry->GetLastVertexNormals(geode2, m_vnormals_color);
    m_vnormals->addChild(geode2, m_display_vertex_normals);
    addChild(m_vnormals.get());
}

void GeneralizedCylinder::SetSectionNormalsColor(const osg::Vec4& color) {
    m_snormals_color = color;
}

void GeneralizedCylinder::SetVertexNormalsColor(const osg::Vec4& color) {
    m_vnormals_color = color;
}

void GeneralizedCylinder::AddPlanarSection(const Circle3D& circle) {

    // add the section
    m_geometry->AddPlanarSection(circle);

    // add the section normal
    osg::Vec3d ctr(circle.center[0], circle.center[1], circle.center[2]);
    osg::Vec3d nrm(circle.normal[0], circle.normal[1], circle.normal[2]);
    nrm *= 2;
    m_snormals->addChild(display_vector3d(ctr, nrm, m_snormals_color), m_display_section_normals);

    // add vertex normals
    osg::Geode* geode = new osg::Geode;
    m_geometry->GetLastVertexNormals(geode, m_vnormals_color);
    m_vnormals->addChild(geode, m_display_vertex_normals);
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

void GeneralizedCylinder::add_to_section_normals(const Circle3D& circle) {

    osg::Vec3d ctr(circle.center[0], circle.center[1], circle.center[2]);
    osg::Vec3d nrm(circle.normal[0], circle.normal[1], circle.normal[2]);
    m_snormals->addChild(display_vector3d(ctr, nrm, m_snormals_color), m_display_section_normals);
}
