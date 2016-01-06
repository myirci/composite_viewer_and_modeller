#include "ComponentGeometryBase.hpp"
#include <iostream>

ComponentGeometryBase::ComponentGeometryBase(const osg::Vec4& color) {

    // general settings for dynamic modification
    setUseDisplayList(false);
    setUseVertexBufferObjects(true);
    setDataVariance(osg::Object::DYNAMIC);

    // set vertex array.
    m_vertices = new osg::Vec3Array;
    setVertexArray(m_vertices.get());

    // set normals array
    m_normals = new osg::Vec3Array;
    setNormalArray(m_normals.get(), osg::Array::BIND_PER_VERTEX);

    // set color
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(color);
    setColorArray(colors, osg::Array::BIND_OVERALL);
}

void ComponentGeometryBase::Update() {
    m_vertices->dirty();
}

void ComponentGeometryBase::SetColor(const osg::Vec4& color) {
    osg::Vec4Array* colors = static_cast<osg::Vec4Array*>(getColorArray());
    colors->back() = color;
    this->Update();
}

void ComponentGeometryBase::Print() const {
    std::cout << "num_vertices: " << m_vertices->size() << std::endl;
}


