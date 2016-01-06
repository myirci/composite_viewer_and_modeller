#ifndef COMPONENT_GEOMETRY_BASE_HPP
#define COMPONENT_GEOMETRY_BASE_HPP

#include <osg/Geometry>

class ComponentGeometryBase : public osg::Geometry {
public:
    ComponentGeometryBase(const osg::Vec4& color);
    virtual void Update();
    virtual void SetColor(const osg::Vec4& color);
    virtual void Print() const;
protected:
    osg::ref_ptr<osg::Vec3Array> m_vertices;   // geometry
    osg::ref_ptr<osg::Vec3Array> m_normals;    // vertex normals for rendering
};

#endif // COMPONENT_GEOMETRY_BASE_HPP
