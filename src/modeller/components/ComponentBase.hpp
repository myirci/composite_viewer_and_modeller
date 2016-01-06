#ifndef COMPONENTBASE_HPP
#define COMPONENTBASE_HPP

#include <osg/Group>

class ComponentBase : public osg::Group {
public:

    ComponentBase(unsigned int component_id);
    unsigned int GetComponentId() const;
    virtual void DisplayVertexNormals(bool flag) = 0;
    virtual void Print() const;

protected:

    unsigned int m_component_id;                // id of the component
    osg::ref_ptr<osg::Vec3dArray> m_anchors;    // anchor points
    osg::ref_ptr<osg::Vec2dArray> m_frame_pts;  // local frame points

};

#endif // COMPONENTBASE_HPP
