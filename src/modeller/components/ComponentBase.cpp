#include <iostream>
#include "ComponentBase.hpp"


ComponentBase::ComponentBase(unsigned int component_id) : m_component_id(component_id) {

    m_frame_pts = new osg::Vec2dArray;
    m_anchors = new osg::Vec3dArray;
}

unsigned int ComponentBase::GetComponentId() const {

    return m_component_id;
}

void ComponentBase::Print() const {

    std::cout << "component_id: " << m_component_id << std::endl;
    if(m_frame_pts.valid()) {
        std::cout << "local frame points:" << std::endl;
        for(auto it = m_frame_pts->begin(); it != m_frame_pts->end(); ++it)
            std::cout << it->x() << " " << it->y() << std::endl;
    }
}
