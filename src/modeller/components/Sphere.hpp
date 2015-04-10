#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "ComponentBase.hpp"

class Sphere : public ComponentBase {
public:
    Sphere(unsigned int component_id, const osg::Vec4& color);
};

#endif // SPHERE_HPP
