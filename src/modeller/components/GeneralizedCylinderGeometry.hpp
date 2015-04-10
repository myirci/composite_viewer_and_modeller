#ifndef GENERALIZED_CYLINDER_GEOMETRY_HPP
#define GENERALIZED_CYLINDER_GEOMETRY_HPP

#include "ComponentGeometryBase.hpp"
#include <memory>
#include <vector>

class Circle3D;

enum class rendering_type : unsigned char {
    planar_sections,
    planar_and_vertical_sections,
    triangle_fan,
    triangle_strip
};

class GeneralizedCylinderGeometry : public ComponentGeometryBase {
protected:
    int m_numpts;                                                   // number of points for each planar section
    rendering_type m_rtype;                                         // rendering type for the generalized cylinder
    std::vector<Circle3D> m_sections;                               // planar sections
    osg::ref_ptr<osg::MatrixdArray> m_transforms;                   // transformation matrix
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_vindices;    // for vertical sections
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_hindices;    // for horizontal sections
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_tindices;    // for triangle strip and triangle fan rendering
public:
    GeneralizedCylinderGeometry(int num_points_per_section, const osg::Vec4& color);
    GeneralizedCylinderGeometry(const Circle3D& base_circle, int num_points_per_section, const osg::Vec4& color);

    void AddPlanarSection(const Circle3D& circle);
    void GetLastVertexNormals(osg::Geode* geode, const osg::Vec4& color);
    void ChangeRenderingType(rendering_type rtype);
    const std::vector<Circle3D>& GetSections() const;
    void Print() const;
protected:
    void update_geometry_and_indices(const Circle3D& circle);
    void rearrange_vertices_and_normals(osg::ref_ptr<osg::Vec3Array>& new_vertices, osg::ref_ptr<osg::Vec3Array>& new_normals) const;
};

#endif // GENERALIZED_CYLINDER_GEOMETRY_HPP
