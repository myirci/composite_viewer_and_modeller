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
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_vindices;    // for vertical sections
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_hindices;    // for horizontal sections
    std::vector<osg::ref_ptr<osg::DrawElementsUInt>> m_tindices;    // for triangle strip and triangle fan rendering
public:
    GeneralizedCylinderGeometry(int num_points_per_section, const osg::Vec4& color);
    GeneralizedCylinderGeometry(const Circle3D& base_circle, int num_points_per_section, const osg::Vec4& color);
    void AddPlanarSection(const Circle3D& circle);
    void GetVertexNormals(size_t section_idx, osg::Vec3Array* vertices);
    void GetFirstCirclePoint(size_t section_idx, osg::Vec3d& pt);
    void ChangeRenderingType(rendering_type rtype);
    const std::vector<Circle3D>& GetSections() const;
    std::vector<Circle3D>& GetSections();
    unsigned int GetNumberOfSections() const;
    void Recalculate();
    void Clear(bool update_flag);
    void Print() const;
protected:
    void update_geometry_and_indices(size_t section_index, const Circle3D& circle);
};

#endif // GENERALIZED_CYLINDER_GEOMETRY_HPP
