#ifndef GENERALIZEDCYLINDER_HPP
#define GENERALIZEDCYLINDER_HPP

#include "ComponentBase.hpp"
#include "GeneralizedCylinderGeometry.hpp"

class GeneralizedCylinder : public ComponentBase {
public:

    // constructors
    GeneralizedCylinder(unsigned int component_id, rendering_type rtype, unsigned int numpoints_per_section = 40, const osg::Vec4& color = osg::Vec4(1.0f,1.0f,0.0f,1.0f));
    GeneralizedCylinder(unsigned int component_id, const Circle3D& base_circle, rendering_type rtype, unsigned int numpoints_per_section = 40, const osg::Vec4& color = osg::Vec4(1.0f,1.0f,0.0f,1.0f));
    void AddPlanarSection(const Circle3D& circle);
    void DisplaySectionNormals(bool flag);
    void DisplayVertexNormals(bool flag) override;
    void ChangeRenderingType(rendering_type rtype);
    void SetSectionNormalsColor(const osg::Vec4& color);
    void SetVertexNormalsColor(const osg::Vec4& color);
    void Update();
    void Clear(bool update_flag);
    void Recalculate();
    void DeleteLastSection();
    const GeneralizedCylinderGeometry* const GetGeometry() const { return m_geometry.get(); }
    GeneralizedCylinderGeometry* GetGeometry()                   { return m_geometry.get(); }
protected:
    osg::ref_ptr<GeneralizedCylinderGeometry> m_geometry;
    osg::ref_ptr<osg::Switch> m_snormals;
    osg::ref_ptr<osg::Switch> m_vnormals;
    osg::Vec4 m_snormals_color;
    osg::Vec4 m_vnormals_color;
    bool m_display_section_normals;
    bool m_display_vertex_normals;
private:
    inline void add_to_section_normals(const Circle3D& circle, unsigned int scale = 1);
    inline void add_to_vertex_normals(size_t section_index);
};

#endif // GENERALIZEDCYLINDER_HPP
