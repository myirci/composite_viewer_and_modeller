#include <iostream>

#include "GeneralizedCylinderGeometry.hpp"
#include "../../osg/OsgUtility.hpp"
#include "../../geometry/Circle3D.hpp"

#include <osg/Geode>
#include <osg/LineWidth>

GeneralizedCylinderGeometry::GeneralizedCylinderGeometry(int num_points_per_section, const osg::Vec4& color) :
    ComponentGeometryBase(color),
    m_numpts(num_points_per_section),
    m_transforms(new osg::MatrixdArray),
    m_rtype(rendering_type::triangle_strip) { }

GeneralizedCylinderGeometry::GeneralizedCylinderGeometry(const Circle3D& base_circle, int num_points_per_section, const osg::Vec4& color) :
    ComponentGeometryBase(color),
    m_sections {base_circle},
    m_transforms(new osg::MatrixdArray),
    m_numpts(num_points_per_section),
    m_rtype(rendering_type::triangle_strip) {

    // update geometry (vertices and normals) and indices
    update_geometry_and_indices(base_circle);
}

const std::vector<Circle3D>& GeneralizedCylinderGeometry::GetSections() const {
    return m_sections;
}

void GeneralizedCylinderGeometry::ChangeRenderingType(rendering_type rtype) {

    // Step-1: Desired rendering type must be different than the current one
    if(m_rtype == rtype) return;

    // Step-2: Modify the geometry if necessary
    if(m_rtype == rendering_type::triangle_fan) {

        // delete the section centers
        auto it1 = m_vertices->begin();
        for(int i = 0; i < m_sections.size(); ++i) {
            std::advance(it1, m_numpts);
            it1 = m_vertices->erase(it1);
        }

        // clear and regenerate the normals
        m_normals->clear();
        for(int i = 0; i < m_sections.size(); ++i) {
            osg::Vec3 ctr(m_sections[i].center[0], m_sections[i].center[1], m_sections[i].center[2]);
            for(int j = 0; j < m_numpts; ++j) {
                osg::Vec3 nrm = m_vertices->at(j + i*m_numpts) - ctr;
                nrm.normalize();
                m_normals->push_back(nrm);
            }
        }
    }

    if(rtype == rendering_type::triangle_fan) {
        // insert the section centers
        auto it1 = m_vertices->begin();
        for(int i = 0; i < m_sections.size(); ++i) {
            i == 0 ? std::advance(it1, m_numpts) : std::advance(it1, m_numpts + 1);
            it1 = m_vertices->insert(it1, osg::Vec3(m_sections[i].center[0], m_sections[i].center[1], m_sections[i].center[2]));
        }

        //clear and regenerate the normals
        m_normals->clear();
        for(int i = 0; i < m_sections.size(); ++i) {
            osg::Vec3 nrm(m_sections[i].normal[0], m_sections[i].normal[1], m_sections[i].normal[2]);
            for(int j = 0; j < m_numpts+1; ++j)
                m_normals->push_back(nrm);
        }
    }

    // Step-3: Clear and regenerate the indices
    for(int i = 0; i < m_hindices.size(); ++i) m_hindices[i]->clear();
    for(int i = 0; i < m_vindices.size(); ++i) m_vindices[i]->clear();
    for(int i = 0; i < m_tindices.size(); ++i) m_tindices[i]->clear();
    m_hindices.clear();
    m_vindices.clear();
    m_tindices.clear();

    switch (rtype) {
    case rendering_type::planar_sections: {
        for(int i = 0; i < m_sections.size(); ++i) {
            m_hindices.push_back(new osg::DrawElementsUInt(GL_LINE_LOOP));
            addPrimitiveSet(m_hindices[i].get());
            for(int j = 0; j < m_numpts; ++j)
                m_hindices[i]->push_back(j + i * m_numpts);
        }
        break;
    }

    case rendering_type::planar_and_vertical_sections: {
        for(int i = 0; i < m_sections.size(); ++i) {
            m_hindices.push_back(new osg::DrawElementsUInt(GL_LINE_LOOP));
            addPrimitiveSet(m_hindices[i].get());

            if(i == 0) {
                for(int j = 0; j < m_numpts; ++j) {
                    m_vindices.push_back(new osg::DrawElementsUInt(GL_LINE_STRIP));
                    addPrimitiveSet((m_vindices.back()).get());
                }
            }

            for(int j = 0; j < m_numpts; ++j) {
                m_hindices[i]->push_back(j + i * m_numpts);
                m_vindices[j]->push_back(j + i * m_numpts);
            }
        }
        break;
    }

    case rendering_type::triangle_fan: {
        for(int i = 0; i < m_sections.size(); ++i) {
            m_tindices.push_back(new osg::DrawElementsUInt(GL_TRIANGLE_FAN));
            m_tindices[i]->push_back(i * m_numpts + i + m_numpts);
            for(int j = 0; j < m_numpts; ++j)
                m_tindices[i]->push_back(j + (m_numpts+1) * i);
            m_tindices[i]->push_back((m_numpts+1)*i);
            addPrimitiveSet(m_tindices[i].get());
        }
        break;
    }

    case rendering_type::triangle_strip: {

        for(int i = 1; i < m_sections.size(); ++i) {
            m_tindices.push_back(new osg::DrawElementsUInt(GL_TRIANGLE_STRIP));
            int start_index = (i-1)*m_numpts;
            for(int j = start_index; j < start_index + m_numpts; ++j) {
                m_tindices[i-1]->push_back(j);
                m_tindices[i-1]->push_back(j + m_numpts);
            }
            m_tindices[i-1]->push_back(start_index);
            m_tindices[i-1]->push_back(start_index + m_numpts);
            addPrimitiveSet(m_tindices[i-1].get());
        }
        break;
    }

    default: {
        std::cout << "ERROR: Unknown render type" << std::endl;
        break;
    }
    }

    m_rtype = rtype;
}

void GeneralizedCylinderGeometry::update_geometry_and_indices(const Circle3D& circle) {

    unsigned int section_index = m_sections.size() - 1;

    // update_geometry
    if(m_rtype == rendering_type::triangle_fan) {
        circle.generate_data(m_vertices, m_numpts);
        m_vertices->push_back(osg::Vec3(circle.center[0], circle.center[1], circle.center[2]));
        for(size_t i = 0; i < m_numpts + 1; ++i)
            m_normals->push_back(osg::Vec3(circle.normal[0], circle.normal[1], circle.normal[2]));
    }
    else {
        if(section_index == 0)
            circle.generate_data(m_vertices, m_normals, m_numpts);
        else {
            osg::ref_ptr<osg::Vec3Array> new_vertices = new osg::Vec3Array;
            osg::ref_ptr<osg::Vec3Array> new_normals = new osg::Vec3Array;
            circle.generate_data(new_vertices, new_normals, m_numpts);

            rearrange_vertices_and_normals(new_vertices, new_normals);
            for(auto it = new_vertices->begin(); it != new_vertices->end(); ++it)
                m_vertices->push_back(*it);
            for(auto it = new_normals->begin(); it != new_normals->end(); ++it)
                m_normals->push_back(*it);
        }
    }

    // update indices
    switch(m_rtype) {

        case rendering_type::planar_sections: {

            m_hindices.push_back(new osg::DrawElementsUInt(GL_LINE_LOOP));
            addPrimitiveSet(m_hindices[section_index].get());
            for(int i = 0; i < m_numpts; ++i)
                m_hindices[section_index]->push_back(i + section_index * m_numpts);
            break;
        }

        case rendering_type::planar_and_vertical_sections: {

            m_hindices.push_back(new osg::DrawElementsUInt(GL_LINE_LOOP));
            addPrimitiveSet(m_hindices[section_index].get());
            if(section_index == 0) {
                for(int i = 0; i < m_numpts; ++i) {
                    m_vindices.push_back(new osg::DrawElementsUInt(GL_LINE_STRIP));
                    addPrimitiveSet((m_vindices.back()).get());
                }
            }

            for(int i = 0; i < m_numpts; ++i) {
                m_hindices[section_index]->push_back(i + section_index * m_numpts);
                m_vindices[i]->push_back(i + section_index * m_numpts);
            }

            break;
        }

        case rendering_type::triangle_fan: {

            m_tindices.push_back(new osg::DrawElementsUInt(GL_TRIANGLE_FAN));
            m_tindices[section_index]->push_back(section_index * m_numpts + section_index + m_numpts);
            for(int i = 0; i < m_numpts; ++i)
                m_tindices[section_index]->push_back(i + (m_numpts+1) * section_index);
            m_tindices[section_index]->push_back((m_numpts+1)*section_index);
            addPrimitiveSet(m_tindices[section_index].get());
            break;
        }

        case rendering_type::triangle_strip: {

            if(section_index == 0) break;
            m_tindices.push_back(new osg::DrawElementsUInt(GL_TRIANGLE_STRIP));
            int start_index = (section_index-1)*m_numpts;

            for(int i = start_index; i < start_index + m_numpts; ++i) {
                m_tindices[section_index-1]->push_back(i);
                m_tindices[section_index-1]->push_back(i + m_numpts);
            }
            m_tindices[section_index-1]->push_back(start_index);
            m_tindices[section_index-1]->push_back(start_index + m_numpts);
            addPrimitiveSet(m_tindices[section_index-1].get());
            break;
        }

        default: {

            std::cout << "ERROR: Unknown render type" << std::endl;
            break;
        }
    }
}

void GeneralizedCylinderGeometry::AddPlanarSection(const Circle3D& section) {

    // step-1: push the given 3D circle into the sections vector
    m_sections.push_back(section);

    // step-2: update geometry (vertices and normals) and indices
    update_geometry_and_indices(section);

    // step-3: calculate the transformation matrix between the last two circles
    if(m_sections.size() > 1) {
        osg::Matrixd mat;
        calculate_transformation_matrix(m_sections[m_sections.size()-2], m_sections[m_sections.size()-1], mat);
        m_transforms->push_back(std::move(mat));
    }
}

void GeneralizedCylinderGeometry::GetLastVertexNormals(osg::Geode* geode, const osg::Vec4& color) {

    if(m_vertices->empty() || m_normals->empty()) {
        std::cout << "ERROR: No vertices and/or normals " << std::endl;
        return;
    }

    size_t numpts;
    m_rtype == rendering_type::triangle_fan ? numpts = m_numpts + 1 : numpts = m_numpts;

    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(color);

    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    size_t size = m_sections.size() - 1;
    size_t idx = 0;
    for(int i = 0; i < numpts; ++i) {
        idx = size * numpts + i;
        vertices->push_back(m_vertices->at(idx));
        vertices->push_back(m_vertices->at(idx) + m_normals->at(idx));
    }

    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setVertexArray(vertices.get());
    geom->setColorArray(colors.get());
    geom->setColorBinding(osg::Geometry::BIND_OVERALL);
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertices->size()));

    geode->addDrawable(geom.release());
    osg::LineWidth* linewidth = new osg::LineWidth();
    linewidth->setWidth(0.3);
    geode->getOrCreateStateSet()->setAttributeAndModes(linewidth, osg::StateAttribute::ON);
    geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

}

void GeneralizedCylinderGeometry::rearrange_vertices_and_normals(osg::ref_ptr<osg::Vec3Array>& new_vertices, osg::ref_ptr<osg::Vec3Array>& new_normals) const {
#ifdef DEBUG
    assert(m_numpts != 0);
    assert(m_vertices->size() == m_normals->size());
#endif // DEBUG

    double min = std::numeric_limits<double>::max();
    double dist = 0;
    osg::Vec3Array::iterator rot1, rot2;
    osg::Vec3Array::iterator it1, it2;

    osg::Vec3d pt(m_vertices->at(m_vertices->size() - m_numpts));

    for(it1 = new_vertices->begin(), it2 = new_normals->begin(); it1 != new_vertices->end(); ++it1, ++it2) {
        dist = squared_distance(pt, *it1);
        if(dist < min) {
            min = dist;
            rot1 = it1;
            rot2 = it2;
        }
    }

   std::rotate(new_vertices->begin(), rot1, new_vertices->end());
   std::rotate(new_normals->begin(), rot2, new_normals->end());
}

void GeneralizedCylinderGeometry::Print() const {
    std::cout << "----------------------------" << std::endl;
    std::cout << "Generalized cylinder" << std::endl;
    ComponentGeometryBase::Print();
    std::cout << "----------------------------" << std::endl;
}

