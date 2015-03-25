#ifndef UIHELPER_HPP
#define UIHELPER_HPP

#include <memory>

#include "../geometry/Primitives.hpp"
#include "../geometry/Ellipse2D.hpp"

#include <osg/Geode>
#include <osg/Array>

class UIHelper {
public:

    UIHelper(osg::Geode* geode);
    void Reset();
    void InitializeMajorAxisDrawing(const osg::Vec2d& pt);
    void InitializeMinorAxisDrawing(const osg::Vec2d& pt);
    void InitializeSpineDrawing(const std::unique_ptr<Ellipse2D>& ellipse);
    void SpinePointCandidate(const osg::Vec2d& pt);
    void AddSpinePoint(const osg::Vec2d& pt);
    void Updatep1(const osg::Vec2d& pt);
    void UpdateDynamicProfile(const std::unique_ptr<Ellipse2D>& ellipse);
    void UpdateBaseEllipse(const std::unique_ptr<Ellipse2D>& elp);

private:

    osg::ref_ptr<osg::Vec2dArray>               m_dprofile_vertices;    // for displaying the dynamic_profile
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_dprofile_arrays;      // draw arrays for m_dprofile_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_base_elp_vertices;    // for displaying the base ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_base_elp_arrays;      // draw arrays for m_base_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_spine_vertices;       // for displaying the base ellipse
    osg::ref_ptr<osg::DrawArrays>               m_spine_array;          // draw arrays for m_base_elp_vertices

    inline osg::Geometry* initialize_dynamic_profile_display();
    inline osg::Geometry* initialize_base_ellipse_display();
    inline osg::Geometry* initialize_spine_display();
};

#endif // UIHELPER_HPP
