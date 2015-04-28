#ifndef UIHELPER_HPP
#define UIHELPER_HPP

#include <memory>
#include <vector>

#include "../geometry/Primitives.hpp"

#include <osg/Geode>
#include <osg/Array>

class Ellipse2D;

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
    void UpdateSweepline(const std::unique_ptr<Ellipse2D>& ellipse);
    void UpdateBaseEllipse(const std::unique_ptr<Ellipse2D>& elp);
    void DisplayLineStrip(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);
    void DisplayLineLoop(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);

private:

    osg::ref_ptr<osg::Vec2dArray>               m_sweepline_vertices;   // for displaying the sweepline
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_sweepline_arrays;     // draw arrays for m_sweepline_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_base_elp_vertices;    // for displaying the base ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_base_elp_arrays;      // draw arrays for m_base_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_spine_vertices;       // for displaying the base ellipse
    osg::ref_ptr<osg::DrawArrays>               m_spine_array;          // draw arrays for m_base_elp_vertices

    std::vector<osg::ref_ptr<osg::Vec2dArray>>  m_proj_vertices;        // for displaying projections
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_proj_arrays;          // draw arrays for m_proj_vertices

    inline osg::Geometry* initialize_sweepline_display();
    inline osg::Geometry* initialize_base_ellipse_display();
    inline osg::Geometry* initialize_spine_display();

    osg::Geode* m_geode;
};

#endif // UIHELPER_HPP
