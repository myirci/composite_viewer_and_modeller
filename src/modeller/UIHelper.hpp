#ifndef UIHELPER_HPP
#define UIHELPER_HPP

#include <memory>
#include <vector>

#include "../geometry/Primitives.hpp"

#include <osg/Geode>
#include <osg/Array>

class Ellipse2D;
class Segment2D;

enum class sweep_curve_type : unsigned char {
    line,
    ellipse
};

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
    void UpdateSweepCurve(const std::unique_ptr<Ellipse2D>& ellipse);
    void UpdateSweepCurve(const std::unique_ptr<Segment2D>& segment);
    void UpdateBaseEllipse(const std::unique_ptr<Ellipse2D>& elp);
    void InitializeFinalEllipseDisplay(const std::unique_ptr<Segment2D>& segment);
    void UpdateFinalEllipse(const std::unique_ptr<Ellipse2D>& elp);
    void DisplayLineStrip(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);
    void DisplayLineLoop(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);
    void DisplayRayCast(const osg::ref_ptr<osg::Vec2dArray>& pts);
    void DeleteLastSpinePoint();

private:

    sweep_curve_type m_sweep_type;
    osg::ref_ptr<osg::Vec2dArray>               m_sweep_ellipse_vertices; // for displaying the sweep ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_sweep_ellipse_arrays;   // draw arrays for sweep ellipse vertices

    osg::ref_ptr<osg::Vec2dArray>               m_sweepline_vertices;     // for displaying the sweepline
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_sweepline_arrays;       // draw arrays for m_sweepline_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_final_ellipse_vertices;  // for displaying the last ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_final_ellipse_arrays;    // draw arrays for last ellipse vertices

    osg::ref_ptr<osg::Vec2dArray>               m_base_elp_vertices;      // for displaying the base ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_base_elp_arrays;        // draw arrays for m_base_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_spine_vertices;         // for displaying the base ellipse
    osg::ref_ptr<osg::DrawArrays>               m_spine_array;            // draw arrays for m_base_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_ray_cast_vertices;      // for displaying the ray_casts
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_ray_cast_arrays;        // draw arrays for m_ray_cast_vertices

    std::vector<osg::ref_ptr<osg::Vec2dArray>>  m_proj_vertices;          // for displaying projections
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_proj_arrays;            // draw arrays for m_proj_vertices

    inline osg::Geometry* initialize_sweepline_display();
    inline osg::Geometry* initialize_sweep_ellipse_display();
    inline osg::Geometry* initialize_base_ellipse_display();
    inline osg::Geometry* initialize_last_ellipse_display();
    inline osg::Geometry* initialize_spine_display();
    inline osg::Geometry* initialize_ray_cast_display();

    osg::Geode* m_geode;
};

#endif // UIHELPER_HPP
