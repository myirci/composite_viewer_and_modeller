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
    void ResetSweepCurve();
    void ResetSecondEllipse();
    void ResetRayCastDisplay();
    void InitializeMajorAxisDrawing(const osg::Vec2d& pt, bool first = true);
    void Updatep1(const osg::Vec2d& pt, bool first = true);
    void InitializeMinorAxisDrawing(const osg::Vec2d& pt, bool first = true);
    void UpdateEllipse(const std::unique_ptr<Ellipse2D>& elp, bool first = true);

    void InitializeAxisDrawing(const std::unique_ptr<Ellipse2D>& ellipse);
    void AxisPointCandidate(const osg::Vec2d& pt);
    void AddAxisPoint(const osg::Vec2d& pt);
    void UpdateSweepCurve(const std::unique_ptr<Ellipse2D>& ellipse);
    void UpdateSweepCurve(const std::unique_ptr<Segment2D>& segment);
    void DisplayLineStrip(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);
    void DisplayLineLoop(const std::vector<osg::Vec2d>& pts, const osg::Vec4& color);
    void DisplayRayCast(const osg::ref_ptr<osg::Vec2dArray>& pts);
    void DeleteLastAxisPoint();

private:

    sweep_curve_type m_sweep_type;
    osg::ref_ptr<osg::Vec2dArray>               m_sweep_ellipse_vertices; // for displaying the sweep ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_sweep_ellipse_arrays;   // draw arrays for sweep ellipse vertices

    osg::ref_ptr<osg::Vec2dArray>               m_sweepline_vertices;     // for displaying the sweepline
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_sweepline_arrays;       // draw arrays for m_sweepline_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_second_elp_vertices;    // for displaying the second ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_second_elp_arrays;      // draw arrays for second ellipse vertices

    osg::ref_ptr<osg::Vec2dArray>               m_first_elp_vertices;     // for displaying the first ellipse
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_first_elp_arrays;       // draw arrays for m_first_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_axis_vertices;          // for displaying the base ellipse
    osg::ref_ptr<osg::DrawArrays>               m_axis_array;             // draw arrays for m_first_elp_vertices

    osg::ref_ptr<osg::Vec2dArray>               m_ray_cast_vertices;      // for displaying the ray_casts
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_ray_cast_arrays;        // draw arrays for m_ray_cast_vertices

    std::vector<osg::ref_ptr<osg::Vec2dArray>>  m_proj_vertices;          // for displaying projections
    std::vector<osg::ref_ptr<osg::DrawArrays>>  m_proj_arrays;            // draw arrays for m_proj_vertices

    inline osg::Geometry* initialize_sweepline_display();
    inline osg::Geometry* initialize_sweep_ellipse_display();
    inline osg::Geometry* initialize_first_ellipse_display();
    inline osg::Geometry* initialize_second_ellipse_display();
    inline osg::Geometry* initialize_axis_display();
    inline osg::Geometry* initialize_ray_cast_display();

    osg::Geode* m_geode;
};

#endif // UIHELPER_HPP
