#ifndef RECTANGLE2D_HPP
#define RECTANGLE2D_HPP

#include "Primitives.hpp"

enum class intersection_zone : unsigned char {
    top_left,
    top,
    top_right,
    left,
    right,
    bottom_left,
    bottom,
    bottom_right
};

enum class segment_type : unsigned char {
    left,
    right,
    top,
    bottom
};

struct Rectangle2D {
public:
    Point2D<int> corner_points[4];
    Rectangle2D(int x_0 = 0, int y_0 = 0, int x_1 = 0, int y_1 = 0);
    bool intersect(Point2D<int>& pt0, Point2D<int>& pt1) const;
private:
    intersection_zone query_zone(const Point2D<int>& pt) const;
    bool perfom_intersection(const Point2D<int>& pt_in, Point2D<int>& pt_out) const;
    bool intersect_with_segment(const Point2D<int>& pt_in, Point2D<int>& pt_out, segment_type stype) const;
    inline bool is_point_inside(const Point2D<int>& pt) const;
};

#endif // RECTANGLE2D_HPP
