#include "Rectangle2D.hpp"

Rectangle2D::Rectangle2D(int x_0, int y_0, int x_1, int y_1) {

    if(y_1 - y_0 <= 0 || x_1 - x_0 <= 0) {
        std::cout << "ERROR: Ractangle2D neagtive/zero width/height" << std::endl;
        return;
    }
    corner_points[0] = Point2D<int>(x_0, y_0);
    corner_points[1] = Point2D<int>(x_1, y_0);
    corner_points[2] = Point2D<int>(x_1, y_1);
    corner_points[3] = Point2D<int>(x_0, y_1);
}

bool Rectangle2D::intersect(Point2D<int>& pt0, Point2D<int>& pt1) const {

    bool pt0_in = is_point_inside(pt0);
    bool pt1_in = is_point_inside(pt1);

    if(pt0_in) {
        if(pt1_in) return true;                             // pt0 & pt1: both in
        else return perfom_intersection(pt0, pt1);          // pt0: in; pt1: out
    }
    else {
        if(pt1_in) return perfom_intersection(pt1, pt0);    // pt0: out; pt1: in
    }
    // pt0 & pt1: both out
    return false;
}

// including the borders of the rectangle
bool Rectangle2D::is_point_inside(const Point2D<int>& pt) const {
    return pt.x >= corner_points[0].x && pt.x <= corner_points[2].x && pt.y >= corner_points[0].y && pt.y <= corner_points[2].y;
}

// precondition: the point resides outside the rectangle, including the borders
intersection_zone Rectangle2D::query_zone(const Point2D<int>& pt) const {

    if(pt.x <= corner_points[0].x) {
        if(pt.y <= corner_points[0].y)      return intersection_zone::top_left;
        else if(pt.y >= corner_points[2].y) return intersection_zone::bottom_left;
        else                                return intersection_zone::left;
    }
    else if(pt.x >= corner_points[2].x) {
        if(pt.y <= corner_points[0].y)      return intersection_zone::top_right;
        else if(pt.y >= corner_points[2].y) return intersection_zone::bottom_right;
        else                                return intersection_zone::right;
    }
    else {
        if(pt.y <= corner_points[0].y)      return intersection_zone::top;
        else                                return intersection_zone::bottom;
    }
}

// Pre_condition: pt_in is inside the rectangle (including the borders), pt_out is outside the rectangle
bool Rectangle2D::perfom_intersection(const Point2D<int>& pt_in, Point2D<int>& pt_out) const {

    switch(query_zone(pt_out)) {
    case intersection_zone::top_left:
        // intersect with the left and top segments and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::left)) return true;
        if(intersect_with_segment(pt_in, pt_out, segment_type::top))  return true;
        break;
    case intersection_zone::top_right:
        // intersect with the right and top segments and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::right)) return true;
        if(intersect_with_segment(pt_in, pt_out, segment_type::top))   return true;
        break;
    case intersection_zone::top:
        // intersect with the top segment and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::top)) return true;
        break;
    case intersection_zone::bottom_left:
        // intersect with the left and bottom segments and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::left))   return true;
        if(intersect_with_segment(pt_in, pt_out, segment_type::bottom)) return true;
        break;
    case intersection_zone::bottom_right:
        // intersect with the right and bottom segments and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::right))  return true;
        if(intersect_with_segment(pt_in, pt_out, segment_type::bottom)) return true;
        break;
    case intersection_zone::bottom:
        // intersect with the bottom segment and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::bottom)) return true;
        break;
    case intersection_zone::left:
        // intersect with the left segment and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::left)) return true;
        break;
    case intersection_zone::right:
        // intersect with the right segment and update pt_out with the intersecting point
        if(intersect_with_segment(pt_in, pt_out, segment_type::right)) return true;
        break;
    deafult:
        break;
    }
    std::cout << "ERROR: Segment does not intersect with the rectangle" << std::endl;
    return false;
}

bool Rectangle2D::intersect_with_segment(const Point2D<int>& pt_in, Point2D<int>& pt_out, segment_type stype) const {

    Vector2D<int> vec1(pt_out - pt_in);
    Vector2D<int> vec2, vec3;
    if(stype == segment_type::left) {
        vec2 = corner_points[0] - corner_points[3]; // left segment
        vec3 = corner_points[3] - pt_in;
    }
    else if(stype == segment_type::right) {
        vec2 = corner_points[2] - corner_points[1]; // right segment
        vec3 = corner_points[1] - pt_in;
    }
    else if(stype == segment_type::top) {
        vec2 = corner_points[1] - corner_points[0]; // top segment
        vec3 = corner_points[0] - pt_in;
    }
    else if(stype == segment_type::bottom) {
        vec2 = corner_points[3] - corner_points[2]; // bottom segment
        vec3 = corner_points[2] - pt_in;
    }
    int denom = vec1.cross(vec2);
    if(denom != 0) {
        float t1 = static_cast<double>(vec3.cross(vec2)) / static_cast<double>(denom);
        if(t1 >= 0 && t1 <= 1) {
            pt_out.x = pt_in.x + static_cast<int>(t1*vec1.x);
            pt_out.y = pt_in.y + static_cast<int>(t1*vec1.y);
            return true;
        }
    }
    return false;
}
