#ifndef LINE2D_HPP
#define LINE2D_HPP

#include <cmath>

struct Line2D {
    Line2D(double m_, double n_) : m(m_), n(n_) { }             // y = mx + n
    Line2D(double a, double b, double c) : m(-a/b), n(-c/b) { } // ax + by + c = 0
    Line2D(double x1, double y1, double x2, double y2) {        // passing two points (x1,y1) and (x2,y2)
        m = (y2-y1) / (x2-x1);
        n = y1 - m*x1;
    }
    double distance_to_line_from_a_point(double x0, double y0) {
        return std::abs(m*x0-y0+n)/std::sqrt(m*m + 1);
    }
    double evaluate(double x, double y) { return (y - m*x - n); }
    double m;
    double n;
};

#endif // LINE2D_HPP
