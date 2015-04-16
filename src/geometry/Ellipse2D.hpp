#ifndef ELLIPSE2D_HPP
#define ELLIPSE2D_HPP

#include "Ellipse2DLight.hpp"

struct Ellipse2D : public Ellipse2DLight {

    Ellipse2D(double smj = 0.0, double smn = 0.0, double rot = 0.0, const osg::Vec2d& center_pt = osg::Vec2d(0,0)) :
        Ellipse2DLight(smj, smn, rot, center_pt), coeff(nullptr) { }

    ~Ellipse2D() {
        if(coeff != nullptr) {
            delete[] coeff;
            coeff = nullptr;
        }
    }
    Ellipse2D& operator=(const Ellipse2D& rhs);
    // from screen-coordinates to projected coordinates
    void calculate_algebraic_equation_in_projected_coordinates(int w, int h, double n, double half_fovy);   
    void rotate(double angle);
    void translate(const osg::Vec2d& translation_vec);
    void update_major_axis(const osg::Vec2d& pt0, const osg::Vec2d& pt1);
    void update_minor_axis(const osg::Vec2d& pt2);

    osg::Vec2d points[4];
    double* coeff;
};

std::ostream& operator<<(std::ostream& out, const Ellipse2D& ellipse);

#endif // ELLIPSE2D_HPP
