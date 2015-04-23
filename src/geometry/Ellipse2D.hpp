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

    void update_major_axis(const osg::Vec2d& pt0, const osg::Vec2d& pt1);
    void update_minor_axis(const osg::Vec2d& pt2);
    void calculate_coefficients_from_parameters();
    void rotate(double angle);
    void translate(const osg::Vec2d& translation_vec);

    osg::Vec2d points[4];

    /* ax^2 + bxy + cy^2 + dx + ey + f = 0
     * coeff[0] : a; coeff[1] : b; coeff[2] : c; coeff[3] : d; coeff[4] : e; coeff[5] : f
    */
    double* coeff;
};

std::ostream& operator<<(std::ostream& out, const Ellipse2D& ellipse);

#endif // ELLIPSE2D_HPP
