#include "Ellipse2D.hpp"

// Not a full assignment: coefficients are not copied
Ellipse2D& Ellipse2D::operator=(const Ellipse2D& rhs) {

    if(this != &rhs) {
        Ellipse2DLight::operator =(rhs);
        for(int i = 0; i < 4; ++i)
            points[i] = rhs.points[i];
    }
    return *this;
}

void Ellipse2D::translate(const osg::Vec2d& translation_vec) {

    center += translation_vec;
    for(int i = 0; i < 4; ++i)
        points[i] += translation_vec;
}

// "angle" is the angle value in radian which is used to rotate the ellipse around its center in
// counter clockwise direction.
void Ellipse2D::rotate(double angle) {

    rot_angle += angle;
    // First translate the center such that center coincides with the origin
    osg::Vec2d current_center = center;
    translate(-current_center);

    double c = cos(angle);
    double s = sin(angle);

    double x, y;
    for(int i = 0; i < 4; ++i) {
        x = c * points[i].x() - s * points[i].y();
        y = s * points[i].x() + c * points[i].y();
        points[i].x() = x;
        points[i].y() = y;
    }

    translate(current_center);
}

void Ellipse2D::update_major_axis(const osg::Vec2d& pt0, const osg::Vec2d& pt1) {

    // update major axis end-points
    points[0] = pt0;
    points[1] = pt1;

    // update center
    center = (points[0] + points[1]) / 2.0;

    // update semi major axis length
    smj_axis = (points[0] - center).length();

    // rotation angle is the angle between the major axis and the the positive x-axis
    osg::Vec2d vec_mj = points[1] - points[0];
    if(vec_mj.y() < 0)
        rot_angle = -std::acos(vec_mj.x()/vec_mj.length());
    else
        rot_angle = std::acos(vec_mj.x()/vec_mj.length());
}

void Ellipse2D::update_minor_axis(const osg::Vec2d& pt2) {

    points[2] = pt2;
    osg::Vec2d vec = center - points[2];
    points[3] = center + vec;
    smn_axis = vec.length();
}

// The calculated ellipse is on the near clipping plane of the camera.
void Ellipse2D::calculate_algebraic_equation_in_projected_coordinates(int w, int h, double n, double half_fovy) {

    /* Ellipse is defined in the logical device coordinates in which the location of the origin is the lower left
     * corner of the viewport and increasing x and y coordinates are towards right and top relatively. The unit of
     * coordinates is pixel.
    */

    /* ax^2 + bxy + cy^2 + dx + ey + f = 0
     * coeff[0] : a; coeff[1] : b; coeff[2] : c; coeff[3] : d; coeff[4] : e; coeff[5] : f
    */
    // If coefficient array has not been defined before, allocate memory
    if(coeff == nullptr)
        coeff = new double[6];

    /*  From viewport coordinates we can calculate the normalized device coordinates
     *
     *  x_ndc = (2/VP_width)*x_vp - x0
     *  y_ndc = (2/VP_height)*y_vp - y0
     *
     *  (x0,y0) is the lower left corner coordinates of the view port on the screen and (x0,y0) = (0,0) in this application.
     *  Furthermore, VP_width and VP_height are always equal to the width and height of the screen. Viewport coordinate frame
     *  located in the middle of the viewport.
     *
     *  1) Convert logical device coordinates to viewport coordinates
     *
     *  x_vp = x_log - VP_width/2
     *  y_vp = y_log - VP_height/2
     *
     *  2) Convert viewport coordinates to normalized device coordinates
     *
     *  x_ndc = (2/VP_width)*x_log - 1 - x0
     *  y_ndc = (2/VP_height)*y_log - 1 - y0
    */
    double consx = 2.0/static_cast<double>(w);
    double consy = 2.0/static_cast<double>(h);

    Point2D<double> points_ndc[4];
    points_ndc[0] = Point2D<double>(consx * points[0].x() - 1, consy * points[0].y() - 1);
    points_ndc[1] = Point2D<double>(consx * points[1].x() - 1, consy * points[1].y() - 1);
    points_ndc[2] = Point2D<double>(consx * points[2].x() - 1, consy * points[2].y() - 1);
    points_ndc[3] = Point2D<double>(consx * points[3].x() - 1, consy * points[3].y() - 1);

    /* Normalized device coordinates are in the [-1,1] range for both x and y coordinates. We can calculate the projected
     * coordinates. Eye (object coordinates in camera space) coordinates are projected (perspective or othographic) on to
     * the near clipping plane and the coordinates of thees points are called projected coordinates.
     * In the case of symmetric viewing volume,
     *
     * x_prj = x_ndc * aspect * n * tan(FOVy/2)
     * y_prj = y_ndc * n * tan(FOVy/2)
     *
    */

    double aspect = static_cast<double>(w)/static_cast<double>(h);
    double c1 = n * tan(half_fovy);
    double c2 = aspect * c1;

    Point2D<double> points_projected[4];
    points_projected[0] = Point2D<double>(points_ndc[0].x*c2, points_ndc[0].y*c1);
    points_projected[1] = Point2D<double>(points_ndc[1].x*c2, points_ndc[1].y*c1);
    points_projected[2] = Point2D<double>(points_ndc[2].x*c2, points_ndc[2].y*c1);
    points_projected[3] = Point2D<double>(points_ndc[3].x*c2, points_ndc[3].y*c1);

    Point2D<double> center_projected((points_projected[0].x + points_projected[1].x)/2.0,
                                     (points_projected[0].y + points_projected[1].y)/2.0);
    double smj = dist(points_projected[0], points_projected[1]) / 2.0;
    double smn = dist(points_projected[2], center_projected);
    double as = smj * smj;
    double bs = smn * smn;
    coeff[0] = 0.5*(as + bs + cos(2*rot_angle)*(bs - as));
    coeff[1] = sin(2*rot_angle)*(bs - as);
    coeff[2] = 0.5*(as + bs - cos(2*rot_angle)*(bs - as));
    coeff[3] = -2*center_projected.x*coeff[0] - center_projected.y*coeff[1];
    coeff[4] = -2*center_projected.y*coeff[2] - center_projected.x*coeff[1];
    coeff[5] = center_projected.x*center_projected.x*coeff[0] + center_projected.x*center_projected.y*coeff[1] + center_projected.y*center_projected.y*coeff[2] - as*bs;
}

std::ostream& operator<<(std::ostream& out, const Ellipse2D& ellipse) {

    out << static_cast<const Ellipse2DLight&>(ellipse);
    for(int i = 0; i < 4; ++i)
        out << "point-" << i << ": " << ellipse.points[i].x() << " " << ellipse.points[i].y() << std::endl;

    if(ellipse.coeff != nullptr) {
        out << "a: " << ellipse.coeff[0]
            << " b: " << ellipse.coeff[1]
            << " c: " << ellipse.coeff[2]
            << " d: " << ellipse.coeff[3]
            << " e: " << ellipse.coeff[4]
            << " f: " << ellipse.coeff[5]
            << std::endl;
    }
    return out;
}
