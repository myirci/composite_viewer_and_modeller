#include "Ellipse2D.hpp"

Ellipse2D::Ellipse2D(const Ellipse2D& other) : Ellipse2DLight(other) {

    for(int i = 0; i < 4; ++i)
        points[i] = other.points[i];

    for(int i = 0; i < 6; ++i)
        coeff[i] = other.coeff[i];

}

Ellipse2D& Ellipse2D::operator=(const Ellipse2D& rhs) {

    if(this != &rhs) {
        Ellipse2DLight::operator =(rhs);
        for(int i = 0; i < 4; ++i)
            points[i] = rhs.points[i];
        for(int i = 0; i < 6; ++i)
            coeff[i] = rhs.coeff[i];
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
    // rot angle [-PI, PI)
    osg::Vec2d vec_mj = points[1] - points[0];
    if(vec_mj.y() < 0)
        rot_angle = -std::acos(vec_mj.x() / vec_mj.length());
    else
        rot_angle = std::acos(vec_mj.x() / vec_mj.length());
}

void Ellipse2D::update_minor_axis(const osg::Vec2d& pt2) {

    points[2] = pt2;
    osg::Vec2d vec = center - points[2];
    points[3] = center + vec;
    smn_axis = vec.length();
}

void Ellipse2D::calculate_coefficients_from_parameters() {

    double as = smj_axis*smj_axis;
    double bs = smn_axis*smn_axis;
    coeff[0] = 0.5 * (as + bs + cos(2*rot_angle) * (bs - as));
    coeff[1] = sin(2*rot_angle) * (bs - as);
    coeff[2] = 0.5 * (as + bs - cos(2*rot_angle) * (bs - as));
    coeff[3] = -2 * center.x() * coeff[0] - center.y() * coeff[1];
    coeff[4] = -2 * center.y() * coeff[2] - center.x() * coeff[1];
    coeff[5] = center.x() * center.x() * coeff[0] + center.x() * center.y() * coeff[1] + center.y() * center.y() * coeff[2] - as * bs;
}

// the coefficents mmust be updated before calling this function
void Ellipse2D::calculate_parameters_from_coeffients() {

    double v1 = coeff[1] * coeff[1] - 4 * coeff[0] * coeff[2];
    double v2 = 2.0 * coeff[2] * coeff[3] - coeff[1] * coeff[4];
    double v3 = 2.0 * coeff[0] * coeff[4] - coeff[1] * coeff[3];
    double v4 = 0.5 * coeff[0] * coeff[4] * coeff[4] +
                0.5 * coeff[2] * coeff[3] * coeff[3] +
                0.5 * coeff[5] * coeff[1] * coeff[1] -
                0.5 * coeff[1] * coeff[3] * coeff[4] -
                2.0 * coeff[0] * coeff[2] * coeff[5];
    double v5 = coeff[0] - coeff[2];
    double v6 = v5*v5 + coeff[1] * coeff[1];
    double v7 = coeff[0] + coeff[2];
    double v8 = coeff[1] / v5;

    center.x() = v2 / v1;
    center.y() = v3 / v1;
    smj_axis = std::sqrt(v4/(0.25*v1*(std::sqrt(v6) - v7)));
    smn_axis = std::sqrt(v4/(0.25*v1*(-std::sqrt(v6) - v7)));

    if(coeff[1] == 0) {
        if(coeff[0] < coeff[2])       rot_angle = 0.0;     // theta = 0
        else if (coeff[0] > coeff[2]) rot_angle = HALF_PI; // theta = 90
        else std::cout << "This is circle dude!" << std::endl;
    }
    else {
        double theta1 = atan(v8);
        if(theta1 < 0) theta1 += PI;
        theta1 /= 2.0;
        double theta2 = theta1 + HALF_PI;
        if(coeff[1] < 0)  // 0 < theta < 90
           if(coeff[0] != coeff[2]) rot_angle = std::min(theta1, theta2);
        else              // 90 < theta < 180
            if(coeff[0] != coeff[2]) rot_angle = std::max(theta1, theta2);
    }

    calculate_axes_end_points();
}

void Ellipse2D::calculate_center_from_coefficients() {

    double v1 = coeff[1] * coeff[1] - 4 * coeff[0] * coeff[2];
    double v2 = 2.0 * coeff[2] * coeff[3] - coeff[1] * coeff[4];
    double v3 = 2.0 * coeff[0] * coeff[4] - coeff[1] * coeff[3];
    center.x() = v2 / v1;
    center.y() = v3 / v1;
}

void Ellipse2D::calculate_semiaxes_from_coefficients() {

    double v1 = coeff[1] * coeff[1] - 4 * coeff[0] * coeff[2];
    double v4 = 0.5 * coeff[0] * coeff[4] * coeff[4] +
            0.5 * coeff[2] * coeff[3] * coeff[3] +
            0.5 * coeff[5] * coeff[1] * coeff[1] -
            0.5 * coeff[1] * coeff[3] * coeff[4] -
            2.0 * coeff[0] * coeff[2] * coeff[5];
    double v5 = coeff[0] - coeff[2];
    double v6 = v5*v5 + coeff[1] * coeff[1];
    double v7 = coeff[0] + coeff[2];
    smj_axis = std::sqrt(v4/(0.25*v1*(std::sqrt(v6) - v7)));
    smn_axis = std::sqrt(v4/(0.25*v1*(-std::sqrt(v6) - v7)));
}

void Ellipse2D::calculate_theta_from_coefficients() {

    double v5 = coeff[0] - coeff[2];
    double v8 = coeff[1] / v5;

    if(coeff[1] == 0) {
        if(coeff[0] < coeff[2])         rot_angle = 0.0;
        else if (coeff[0] > coeff[2])   rot_angle = HALF_PI;
        else std::cout << "This is circle dude!" << std::endl;
    }
    else {
        double theta1 = atan(v8);
        if(theta1 < 0) theta1 += PI;
        theta1 /= 2.0;
        double theta2 = theta1 + HALF_PI;
        if(coeff[1] < 0)  // 0 < theta < 90
           if(coeff[0] != coeff[2]) rot_angle = std::min(theta1, theta2);
        else              // 90 < theta < 180
            if(coeff[0] != coeff[2]) rot_angle = std::max(theta1, theta2);
    }
}

void Ellipse2D::calculate_axes_end_points() {

    double cos_theta = cos(rot_angle);
    double sin_theta = sin(rot_angle);

    osg::Vec2d smj_vec = osg::Vec2d(cos_theta, sin_theta) * smj_axis;
    points[0] = center + smj_vec;
    points[1] = center - smj_vec;

    osg::Vec2d smn_vec = osg::Vec2d(-sin_theta, cos_theta) * smn_axis;
    points[2] = center + smn_vec;
    points[3] = center - smn_vec;
}

void Ellipse2D::get_tangent_points(const osg::Vec2d& dir, osg::Vec2d& left, osg::Vec2d& right) const {

    // check the direction vector
    const float EPSILON = 0.0000001;
    if(std::abs(dir.x()) < EPSILON && std::abs(dir.y()) < EPSILON) {
        std::cerr << "Zero direction vector: " << dir.x() << " " << dir.y() << std::endl;
        return;
    }

    // if the ellipse axes and the coordinate axes are parallel
    double k1 = 4*coeff[0]*coeff[2] - coeff[1]*coeff[1];
    double k2 = 2*coeff[2]*coeff[3] - coeff[1]*coeff[4];
    double k3 = 4*coeff[2]*coeff[5] - coeff[4]*coeff[4];
    double k4 = coeff[1]*coeff[1]*coeff[5] - coeff[1]*coeff[3]*coeff[4] + coeff[2]*coeff[3]*coeff[3];
    double k5 = coeff[0]*dir.x()*dir.x() + coeff[1]*dir.x()*dir.y() + coeff[2]*dir.y()*dir.y();
    // general case
    if(std::abs(dir.x()) > EPSILON && std::abs(dir.y()) > EPSILON) {

        double a = k1*k5;
        double b = 2*k2*k5;
        double c = k4*dir.x()*dir.x() + k3*dir.y()*(coeff[2]*dir.y() + coeff[1]*dir.x());

        double delta_sqrt = std::sqrt(b*b - 4*a*c);
        left.x()  = (-b + delta_sqrt) / (2*a);
        right.x() = (-b - delta_sqrt) / (2*a);

        double k6 = dir.y()*coeff[1] + 2*dir.x()*coeff[0];
        double k7 = dir.y()*coeff[4] + dir.x()*coeff[3];
        double k8 = dir.x()*coeff[1] + 2*dir.y()*coeff[2];
        left.y()  = -(left.x()*k6 + k7) / (k8);
        right.y()  = -(right.x()*k6 + k7) / (k8);
    }
    else {
        if(std::abs(coeff[1]) < EPSILON) {
            if(std::abs(dir.x()) < EPSILON) {      // vertical tangent points
                if(coeff[0] > coeff[2]) { // theta = 90
                    left = points[2];
                    right = points[3];
                }
                else { // theta = 0
                    left = points[0];
                    right = points[1];
                }
            }
            else if(std::abs(dir.y()) < EPSILON) { // horizontal tangent points
                if(coeff[0] > coeff[2]) { // theta = 90
                    left = points[0];
                    right = points[1];
                }
                else { // theta = 0
                    left = points[2];
                    right = points[3];
                }
            }
        }
        else {
            if(std::abs(dir.x()) < EPSILON) {
                // vertical tangent points
                double a = k1;
                double b = 2*k2;
                double c = k3;

                double delta_sqrt = std::sqrt(b*b - 4*a*c);
                left.x()  = (-b + delta_sqrt) / (2*a);
                right.x() = (-b - delta_sqrt) / (2*a);
                left.y()  = (-coeff[1]*left.x() - coeff[4] ) / (2*coeff[2]);
                right.y() = (-coeff[1]*right.x() - coeff[4]) / (2*coeff[2]);
            }
            else {
                // horizontal tangent points
                double a = k1*coeff[0];
                double b = 2*k2*coeff[0];
                double c = k4;

                double delta_sqrt = std::sqrt(b*b - 4*a*c);
                left.x()  = (-b + delta_sqrt) / (2*a);
                right.x() = (-b - delta_sqrt) / (2*a);
                left.y()  = (-2*coeff[0]*left.x() - coeff[3] ) / coeff[1];
                right.y() = (-2*coeff[0]*right.x() - coeff[3]) / coeff[1];
            }
        }
    }

    // arrange points according to their orientation with respect to the line passing through the center of the ellipse.
    osg::Vec2d vec = left - center;
    if((dir.x() * vec.y() - dir.y() * vec.x()) < 0)
        std::swap(left, right);
}

std::ostream& operator<<(std::ostream& out, const Ellipse2D& ellipse) {

    out << static_cast<const Ellipse2DLight&>(ellipse);
    for(int i = 0; i < 4; ++i)
        out << "point-" << i << ": " << ellipse.points[i].x() << " " << ellipse.points[i].y() << std::endl;

    out << " a: " << ellipse.coeff[0]
        << " b: " << ellipse.coeff[1]
        << " c: " << ellipse.coeff[2]
        << " d: " << ellipse.coeff[3]
        << " e: " << ellipse.coeff[4]
        << " f: " << ellipse.coeff[5]
        << std::endl;
    return out;
}
