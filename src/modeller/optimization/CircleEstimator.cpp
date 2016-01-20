#include "CircleEstimator.hpp"
#include "../../geometry/Ellipse2D.hpp"
#include "../../geometry/Circle3D.hpp"
#include "../../geometry/Segment2D.hpp"
#include "../../geometry/Ray3D.hpp"
#include "../../geometry/Plane3D.hpp"
#include "ExtractPlaneNormals.hpp"
#include "../ProjectionParameters.hpp"
#include "../../osg/OsgUtility.hpp"

// PUBLIC METHODS

/*
// Method-1: Needed to solve a nonlinear optimization problem for estimating the orientation of the circles
void CircleEstimator::estimate_3d_circles_with_fixed_depth_method1(const Ellipse2D& ellipse, Circle3D* circles, ProjectionParameters const * const pp, double desired_depth) {

    if(desired_depth > -pp->near || desired_depth < -pp->far) {
        std::cout << "ERROR: desired depth must be within [-near, -far]: [" << -pp->near << ", " << -pp->far << "]" << std::endl;
    }

    estimate_unit_3d_circles_method1(ellipse, circles, pp);
    for(int i = 0; i < 2; ++i) {
        if(circles[i].center[2] > 0) {
            circles[i].center = - circles[i].center;
            circles[i].normal = - circles[i].normal;
        }
        if(desired_depth != circles[i].center[2]) {
            circles[i].radius = circles[i].radius * (desired_depth/circles[i].center[2]);
        }
        else {
            continue;
        }
        circles[i].center *= circles[i].radius;
    }
}

void CircleEstimator::estimate_3d_circles_with_fixed_radius_method1(const Ellipse2D& ellipse, Circle3D* circles, ProjectionParameters const * const pp, double desired_radius) {

    estimate_unit_3d_circles_method1(ellipse, circles, pp);
    for(int i = 0; i < 2; ++i) {
        circles[i].center *= desired_radius;
        circles[i].radius = desired_radius;
        if(circles[i].center[2] > -pp->near) {
            circles[i].center = - circles[i].center;
            circles[i].normal = - circles[i].normal;
        }
    }
}
*/

// Method-2: Based on analytical solution
int CircleEstimator::estimate_3d_circles_with_fixed_depth(const Ellipse2D& ellipse, Circle3D* circles, const ProjectionParameters* const pp, double desired_depth) {

    if(desired_depth > -pp->near || desired_depth < -pp->far)
        std::cout << "ERROR: desired depth must be within [-near, -far]: [" << -pp->near << ", " << -pp->far << "]" << std::endl;

    int count = estimate_unit_3d_circles(ellipse, circles, pp);
    for(int i = 0; i < count; ++i) {
        if(desired_depth != circles[i].center[2])
            circles[i].radius = desired_depth / circles[i].center[2];
        else
            continue;
        circles[i].center *= circles[i].radius;
    }
    return count;
}

int CircleEstimator::estimate_3d_circles_with_fixed_radius(const Ellipse2D& ellipse, Circle3D* circles, const ProjectionParameters* const pp, double desired_radius) {

    int count = estimate_unit_3d_circles(ellipse, circles, pp);
    for(int i = 0; i < count; ++i) {
        circles[i].center *= desired_radius;
        circles[i].radius = desired_radius;
    }
    return count;
}

void CircleEstimator::estimate_3d_circles_under_orthographic_projection(const Ellipse2D& ellipse, Circle3D& circle, double near) {

    // 2d semi-major and semi-minor axis vectors
    osg::Vec2d smj_vec_2d = ellipse.points[1] - ellipse.center;
    osg::Vec2d smn_vec_2d = ellipse.points[2] - ellipse.center;

    // 3d semi-major axis: the given ellipse must be on the image plane (z = k). Thus, the z-ccordinate of the 3d semi-major
    // axis must be zero.
    osg::Vec3 smj_vec_3d(smj_vec_2d.x(), smj_vec_2d.y(), 0);

    /* tilted 3d semi-minor axis vector must have the same magnitude with the 3d semi-major axis. In otder to have this,
     *
     * smj_vec_2d = (ux,uy)  smj_vec_3d = (ux,uy,0)
     * smn_vec_2d = (vx,vy)  smn_vec_3d = (vx,vy,0)  tilted_smn_vec_3d = (vx,vy,z)
     *
     * length²(tilted_smn_vec_3d) = length²(smj_vec_3d)
     * vx² + vy² + z² = ux² + uy²
     * z² = ux² + uy² - vx² - vy²
     * z = +/- sqrt(ux² + uy² - vx² - vy²)
     *
     * we select the positive z, which means tilt of the minor axis vector is always towards the camera. (Camera is located at
     * the origin, image plane is located at the negative side of the z-axis and z-coordinate of the tilted semi_minor axis is
     * selected as positive)
     *
     */
    osg::Vec3 tilted_smn_vec_3d(smn_vec_2d.x(), smn_vec_2d.y(), std::sqrt(smj_vec_2d*smj_vec_2d - smn_vec_2d*smn_vec_2d));

    // orientation of the circle : normal of the circle
    osg::Vec3 n = smj_vec_3d ^ tilted_smn_vec_3d;
    n.normalize();
    circle.normal[0] = n.x();
    circle.normal[1] = n.y();
    circle.normal[2] = n.z();

    // radius of the circle is equal to the length of the smj_vec_2d = smj_vec_3d = tilted_smn_vec_3d = ellipse.smj_axis on the
    // near plane (image plane).
    circle.radius =  ellipse.smj_axis;

    // center of the ellipse on the image plane is (ellipse.center.x(), ellipse.center.y(), near)
    circle.center[0] = ellipse.center.x();
    circle.center[1] = ellipse.center.y();
    circle.center[2] = near;
}


// PRIVATE METHODS
/*
void CircleEstimator::estimate_unit_3d_circles_method1(const Ellipse2D& ellipse, Circle3D* circles, const ProjectionParameters* const pp) {

    // we need to negate the near value, because in opengl near and far values are positive. We need the actual value
    // in mathemetical calculations.

    double near = -pp->near;
    Eigen::Matrix3d M;
    M << ellipse.coeff[0],          ellipse.coeff[1]/2.0,      ellipse.coeff[3]/(2*near),
         ellipse.coeff[1]/2.0,      ellipse.coeff[2],          ellipse.coeff[4]/(2*near),
         ellipse.coeff[3]/(2*near), ellipse.coeff[4]/(2*near), ellipse.coeff[5]/(near*near);

    Eigen::Vector3d n1, n2;
    solve_for_normals(M, n1, n2);

    Eigen::Matrix3d E1, E2;
    construct_change_of_basis_matrix(E1, n1);
    construct_change_of_basis_matrix(E2, n2);

    Eigen::Matrix3d Mp1 = E1.transpose()*M*E1;
    Eigen::Matrix3d Mp2 = E2.transpose()*M*E2;

    int sign = (Mp1(0,0) >= 0 ? 1 : -1);
    double c1(0.0), c2(0.0);
    c1 = sign / std::sqrt(Mp1(0,2)*Mp1(0,2) + Mp1(1,2)*Mp1(1,2) - Mp1(2,2)*Mp1(0,0));
    sign = (Mp2(0,0) >= 0 ? 1 : -1);
    c2 = sign / std::sqrt(Mp2(0,2)*Mp2(0,2) + Mp2(1,2)*Mp2(1,2) - Mp2(2,2)*Mp2(0,0));

    circles[0].center = Eigen::Vector3d(-Mp1(0,2), -Mp1(1,2), Mp1(0,0));
    circles[1].center = Eigen::Vector3d(-Mp2(0,2), -Mp2(1,2), Mp2(0,0));
    circles[0].center *= c1;
    circles[0].radius = 1.0;
    circles[1].center *= c2;
    circles[1].radius = 1.0;
    circles[0].center = E1*circles[0].center;
    circles[1].center = E2*circles[1].center;
    circles[0].normal = n1;
    circles[1].normal = n2;
}

void CircleEstimator::construct_change_of_basis_matrix(Eigen::Matrix3d& mat, const Eigen::Vector3d& vec2) {

    Eigen::Vector3d vec0, vec1;
    if(std::abs(vec2(0)) < std::abs(vec2(1)) && std::abs(vec2(0)) < std::abs(vec2(2))) {
        vec0(0) = 0;
        vec0(1) = vec2(2);
        vec0(2) = -vec2(1);
    }
    else if(std::abs(vec2(1)) < std::abs(vec2(0)) && std::abs(vec2(1)) < std::abs(vec2(2))) {
        vec0(0) = -vec2(2);
        vec0(1) = 0;
        vec0(2) = vec2(0);
    }
    else {
        vec0(0) = vec2(1);
        vec0(1) = -vec2(0);
        vec0(2) = 0;
    }
    vec1 = vec2.cross(vec0);
    vec0.normalize();
    vec1.normalize();

    mat(0, 0) = vec0(0);
    mat(1, 0) = vec0(1);
    mat(2, 0) = vec0(2);

    mat(0, 1) = vec1(0);
    mat(1, 1) = vec1(1);
    mat(2, 1) = vec1(2);

    mat(0, 2) = vec2(0);
    mat(1, 2) = vec2(1);
    mat(2, 2) = vec2(2);
}
*/

int CircleEstimator::estimate_unit_3d_circles(const Ellipse2D& ellipse, Circle3D* circles, const ProjectionParameters *const pp) {

    // Step-1: Construct the associated quadratic form matrix of the 3D cone.
    /*
     * 'ellipse' is the intersection of the 3D cone (whose vertex is at origin) and a plane z = k.
     * The plane is the near clipping plane, since the ellipse is on the near clipping plane.
     * k = -near (recall that near > 0).
     *
     * The cone is constructed in XYZ coordinate frame.
     */

    double near = -pp->near;
    Eigen::Matrix3d Q;
    Q << ellipse.coeff[0],          ellipse.coeff[1]/2.0,      ellipse.coeff[3]/(2*near),
         ellipse.coeff[1]/2.0,      ellipse.coeff[2],          ellipse.coeff[4]/(2*near),
         ellipse.coeff[3]/(2*near), ellipse.coeff[4]/(2*near), ellipse.coeff[5]/(near*near);

    // Step-2: Find the eigenvalues and eigenvectors of the matrix Q.
    /*
     * A real symmetric matrix is self-adjoint. Thus we use the SelfAdjointEigenSolver.
     * Furthermore, a real symmetric matrix has real eigenvalues. The obtained eigenvalues
     * are ordered from smallest to biggest. Since Q represents a cone, one of the eigenvalues
     * must have a different sign than the other two. Multipliying Q with -1 cahnges the signs
     * of eigenvalues. One can obtain two positif and one negatif eigenvalues.
    */

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(Q);
    if(eigensolver.info() != Eigen::Success)
        std::cout << "ERROR: Eigen solver is not successful!" << std::endl;

    if(!check_eigenvalue_constraints(eigensolver.eigenvalues())) {
        Q *= -1;
        eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(Q);
        if(!check_eigenvalue_constraints(eigensolver.eigenvalues()))
            std::cout << "ERROR: Eigenvales does not macth with the ellipse constraints" << std::endl;
    }

    // Step-3: Construct the change of variables matrix from the eigenvectors of Q
    /*
    * From three eigenvectors we can construct 6 different P matrix for eigen decomposition of Q.
    * Q = PDP^T. Determinant of three of these matrices are equal to 1, and determinant of the rest is
    * equal to -1. The change of variable matrix (P) must be a rotation. That is: det(P) must be equal to 1.
    * Moreover, we know that the negative eigenvalue is related with the z coordinate. Thus, we put the eigenvector
    * related to negative eigenvalue to the 3rd column of P. The other two columns are selected based on the
    * determinant = 1 rule.
    *
    * lambda_1: coefficent of x
    * lambda_2: coefficent of y
    * lambda_3: coefficent of z
    *
    * D = diag(lambda_1, lambda_2, lambda_3)
    *
    * After change of variables is applied, we are in the X'Y'Z' coordinate frame.
    *
    */

    Eigen::Matrix3d P;
    P.col(0) << eigensolver.eigenvectors().col(1);
    double lambda_1 = eigensolver.eigenvalues()(1);

    P.col(1) << eigensolver.eigenvectors().col(2);
    double lambda_2 = eigensolver.eigenvalues()(2);

    P.col(2) << eigensolver.eigenvectors().col(0);
    double lambda_3 = eigensolver.eigenvalues()(0);

    if(P.determinant() < 0) {
        P.col(0) << eigensolver.eigenvectors().col(2);
        P.col(1) << eigensolver.eigenvectors().col(1);
        std::swap(lambda_1, lambda_2);
    }

    // Step-4: Calculate the normals and the centers of the 3D circles
    /*
     * The unit normal of the circle is (a,b,c) in X'Y'Z' coordinate frame. There exists
     * four possible circles whose projection matches with the given elllipse. Two pairs
     * of these circles are symmetric with respect to the origin.
     *
    */

    int count = 0;

    if(lambda_1 > lambda_2) {

        // a = +/- k1,  b = 0, c = +/- k2
        double k1 = std::sqrt((lambda_1 - lambda_2)/(lambda_1 - lambda_3));
        double k2 = std::sqrt((lambda_2 - lambda_3)/(lambda_1 - lambda_3));
        double k3 = 1.0 / std::sqrt(-lambda_1*lambda_3);

        // calculate the first center in the x'y'z' coordinate frame for (a > 0 , c > 0)
        Eigen::Vector3d ctr = k3 * Eigen::Vector3d(lambda_3*k1, 0, lambda_1*k2);
        circles[0].center = P * ctr;
        circles[0].normal = P * Eigen::Vector3d(k1, 0, k2);
        circles[0].radius = 1;
        if(circles[0].center[2] > 0) {
            // (a < 0 , c < 0)
            circles[0].center = - circles[0].center;
            circles[0].normal = - circles[0].normal;
        }

        // calculate the second center in the x'y'z' coordinate frame for (a > 0 , c < 0)
        ctr[2] = -ctr[2];
        circles[1].center = P * ctr;
        circles[1].normal = P * Eigen::Vector3d(k1, 0, -k2);
        circles[1].radius = 1;
        if(circles[1].center[2] > 0) {
            // (a < 0 , c > 0)
            circles[1].center = - circles[1].center;
            circles[1].normal = - circles[1].normal;
        }

        count = 2;
    }
    else if(lambda_2 > lambda_1) {

        // a = 0,  b = +/- k1, c = +/- k2
        double k1 = std::sqrt((lambda_2 - lambda_1)/(lambda_2 - lambda_3));
        double k2 = std::sqrt((lambda_1 - lambda_3)/(lambda_2 - lambda_3));
        double k3 = 1.0 / std::sqrt(-lambda_2*lambda_3);

        // calculate the first center in the x'y'z' coordinate frame for (b > 0 , c > 0)
        Eigen::Vector3d ctr = k3 * Eigen::Vector3d(0, lambda_3*k1, lambda_2*k2);
        circles[0].center = P * ctr;
        circles[0].normal = P * Eigen::Vector3d(0, k1, k2);
        circles[0].radius = 1;
        if(circles[0].center[2] > 0) {
            // (b < 0 , c < 0)
            circles[0].center = - circles[0].center;
            circles[0].normal = - circles[0].normal;
        }

        // calculate the second center in the x'y'z' coordinate frame for (b > 0 , c < 0)
        ctr[2] = -ctr[2];
        circles[1].center = P * ctr;
        circles[1].normal = P * Eigen::Vector3d(0, k1, -k2);
        circles[1].radius = 1;
        if(circles[1].center[2] > 0) {
            // (b < 0 , c > 0)
            circles[1].center = - circles[1].center;
            circles[1].normal = - circles[1].normal;
        }
        count = 2;
    }
    else {
        circles[0].center = P * Eigen::Vector3d(0, 0, std::sqrt(-lambda_1/lambda_3));
        circles[0].normal = P * Eigen::Vector3d(0, 0, 1);
        circles[0].radius = 1;
        count = 1;
    }
    return count;
}

bool CircleEstimator::check_eigenvalue_constraints(const Eigen::Vector3d& eigenvalues) {

    if(eigenvalues(0) == 0 || eigenvalues(1) == 0 || eigenvalues(2) == 0) {
        std::cout << "ERROR: Eigenvalue equal to zero" << std::endl;
        return false;
    }
    unsigned char pos(0), neg(0);
    eigenvalues(0) > 0 ? ++pos : ++neg;
    eigenvalues(1) > 0 ? ++pos : ++neg;
    eigenvalues(2) > 0 ? ++pos : ++neg;
    return (pos == 2 && neg == 1) ? true : false;
}

/*
 * The normal of the circle (circle.normal) has to be set in advance.
 * The depth of the circle (circle.center[2]) has to be set to the desired depth in advance.
*/
void CircleEstimator::estimate_3d_circle_from_major_axis_when_circle_depth_is_fixed(const Segment2D& seg, double near, Circle3D& circle) {

    double desired_depth = circle.center[2];
    estimate_unit_3d_circle_from_major_axis(seg, near, circle);
    circle.radius = desired_depth / circle.center[2];
    circle.center *= circle.radius;
}

/*
 * The normal of the circle (circle.normal) has to be set in advance.
 * The radius of the circle (circle.radius) has to be set to the desired radius in advance.
*/
void CircleEstimator::estimate_3d_circle_from_major_axis_when_circle_radius_is_fixed(const Segment2D& seg, double near, Circle3D& circle) {

    double desired_radius = circle.radius;
    estimate_unit_3d_circle_from_major_axis(seg, near, circle);
    circle.center *= desired_radius;
    circle.radius = desired_radius;
}

void CircleEstimator::estimate_3d_circles_using_orthogonality_constraint(const Ellipse2D& ellipse, double near, Circle3D* circles, bool use_forth_pt) {

    osg::Vec2d pt4 = ellipse.points[2] - ellipse.points[3];
    pt4.normalize();
    std::vector<osg::Vec3d> lframe3d {osg::Vec3d(ellipse.points[2], near),
                                      osg::Vec3d(ellipse.points[0], near),
                                      osg::Vec3d(ellipse.points[1], near)};
    if(use_forth_pt) {
        lframe3d.push_back(osg::Vec3d(ellipse.points[3], near));
    }
    else {
        lframe3d.push_back(osg::Vec3d(ellipse.points[2] +  pt4, near));
    }

    double f00 = lframe3d[0] * lframe3d[0];
    double f01 = lframe3d[0] * lframe3d[1];
    double f02 = lframe3d[0] * lframe3d[2];
    double f03 = lframe3d[0] * lframe3d[3];
    double f12 = lframe3d[1] * lframe3d[2];
    double f13 = lframe3d[1] * lframe3d[3];
    double f23 = lframe3d[2] * lframe3d[3];

    double a = f01*f02*f13 + f01*f03*f12 - f00*f12*f13 - f01*f01*f23;
    double b = 2*f01*(f00*f23 - f02*f03);
    double c = f00*(f02*f03 - f00*f23);

    // Z1 = k1*Z0, Z2 = k2*Z0, Z3 = k3*Z0
    double dlt = b*b - 4*a*c;
    double k1_1 = (-b + std::sqrt(dlt)) / (2*a);
    double k2_1 = (k1_1*f01 - f00) / (k1_1*f12 - f02);
    double k3_1 = (k1_1*f01 - f00) / (k1_1*f13 - f03);
    double k1_2 = (-b - std::sqrt(dlt)) / (2*a);
    double k2_2 = (k1_2*f01 - f00) / (k1_2*f12 - f02);
    double k3_2 = (k1_2*f01 - f00) / (k1_2*f13 - f03);

    double Z0 = near;

    osg::Vec3d P0(Z0 * lframe3d[0].x() / near, Z0 * lframe3d[0].y() / near, Z0);

    double Z1_1 = k1_1 * Z0;
    osg::Vec3d P1_1(Z1_1 * lframe3d[1].x() / near, Z1_1 * lframe3d[1].y() / near, Z1_1);
    double Z1_2 = k1_2 * Z0;
    osg::Vec3d P1_2(Z1_2 * lframe3d[1].x() / near, Z1_2 * lframe3d[1].y() / near, Z1_2);

    double Z2_1 = k2_1 * Z0;
    osg::Vec3d P2_1(Z2_1 * lframe3d[2].x() / near, Z2_1 * lframe3d[2].y() / near, Z2_1);
    double Z2_2 = k2_2 * Z0;
    osg::Vec3d P2_2(Z2_2 * lframe3d[2].x() / near, Z2_2 * lframe3d[2].y() / near, Z2_2);

    double Z3_1 = k3_1 * Z0;
    osg::Vec3d P3_1(Z3_1 * lframe3d[3].x() / near, Z3_1 * lframe3d[3].y() / near, Z3_1);
    double Z3_2 = k3_2 * Z0;
    osg::Vec3d P3_2(Z3_2 * lframe3d[3].x() / near, Z3_2 * lframe3d[3].y() / near, Z3_2);

    osg::Vec3d ctr = (P1_1 + P2_1) / 2.0;
    circles[0].center[0] = ctr.x();
    circles[0].center[1] = ctr.y();
    circles[0].center[2] = ctr.z();

    osg::Vec3d nrm = (P3_1 - P0);
    nrm.normalize();
    circles[0].normal[0] = nrm.x();
    circles[0].normal[1] = nrm.y();
    circles[0].normal[2] = nrm.z();

    circles[0].radius = (ctr - P0).length();

    ctr = (P1_2 + P2_2) / 2.0;
    circles[1].center[0] = ctr.x();
    circles[1].center[1] = ctr.y();
    circles[1].center[2] = ctr.z();

    nrm = (P3_2 - P0);
    nrm.normalize();
    circles[1].normal[0] = nrm.x();
    circles[1].normal[1] = nrm.y();
    circles[1].normal[2] = nrm.z();
    circles[1].radius = (ctr - P0).length();
}

/*
 * The normal of the circle (circle.normal) has to be set in advance
*/
void CircleEstimator::estimate_unit_3d_circle_from_major_axis(const Segment2D& seg, double near, Circle3D& circle) {

    osg::Vec3d v1(seg.pt1, near);
    osg::Vec3d v2(seg.pt2, near);
    osg::Vec3d diff = v2 - v1;

    osg::Vec3d n0(circle.normal[0], circle.normal[1], circle.normal[2]);
    osg::Vec3d n1 = n0 ^ (v2 ^ v1);
    osg::Vec3d n2(-near*diff.x(), -near*diff.y(), v1 * diff);
    osg::Vec3d n3(n2.x(), n2.y(), v2 * diff);

    double par1 = n0 * v1;
    double par2 = n0 * v2;
    double par3 = n0 * n0;
    double par4 = n0 * n2;
    double par5 = n0 * n3;

    double k1 = (par2 * par2 * (v1 * v1) - par1 * par1 * (v2 * v2) ) / (2 * near * par1 );
    double k2 = std::sqrt(par3 * (n2 * n2) - par4 * par4);
    double k3 = std::sqrt(par3 * (n3 * n3) - par5 * par5);

    osg::Vec3d cam_centre(0,0,0);
    Ray3D ray1(cam_centre, v1);
    Ray3D ray2(cam_centre, v2);
    if(!is_intersecting(ray2, Plane3D(n2.x(), n2.y(), n2.z(), k2))) k2 *= -1;
    if(!is_intersecting(ray1, Plane3D(n3.x(), n3.y(), n3.z(), k3))) k3 *= -1;

    // Solve the plane equations parametrically
    Eigen::Matrix3d mat0;
    mat0 << n1.x(), n1.y(), n1.z(),
            n2.x(), n2.y(), n2.z(),
            n3.x(), n3.y(), n3.z();

    Eigen::Matrix3d mat1;
    mat1 << 0, n1.y(), n1.z(),
            k2, n2.y(), n2.z(),
            k3, n3.y(), n3.z();

    Eigen::Matrix2d mat2;
    mat2 << n2.y(), n2.z(),
            n3.y(), n3.z();

    Eigen::Matrix3d mat3;
    mat3 << n1.x(), 0, n1.z(),
            n2.x(), k2, n2.z(),
            n3.x(), k3, n3.z();

    Eigen::Matrix2d mat4;
    mat4 << n2.x(), n2.z(),
            n3.x(), n3.z();

    Eigen::Matrix3d mat5;
    mat5 << n1.x(), n1.y(), 0,
            n2.x(), n2.y(), k2,
            n3.x(), n3.y(), k3;

    // X = s1 + s2*Z2, Y = s3 + s4*Z2, Z = s5
    double D = mat0.determinant();
    double s1 = -mat1.determinant() / D;
    double s2 = -k1 * mat2.determinant() / D;
    double s3 = -mat3.determinant() / D ;
    double s4 = k1 * mat4.determinant() / D;
    double s5 = -mat5.determinant() / D;

    double par6 = s2 - v2.x() / v2.z();
    double par7 = s4 - v2.y() / v2.z();
    double a = par6*par6 + par7*par7 + 1;
    double b = 2*(s1*par6 + s3*par7 - s5);

    /*
    double c = s1*s1 + s3*s3 + s5*s5 - 1;
    double disc1 = b*b - 4*a*c;
    double disc2 = 4*(a - (s3*(s2 - t1) - s1*(s4 - t2))*(s3*(s2 - t1) - s1*(s4 - t2)) -
                          (s1 + s5*(s2 - t1))*(s1 + s5*(s2 - t1)) - (s3 + s5*(s4 - t2))*(s3 + s5*(s4 - t2)));
    std::cout << "s1: " << s1 << std::endl;
    std::cout << "s2: " << s2 << std::endl;
    std::cout << "s3: " << s3 << std::endl;
    std::cout << "s4: " << s4 << std::endl;
    std::cout << "s5: " << s5 << std::endl;
    std::cout << "Discriminant-1: " << disc1 << std::endl;
    std::cout << "Discriminant-2: " << disc2 << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "var1: " << (s3*(s2 - t1) - s1*(s4 - t2))*(s3*(s2 - t1) - s1*(s4 - t2)) << std::endl;
    std::cout << "var2: " << (s1 + s5*(s2 - t1))*(s1 + s5*(s2 - t1)) << std::endl;
    std::cout << "var3: " << (s3 + s5*(s4 - t2))*(s3 + s5*(s4 - t2)) << std::endl;
    */

    circle.radius = 1;
    double Z2 = -b / (2*a);
    circle.center[0] = s1 + s2 * Z2;
    circle.center[1] = s3 + s4 * Z2;
    circle.center[2] = s5;

}
