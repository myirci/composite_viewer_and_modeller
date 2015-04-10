#include "CircleEstimator.hpp"
#include "../../geometry/Ellipse2D.hpp"
#include "../../geometry/Circle3D.hpp"
#include "ExtractPlaneNormals.hpp"
#include "../PersProjParam.hpp"

// PUBLIC METHODS

/*
// Method-1: Needed to solve a nonlinear optimization problem for estimating the orientation of the circles
void CircleEstimator::estimate_3d_circles_with_fixed_depth_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_depth) {

    if(desired_depth > -ppp->near || desired_depth < -ppp->far) {
        std::cout << "ERROR: desired depth must be within [-near, -far]: [" << -ppp->near << ", " << -ppp->far << "]" << std::endl;
    }

    estimate_unit_3d_circles_method1(ellipse, circles, ppp);
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

void CircleEstimator::estimate_3d_circles_with_fixed_radius_method1(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp, double desired_radius) {

    estimate_unit_3d_circles_method1(ellipse, circles, ppp);
    for(int i = 0; i < 2; ++i) {
        circles[i].center *= desired_radius;
        circles[i].radius = desired_radius;
        if(circles[i].center[2] > -ppp->near) {
            circles[i].center = - circles[i].center;
            circles[i].normal = - circles[i].normal;
        }
    }
}
*/

// Method-2: Based on analytical solution
int CircleEstimator::estimate_3d_circles_with_fixed_depth(const Ellipse2D& ellipse, Circle3D* circles, const PersProjParam* const ppp, double desired_depth) {

    if(desired_depth > -ppp->near || desired_depth < -ppp->far)
        std::cout << "ERROR: desired depth must be within [-near, -far]: [" << -ppp->near << ", " << -ppp->far << "]" << std::endl;

    int count = estimate_unit_3d_circles(ellipse, circles, ppp);
    for(int i = 0; i < count; ++i) {
        if(desired_depth != circles[i].center[2])
            circles[i].radius = desired_depth/circles[i].center[2];
        else
            continue;
        circles[i].center *= circles[i].radius;
    }

    return count;
}

int CircleEstimator::estimate_3d_circles_with_fixed_radius(const Ellipse2D& ellipse, Circle3D* circles, const PersProjParam* const ppp, double desired_radius) {

    int count = estimate_unit_3d_circles(ellipse, circles, ppp);
    for(int i = 0; i < count; ++i) {
        circles[i].center *= desired_radius;
        circles[i].radius = desired_radius;
    }
    return count;
}


void CircleEstimator::estimate_3d_circles_under_orthographic_projection(const Ellipse2D& ellipse, Circle3D* circles, PersProjParam const * const ppp) {

    osg::Vec2d u2 = ellipse.points[1] - ellipse.center;
    osg::Vec2d v2 = ellipse.points[2] - ellipse.center;

    osg::Vec3 u3(u2.x(), u2.y(), 0);
    osg::Vec3 v3(v2.x(), v2.y(), std::sqrt(u2*u2 - v2*v2));
    osg::Vec3 n = u3 ^ v3;
    circles[0].normal[0] = n.x();
    circles[0].normal[1] = n.y();
    circles[0].normal[2] = n.z();
    circles[0].center[0] = ellipse.center.x();
    circles[0].center[1] = ellipse.center.y();
    circles[0].center[2] = ppp->near;

    circles[1].normal[0] = n.x();
    circles[1].normal[1] = n.y();
    circles[1].normal[2] = n.z();
    circles[1].center[0] = ellipse.center.x();
    circles[1].center[1] = ellipse.center.y();
    circles[1].center[2] = ppp->near;
}


// PRIVATE METHODS
/*
void CircleEstimator::estimate_unit_3d_circles_method1(const Ellipse2D& ellipse, Circle3D* circles, const PersProjParam* const ppp) {

    // we need to negate the near value, because in opengl near and far values are positive. We need the actual value
    // in mathemetical calculations.

    double near = -ppp->near;
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

int CircleEstimator::estimate_unit_3d_circles(const Ellipse2D& ellipse, Circle3D* circles, const PersProjParam *const ppp) {

     // Step-1: Construct the associated quadratic form matrix of the 3D cone.
     /*
     * 'ellipse' is the intersection of the 3D cone (whose vertex is at origin) and a plane z = k.
     * The plane is the near clipping plane, since the ellipse is on the near clipping plane.
     * k = -near (recall that near > 0).
     *
     * The cone is constructed in XYZ coordinate frame.
     */

    double near = -ppp->near;
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
