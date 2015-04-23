#include "ComponentSolver.hpp"
#include "../../geometry/Circle3D.hpp"
#include "../../utility/Utility.hpp"
#include "../CoordinateTransformations.hpp"
#include <iostream>

ComponentSolver::ComponentSolver(const std::shared_ptr<CoordinateTransformations>& ppp) : m_lfpts(nullptr), m_ppp(ppp) {

    base_circles = new Circle3D[2];
    base_circles[0].radius = 1.0;
    base_circles[1].radius = 1.0;
}

ComponentSolver::~ComponentSolver() {
    delete[] base_circles;
}

Circle3D* ComponentSolver::GetBaseCircles() {
    return base_circles;
}

void ComponentSolver::SolveDepthValues() const {

    if(!m_lfpts.valid()) {
        std::cout << "Coordinate frame points have not been updated" << std::endl;
        return;
    }

    // For base_circle_1:
    double k = std::sqrt(-F(1, 2) * F(1, 3) * F(2, 3));
    double a1 = (1/F(2, 3)) * k; // z1 = z0 +- a1
    double a2 = (1/F(1, 3)) * k; // z2 = z0 +- a2
    double a3 = (1/F(1, 2)) * k; // z3 = z0 +- a3

    osg::Vec3d cbt, ctp, P0, P1, P2, P3;

    auto z01 = quadratic_root_finder_x(base_circles[0].center[0], a1, a2);
    calculate_P0(P0, z01.first);
    calculate_P1(P1, z01.first + a1);
    calculate_P2(P2, z01.first + a2);
    calculate_P3(P3, z01.first + a3);
    calculate_bottom_center_point(cbt, z01.first+a1, z01.first+a2);
    calculate_top_center_point(cbt, ctp, z01.first, z01.first+a3);

    calculate_P0(P0, z01.second);
    calculate_P1(P1, z01.second + a1);
    calculate_P2(P2, z01.second + a2);
    calculate_P3(P3, z01.second + a3);
    calculate_bottom_center_point(cbt, z01.second+a1, z01.second+a2);
    calculate_top_center_point(cbt, ctp, z01.second, z01.second+a3);

    auto z04 = quadratic_root_finder_x(base_circles[0].center[0], -a1, -a2);
    calculate_P0(P0, z04.first);
    calculate_P1(P1, z04.first - a1);
    calculate_P2(P2, z04.first - a2);
    calculate_P3(P3, z04.first - a3);
    calculate_bottom_center_point(cbt, z04.first-a1, z04.first-a2);
    calculate_top_center_point(cbt, ctp, z04.first, z04.first-a3);

    calculate_P0(P0, z04.second);
    calculate_P1(P1, z04.second - a1);
    calculate_P2(P2, z04.second - a2);
    calculate_P3(P3, z04.second - a3);
    calculate_bottom_center_point(cbt, z04.second-a1, z04.second-a2);
    calculate_top_center_point(cbt, ctp, z04.second, z04.second-a3);

    auto z05 = quadratic_root_finder_x(base_circles[1].center[0], a1, a2);
    calculate_P0(P0, z05.first);
    calculate_P1(P1, z05.first + a1);
    calculate_P2(P2, z05.first + a2);
    calculate_P3(P3, z05.first + a3);
    calculate_bottom_center_point(cbt, z05.first+a1, z05.first+a2);
    calculate_top_center_point(cbt, ctp, z05.first, z05.first+a3);

    calculate_P0(P0, z05.second);
    calculate_P1(P1, z05.second + a1);
    calculate_P2(P2, z05.second + a2);
    calculate_P3(P3, z05.second + a3);
    calculate_bottom_center_point(cbt, z05.second+a1, z05.second+a2);
    calculate_top_center_point(cbt, ctp, z05.second, z05.second+a3);

    auto z08 = quadratic_root_finder_x(base_circles[1].center[0], -a1, -a2);
    calculate_P0(P0, z08.first);
    calculate_P1(P1, z08.first - a1);
    calculate_P2(P2, z08.first - a2);
    calculate_P3(P3, z08.first - a3);
    calculate_bottom_center_point(cbt, z08.first-a1, z08.first-a2);
    calculate_top_center_point(cbt, ctp, z08.first, z08.first-a3);

    calculate_P0(P0, z08.second);
    calculate_P1(P1, z08.second - a1);
    calculate_P2(P2, z08.second - a2);
    calculate_P3(P3, z08.second - a3);
    calculate_bottom_center_point(cbt, z08.second-a1, z08.second-a2);
    calculate_top_center_point(cbt, ctp, z08.second, z08.second-a3);

 /*
    auto z01 = quadratic_root_finder_y(base_circles[0].center[1], a1, a2);
    calculate_bottom_center_point(cbt, z01.first+a1, z01.first+a2);
    calculate_bottom_center_point(cbt, z01.second+a1, z01.second+a2);
    auto z04 = quadratic_root_finder_y(base_circles[0].center[1], -a1, -a2);
    calculate_bottom_center_point(cbt, z04.first-a1, z04.first-a2);
    calculate_bottom_center_point(cbt, z04.second-a1, z04.second-a2);
    auto z05 = quadratic_root_finder_y(base_circles[1].center[1], a1, a2);
    calculate_bottom_center_point(cbt, z05.first+a1, z05.first+a2);
    calculate_bottom_center_point(cbt, z05.second+a1, z05.second+a2);
    auto z08 = quadratic_root_finder_y(base_circles[1].center[1], -a1, -a2);
    calculate_bottom_center_point(cbt, z08.first-a1, z08.first-a2);
    calculate_bottom_center_point(cbt, z08.second-a1, z08.second-a2);

*/
    /*
    auto z01 = quadratic_root_finder_z(base_circles[0].center[2], a1, a2);
    calculate_bottom_center_point(cbt, z01.first+a1, z01.first+a2);
    calculate_top_center_point(cbt, ctp, z01.first+a1, z01.first+a3);
    calculate_bottom_center_point(cbt, z01.second+a1, z01.second+a2);
    calculate_top_center_point(cbt, ctp, z01.second+a1, z01.second+a3);
    auto z04 = quadratic_root_finder_z(base_circles[0].center[2], -a1, -a2);
    calculate_bottom_center_point(cbt, z04.first-a1, z04.first-a2);
    calculate_top_center_point(cbt, ctp, z04.first-a1, z04.first-a3);
    calculate_bottom_center_point(cbt, z04.second-a1, z04.second-a2);
    calculate_top_center_point(cbt, ctp, z04.second-a1, z04.second-a3);
    auto z05 = quadratic_root_finder_z(base_circles[1].center[2], a1, a2);
    calculate_bottom_center_point(cbt, z05.first+a1, z05.first+a2);
    calculate_top_center_point(cbt, ctp, z05.first+a1, z05.first+a3);
    calculate_bottom_center_point(cbt, z05.second+a1, z05.second+a2);
    calculate_top_center_point(cbt, ctp, z05.second+a1, z05.second+a3);
    auto z08 = quadratic_root_finder_z(base_circles[1].center[2], -a1, -a2);
    calculate_bottom_center_point(cbt, z08.first-a1, z08.first-a2);
    calculate_top_center_point(cbt, ctp, z08.first-a1, z01.first-a3);
    calculate_bottom_center_point(cbt, z08.second-a1, z08.second-a2);
    calculate_top_center_point(cbt, ctp, z08.second-a1, z08.second-a3);
    */
}

double ComponentSolver::F(int i, int j) const {

    return m_ppp->far * m_ppp->far *
           ((m_lfpts->at(0).x()*m_lfpts->at(0).x() + m_lfpts->at(0).y()*m_lfpts->at(0).y()) - // E00
            (m_lfpts->at(0).x()*m_lfpts->at(i).x() + m_lfpts->at(0).y()*m_lfpts->at(i).y()) - // E0i
            (m_lfpts->at(0).x()*m_lfpts->at(j).x() + m_lfpts->at(0).y()*m_lfpts->at(j).y()) + // E0j
            (m_lfpts->at(i).x()*m_lfpts->at(j).x() + m_lfpts->at(i).y()*m_lfpts->at(j).y())   // Eij
           );
}

std::pair<double, double> ComponentSolver::quadratic_root_finder_x(double x, double a1, double a2) const {

    double g = m_ppp->far + m_ppp->near;
    double a = 2*x;
    double b = a*(a1+a2+2*g) - m_ppp->far*(m_lfpts->at(1).x() + m_lfpts->at(2).x());
    double c = a*(a1*a2+g*(a1+a2+g)) - m_ppp->far*(m_lfpts->at(1).x()*(a2+g) + m_lfpts->at(2).x()*(a1+g));
    double sqrt_delta = std::sqrt(b*b - 4*a*c);
    return std::make_pair((-b + sqrt_delta) / (2*a), (-b - sqrt_delta) / (2*a));
}

std::pair<double, double> ComponentSolver::quadratic_root_finder_y(double y, double a1, double a2) const {

    double g = m_ppp->far + m_ppp->near;
    double a = 2*y;
    double b = a*(a1+a2+2*g) - m_ppp->far*(m_lfpts->at(1).y() + m_lfpts->at(2).y());
    double c = a*(a1*a2+g*(a1+a2+g)) - m_ppp->far*(m_lfpts->at(1).y()*(a2+g) + m_lfpts->at(2).y()*(a1+g));
    double sqrt_delta = std::sqrt(b*b - 4*a*c);
    return std::make_pair((-b + sqrt_delta) / (2*a), (-b - sqrt_delta) / (2*a));
}

std::pair<double, double> ComponentSolver::quadratic_root_finder_z(double z, double a1, double a2) const {

    double g = m_ppp->far + m_ppp->near;
    double a = 2*z;
    double b = a*(a1+a2+2*g) + 2*m_ppp->far*m_ppp->near;
    double c = a*(a1*a2+g*(a1+a2+g)) + m_ppp->far*m_ppp->near*(a1+a2+2*g);
    double sqrt_delta = std::sqrt(b*b - 4*a*c);
    return std::make_pair((-b + sqrt_delta) / (2*a), (-b - sqrt_delta) / (2*a));
}

void ComponentSolver::calculate_bottom_center_point(osg::Vec3d& cbt, double z1, double z2) const {

    double g = m_ppp->far + m_ppp->near;
    cbt.set(0.5*m_ppp->far*(m_lfpts->at(1).x()/(z1+g) + m_lfpts->at(2).x()/(z2 + g)),
            0.5*m_ppp->far*(m_lfpts->at(1).y()/(z1+g) + m_lfpts->at(2).y()/(z2 + g)),
           -0.5*m_ppp->far*m_ppp->near*(1/(z1+g) + 1/(z2+g)));
}

void ComponentSolver::calculate_top_center_point(const osg::Vec3d& cbt, osg::Vec3d& ctp, double z0, double z3) const {

    double g = m_ppp->far + m_ppp->near;
    ctp = cbt + osg::Vec3d(((m_ppp->far*m_lfpts->at(3).x())/(z3+g)) - ((m_ppp->far*m_lfpts->at(0).x())/(z0+g)),
                           ((m_ppp->far*m_lfpts->at(3).y())/(z3+g)) - ((m_ppp->far*m_lfpts->at(0).y())/(z0+g)),
                           ((m_ppp->far*m_ppp->near)/(z0+g)) - ((m_ppp->far*m_ppp->near)/(z3+g)));
}

void ComponentSolver::calculate_P0(osg::Vec3d& P0, double z0) const {

    P0[0] = (m_ppp->far  * m_lfpts->at(0).x()) / (z0 + m_ppp->far + m_ppp->near);
    P0[1] = (m_ppp->far  * m_lfpts->at(0).y()) / (z0 + m_ppp->far + m_ppp->near);
    P0[2] = (-m_ppp->far * m_ppp->near)        / (z0 + m_ppp->far + m_ppp->near);
}

void ComponentSolver::calculate_P1(osg::Vec3d& P1, double z1) const {

    P1[0] = (m_ppp->far  * m_lfpts->at(1).x()) / (z1 + m_ppp->far + m_ppp->near);
    P1[1] = (m_ppp->far  * m_lfpts->at(1).y()) / (z1 + m_ppp->far + m_ppp->near);
    P1[2] = (-m_ppp->far * m_ppp->near)        / (z1 + m_ppp->far + m_ppp->near);
}

void ComponentSolver::calculate_P2(osg::Vec3d& P2, double z2) const {

    P2[0] = (m_ppp->far*m_lfpts->at(2).x())/(z2+m_ppp->far+m_ppp->near);
    P2[1] = (m_ppp->far*m_lfpts->at(2).y())/(z2+m_ppp->far+m_ppp->near);
    P2[2] = (-m_ppp->far*m_ppp->near)/(z2+m_ppp->far+m_ppp->near);
}

void ComponentSolver::calculate_P3(osg::Vec3d& P3, double z3) const {

    P3[0] = (m_ppp->far*m_lfpts->at(3).x())/(z3+m_ppp->far+m_ppp->near);
    P3[1] = (m_ppp->far*m_lfpts->at(3).y())/(z3+m_ppp->far+m_ppp->near);
    P3[2] = (-m_ppp->far*m_ppp->near)/(z3+m_ppp->far+m_ppp->near);
}

double ComponentSolver::test_orthogonality(const osg::Vec3d& P0, const osg::Vec3d& P1, const osg::Vec3d& P2, const osg::Vec3d& P3) const {
    return std::abs((P1 - P0)*(P2 - P0)) + std::abs((P1 - P0)*(P3 - P0)) + std::abs((P2 - P0)*(P3 - P0));
}

double ComponentSolver::test_orthogonality(const osg::Vec3d& cbt, const osg::Vec3d& ctp, int base) const {
    osg::Vec3d normal(base_circles[base].normal[0], base_circles[base].normal[1], base_circles[base].normal[2]);
    osg::Vec3d spine = ctp - cbt;
    double angle = std::acos((spine * normal) /(spine.length() * normal.length()));
    return rad2deg(angle);
}
