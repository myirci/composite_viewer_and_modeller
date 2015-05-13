#include <iostream>
#include <cmath>
#include "ComponentSolver.hpp"


ComponentSolver::ComponentSolver(double near_) : m_lframe_proj(new osg::Vec2dArray(4)), near(near_) {

    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
}

double ComponentSolver::func(int i, int j) {
    // n^2 + xi*xj + yi*yj
    return near*near +
           m_lframe_proj->at(i).x() * m_lframe_proj->at(j).x() +
           m_lframe_proj->at(i).y() * m_lframe_proj->at(j).y();
}

void ComponentSolver::calculate_coefficients() {

    double f00 = func(0,0);
    double f01 = func(0,1);
    double f02 = func(0,2);
    double f03 = func(0,3);
    double f12 = func(1,2);
    double f13 = func(1,3);
    double f23 = func(2,3);

    std::cout.precision(16);
    // a*Z1^2 + b*Z1*Z0 + c*Z0^2 = 0
    double a = f01*f02*f13 + f01*f03*f12 - f00*f12*f13 - f01*f01*f23;
    double b = 2*f01*(f00*f23 - f02*f03);
    double c = f00*(f02*f03 - f00*f23);

    // Z1 = k1*Z0
    double dlt = b*b - 4*a*c;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "discriminant: " << dlt << std::endl;
    double k1_1 = (-b + std::sqrt(dlt)) / (2*a);
    double k1_2 = (-b - std::sqrt(dlt)) / (2*a);
    std::cout << "k11: " << k1_1 << std::endl;
    std::cout << "k12: " << k1_2 << std::endl;
    coeff[0] = k1_1; // or k1_2
    coeff[1] = (coeff[0]*f01 - f00) / (coeff[0]*f12 - f02);
    coeff[2] = (coeff[0]*f01 - f00) / (coeff[0]*f13 - f03);

    // C1 = (P1 + P2)/2
    C1v[0] = (coeff[0] * m_lframe_proj->at(1).x() + coeff[1] * m_lframe_proj->at(2).x()) / 2*near;
    C1v[1] = (coeff[0] * m_lframe_proj->at(1).y() + coeff[1] * m_lframe_proj->at(2).y()) / 2*near;
    C1v[2] = (coeff[0] + coeff[1]) / 2.0;

    // C2 = C1 + P3 - P0
    C2v[0] = C1v[0] + (coeff[2] * m_lframe_proj->at(3).x() - m_lframe_proj->at(0).x()) / near;
    C2v[1] = C1v[1] + (coeff[2] * m_lframe_proj->at(3).y() - m_lframe_proj->at(0).y()) / near;
    C2v[2] = C1v[2] + coeff[2] - 1;
}

osg::Vec2dArray* ComponentSolver::GetLocalCoordinateFrameProjections() {

    return m_lframe_proj.get();
}

void ComponentSolver::Solve(double& Z0, double* C1, double* C2) {

    // print the initial values:
    std::cout << "Initial Z0: " << Z0 << std::endl;
    std::cout << "Initial C1: " << C1[0] << " " << C1[1] << " " << C1[2] << std::endl;
    std::cout << "Initial C2: " << C2[0] << " " << C2[1] << " " << C2[2] << std::endl;

    calculate_coefficients();
    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 2, 1, 3, 3>(new CostFunctor(C1v, C2v));
    problem.AddResidualBlock(cost_function, NULL, &Z0, C1, C2);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    // print the final values
    std::cout << "Final Z0: " << Z0 << std::endl;
    std::cout << "Final C1: " << C1[0] << " " << C1[1] << " " << C1[2] << std::endl;
    std::cout << "Final C2: " << C2[0] << " " << C2[1] << " " << C2[2] << std::endl;

    P0 = osg::Vec3d(Z0 * m_lframe_proj->at(0).x() / near, Z0 * m_lframe_proj->at(0).y() / near, Z0);
    P1 = osg::Vec3d(coeff[0] * Z0 * m_lframe_proj->at(1).x() / near, coeff[0] * Z0 * m_lframe_proj->at(1).y() / near, coeff[0] * Z0);
    P2 = osg::Vec3d(coeff[1] * Z0 * m_lframe_proj->at(2).x() / near, coeff[1] * Z0 * m_lframe_proj->at(2).y() / near, coeff[1] * Z0);
    P3 = osg::Vec3d(coeff[2] * Z0 * m_lframe_proj->at(3).x() / near, coeff[2] * Z0 * m_lframe_proj->at(3).y() / near, coeff[2] * Z0);

    std::cout << "P0: " << P0.x() << " " << P0.y() << " " << P0.z() << std::endl;
    std::cout << "P1: " << P1.x() << " " << P1.y() << " " << P1.z() << std::endl;
    std::cout << "P2: " << P2.x() << " " << P2.y() << " " << P2.z() << std::endl;
    std::cout << "P3: " << P3.x() << " " << P2.y() << " " << P3.z() << std::endl;

    osg::Vec3d calculated_c1 = (P1 + P2) / 2.0;
    std::cout << "Calculated_C1: " << calculated_c1.x() << " " << calculated_c1.y() << " " << calculated_c1.z() << std::endl;
    std::cout << "Distance Calculated_C1 to P0: " << (P0-calculated_c1).length() << std::endl;
    std::cout << "Distance Calculated_C1 to P1: " << (P1-calculated_c1).length() << std::endl;
    std::cout << "Distance Calculated_C1 to P2: " << (P2-calculated_c1).length() << std::endl;

    osg::Vec3d optimized_c1(C1[0],C1[1],C1[2]);
    std::cout << "Distance optimized_c1 to P0: " << (P0-optimized_c1).length() << std::endl;
    std::cout << "Distance optimized_c1 to P1: " << (P1-optimized_c1).length() << std::endl;
    std::cout << "Distance optimized_c1 to P2: " << (P2-optimized_c1).length() << std::endl;

    std::cout << "orthogonality constraint-1: " << (P1-P0) * (P2-P0) << std::endl;
    std::cout << "orthogonality constraint-2: " << (P2-P0) * (P3-P0) << std::endl;
    std::cout << "orthogonality constraint-3: " << (P1-P0) * (P3-P0) << std::endl;
}
