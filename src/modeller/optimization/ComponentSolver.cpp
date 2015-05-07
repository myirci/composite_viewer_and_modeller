#include <iostream>
#include <cmath>
#include "ComponentSolver.hpp"


ComponentSolver::ComponentSolver(double near_) : m_lframe_proj(new osg::Vec2dArray(4)), near(near_), squared_near(near_ * near_) {

    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
}

double ComponentSolver::func(int i, int j) {
    return squared_near + m_lframe_proj->at(i).x() * m_lframe_proj->at(j).x() + m_lframe_proj->at(i).y() * m_lframe_proj->at(j).y();
}

void ComponentSolver::calculate_coefficients() {

    double f00 = func(0,0);
    double f01 = func(0,1);
    double f02 = func(0,2);
    double f03 = func(0,3);
    double f12 = func(1,2);
    double f13 = func(1,3);
    double f23 = func(2,3);

    double a = f00*f12*f13 + f01*f01*f23 - f01*f02*f13 - f01*f03*f12;
    double b = 2*f01*(f02*f03 - f00*f23);
    double c = f00*(f00*f23 - f02*f03);

    double dlt = b*b - 4*a*c;
    std::cout << "discriminant: " << dlt << std::endl;
    double d1 = (-b + std::sqrt(dlt)) / (2*a);
    double d2 = (-b - std::sqrt(dlt)) / (2*a);
    std::cout << "d1: " << d1 << std::endl;
    std::cout << "d2: " << d2 << std::endl;
    coeff[0] = d1; // or d2
    coeff[2] = (coeff[0]*f01 - f00) / (coeff[0]*f13 - f03);
    coeff[1] = (coeff[2]*f03 - f00) / (coeff[2]*f23 - f02);

    C1v[0] = (coeff[0] * m_lframe_proj->at(1).x() + coeff[1] * m_lframe_proj->at(2).x()) / 2*near;
    C1v[1] = (coeff[0] * m_lframe_proj->at(1).y() + coeff[1] * m_lframe_proj->at(2).y()) / 2*near;
    C1v[2] = (coeff[0] + coeff[1]) / 2.0;

    C2v[0] = C1v[0] + (coeff[2] * m_lframe_proj->at(3).x() - m_lframe_proj->at(0).x()) / near;
    C2v[1] = C1v[1] + (coeff[2] * m_lframe_proj->at(3).y() - m_lframe_proj->at(0).y()) / near;
    C2v[2] = C1v[2] + coeff[2] - 1;
}

void ComponentSolver::Solve(double& Z0, double* C1, double* C2) {

    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 1, 1, 3, 3>(new CostFunctor(C1v, C2v));
    problem.AddResidualBlock(cost_function, NULL, &Z0, C1, C2);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
}
