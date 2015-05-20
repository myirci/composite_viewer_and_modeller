#include <iostream>
#include "ComponentSolver.hpp"

void ComponentSolver::SolveForSingleCircle(osg::Vec2dArray const * const projections, const osg::Vec3d& ctr_prev, osg::Vec3d& ctr, osg::Vec3d& depths) {

    double ctrArr[3]   = { ctr.x(), ctr.y(), ctr.z()          };
    double depthArr[3] = { depths.x(), depths.y(), depths.z() };

    // print the initial values:
    std::cout << "Initial Center " << ctrArr[0] << " " << ctrArr[1] << " " << ctrArr[2] << std::endl;
    std::cout << "Initial depth of P1: " << depthArr[0] << std::endl;
    std::cout << "Initial depth of Cp: " << depthArr[1] << std::endl;
    std::cout << "Initial depth of P2: " << depthArr[2] << std::endl;

    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 2, 3, 3>(new CostFunctor(projections, ctr_prev, n));
    problem.AddResidualBlock(cost_function, NULL, ctrArr, depthArr);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    // print the final values
    std::cout << "Final Center " << ctrArr[0] << " " << ctrArr[1] << " " << ctrArr[2] << std::endl;
    std::cout << "Final depth of P1: " << depthArr[0] << std::endl;
    std::cout << "Final depth of Cp: " << depthArr[1] << std::endl;
    std::cout << "Final depth of P2: " << depthArr[2] << std::endl;

    ctr.x() = ctrArr[0];
    ctr.y() = ctrArr[1];
    ctr.z() = ctrArr[2];
    depths.x() = depthArr[0];
    depths.y() = depthArr[1];
    depths.z() = depthArr[2];
}
