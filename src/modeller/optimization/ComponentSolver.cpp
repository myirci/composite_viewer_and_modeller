#include <iostream>
#include "ComponentSolver.hpp"
#include "../../geometry/Plane3D.hpp"
#include "../../geometry/Ray3D.hpp"
#include "../../utility/Utility.hpp"
#include "../components/GeneralizedCylinderGeometry.hpp"
#include "../components/GeneralizedCylinder.hpp"

void ComponentSolver::SolveForSingleCircle(osg::Vec2dArray const * const proj, Circle3D& circle) {

    // initialize the optimization parameters
    double center[3] = { circle.center[0], circle.center[1], circle.center[2] };
    double depths[2];
    Plane3D cplane(osg::Vec3d(circle.normal[0], circle.normal[1], circle.normal[2]),
                   osg::Vec3d(circle.center[0], circle.center[1], circle.center[2]));
    osg::Vec3d p;
    intersect_ray_and_plane(Ray3D(osg::Vec3d(proj->at(0), n)), cplane, p);
    depths[0] = p.z();
    intersect_ray_and_plane(Ray3D(osg::Vec3d(proj->at(1), n)), cplane, p);
    depths[1] = p.z();

    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor_1, 4, 3, 2>(new CostFunctor_1(proj, circle, n));
    problem.AddResidualBlock(cost_function, NULL, center, depths);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    // update the circle
    circle.center[0] = center[0];
    circle.center[1] = center[1];
    circle.center[2] = center[2];
    osg::Vec3d P1(depths[0]*proj->at(0).x()/n, depths[0]*proj->at(0).y()/n, depths[0]);
    osg::Vec3d P2(depths[1]*proj->at(1).x()/n, depths[1]*proj->at(1).y()/n, depths[1]);
    osg::Vec3d C(center[0], center[1], center[2]);
    osg::Vec3d nrm = ((C - P1) ^ (C - P2));
    nrm.normalize();
    circle.normal[0] = nrm.x();
    circle.normal[1] = nrm.y();
    circle.normal[2] = nrm.z();
}

void ComponentSolver::SolveGeneralizedCylinder(GeneralizedCylinder* gcyl) {

    ceres::Problem problem;
    // initialize the lambdas

    GeneralizedCylinderGeometry* geom = gcyl->GetGeometry();
    std::vector<double> lambdas(geom->GetNumberOfSections() - 1, 1.0);

    for(int i = 0; i < lambdas.size(); ++i)
        std::cout << lambdas[i] << std::endl;
    std::cout << "----------------" << std::endl;

    ceres::DynamicAutoDiffCostFunction<CostFunctor_2, 4>* cost_function =
            new ceres::DynamicAutoDiffCostFunction<CostFunctor_2, 4>(new CostFunctor_2(geom->GetSections()));

    std::vector<double*> parameters;
    for(int i = 0; i < lambdas.size(); ++i) {
        parameters.push_back(&lambdas[i]);
        cost_function->AddParameterBlock(1);
    }
    cost_function->SetNumResiduals(geom->GetNumberOfSections() - 1);

    problem.AddResidualBlock(cost_function, NULL, parameters);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    for(int i = 0; i < lambdas.size(); ++i)
        std::cout << lambdas[i] << std::endl;

    // update the generalized cylinder
    std::vector<Circle3D>& sections = geom->GetSections();
    for(int i = 1; i < sections.size(); ++i) {
        sections[i].center *= lambdas[i-1];
        sections[i].radius *= lambdas[i-1];
    }
    gcyl->Recalculate();
}

void ComponentSolver::SolveDepth(const Circle3D& C0, Circle3D& C1) {

    // initialize the optimization parameters
    double s = 1;
    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor_Depth, 1, 1>(new CostFunctor_Depth(C0, C1));
    problem.AddResidualBlock(cost_function, NULL, &s);
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    // scale the circle
    C1.center *= s;
    C1.radius *= s;
}



