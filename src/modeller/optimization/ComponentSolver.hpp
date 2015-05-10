#ifndef COMPONENT_SOLVER_HPP
#define COMPONENT_SOLVER_HPP

#include <osg/Array>
#include <ceres/ceres.h>

struct CostFunctor {

    CostFunctor(double* c1_, double* c2_) : c1(c1_), c2(c2_) { }

    template <typename T>
    bool operator()(const T* const Z0, const T* const C1, const T* const C2, T* residual) const {
        residual[0] = (T(c1[0]) * (*Z0) - C1[0]) * (T(c1[0]) * (*Z0) - C1[0]) +
                      (T(c1[1]) * (*Z0) - C1[1]) * (T(c1[1]) * (*Z0) - C1[1]) +
                      (T(c1[2]) * (*Z0) - C1[2]) * (T(c1[2]) * (*Z0) - C1[2]) +
                      (T(c2[0]) * (*Z0) - C2[0]) * (T(c2[0]) * (*Z0) - C2[0]) +
                      (T(c2[1]) * (*Z0) - C2[1]) * (T(c2[1]) * (*Z0) - C2[1]) +
                      (T(c2[2]) * (*Z0) - C2[2]) * (T(c2[2]) * (*Z0) - C2[2]);
        return true;
    }
    double* c1;
    double* c2;
};

class ComponentSolver {
public:
    ComponentSolver(double near_);
    void Solve(double& Z0, double* C1, double* C2);
private:
    osg::ref_ptr<osg::Vec2dArray> m_lframe_proj;  // projection of the local frame points:
    // m_lframe_proj->at(0) = P0 : origin - third click,
    // m_lframe_proj->at(1) = P1 : first click,
    // m_lframe_proj->at(2) = P1 : second click,
    // m_lframe_proj->at(3) = P3 : fourth click

    double coeff[3];       // Z1 = coeff[0]*Z0, Z2 = coeff[1]*Z0,  Z3 = coeff[3]*Z0
    double C1v[3];         // C1 = Z0*C1v
    double C2v[3];         // C2 = Z0*C2v
    double near;
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    inline double func(int i, int j);
    void calculate_coefficients();
};

#endif // COMPONENT_SOLVER_HPP
