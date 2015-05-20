#ifndef COMPONENT_SOLVER_HPP
#define COMPONENT_SOLVER_HPP

#include "OptimizationUtility.hpp"
#include <osg/Array>
#include <ceres/ceres.h>

struct CostFunctor {
    CostFunctor(osg::Vec2dArray const * const projections, const osg::Vec3d& previous_ctr, double near) :
        m_projections(new osg::Vec2dArray(projections->size())), m_center(previous_ctr)  {

        for(int i = 0; i < projections->size(); ++i)
            m_projections->at(i) = projections->at(i) / near;
    }

    template <typename T>
    bool operator()(const T* const ctr, const T* const depths, T* residuals) const {

        // depths[0]: depth of Pi1, depths[1]: depth of Ci', depths[2]: depth of Pi2
        Vector3D<T> Cvec(ctr[0], ctr[1], ctr[2]);
        Vector3D<T> P1vec(depths[0] * T(m_projections->at(0).x()), depths[0] * T(m_projections->at(0).y()), depths[0]);
        Vector3D<T> Cpvec(depths[1] * T(m_projections->at(1).x()), depths[1] * T(m_projections->at(1).y()), depths[1]);
        Vector3D<T> P2vec(depths[2] * T(m_projections->at(2).x()), depths[2] * T(m_projections->at(2).y()), depths[2]);
        residuals[0] = (Cvec - Cpvec).squared_norm() + (P1vec - Cpvec).norm() * (P2vec - Cpvec).norm() - (P1vec - Cvec).norm() * (P2vec - Cvec).norm();

        Vector3D<T> prev_ctr(T(m_center.x()), T(m_center.y()), T(m_center.z()));
        residuals[1] = ((P1vec - Cvec).cross(P2vec - Cvec)).cross(Cvec - prev_ctr).norm();

        return true;
    }

    private:
    // projections: m_projections->at(0): pi1, m_projections->at(1): ellipse center, m_projections->at(2): pi2
    osg::ref_ptr<osg::Vec2dArray> m_projections;
    // center of previous circle
    osg::Vec3d m_center;
};

class ComponentSolver {
public:
    ComponentSolver(double near) : n(near) {
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = true;
    }
    void SolveForSingleCircle(osg::Vec2dArray const * const projections, const osg::Vec3d& ctr_prev, osg::Vec3d& ctr, osg::Vec3d& depths);
    void SolveForAllCircles();
private:
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    double n;
};

#endif // COMPONENT_SOLVER_HPP
