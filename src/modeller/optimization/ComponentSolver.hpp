#ifndef COMPONENT_SOLVER_HPP
#define COMPONENT_SOLVER_HPP

#include "OptimizationUtility.hpp"
#include "../../geometry/Circle3D.hpp"
#include <osg/Array>
#include <ceres/ceres.h>

class GeneralizedCylinder;

struct CostFunctor_1 {

    CostFunctor_1(osg::Vec2dArray const * const proj, Circle3D& circle, double near) :
        m_proj(proj), m_circle(circle), n(near) { }

    template <typename T>
    bool operator()(const T* const ctr, const T* const depths, T* residuals) const {

        Vector3D<T> C(ctr[0], ctr[1], ctr[2]);
        Vector3D<T> P1(depths[0] * T(m_proj->at(0).x()), depths[0] * T(m_proj->at(0).y()), depths[0]);
        Vector3D<T> P2(depths[1] * T(m_proj->at(1).x()), depths[1] * T(m_proj->at(1).y()), depths[1]);
        Vector3D<T> N = (C-P1).cross(C-P2);
        Vector3D<T> Cp(T(m_circle.center[0]), T(m_circle.center[1]), T(m_circle.center[2]));
        T zc = T(2)*depths[0]*depths[1] / (depths[0] + depths[1]);
        Vector3D<T> Pc(zc * T((m_proj->at(0).x() + m_proj->at(1).x())/2), T((m_proj->at(0).y() + m_proj->at(1).y())/2), zc);

        T r = T(m_circle.radius);

        Vector3D<T> vec1(T(m_proj->at(0).x()), T(m_proj->at(0).x()), T(n));
        Vector3D<T> vec2(T(m_proj->at(1).x()), T(m_proj->at(1).x()), T(n));

        // residual for bisector plane
        // residuals[0] = (P1-P2).dot(C - (P1+P2)/T(2.0));
        residuals[0] = (T(m_proj->at(1).x()) * depths[1] - T(m_proj->at(0).x()) * depths[0]) * ctr[0] +
                       (T(m_proj->at(1).y()) * depths[1] - T(m_proj->at(0).y()) * depths[0]) * ctr[1] +
                       T(n)*(depths[1] - depths[0]) * ctr[2] +
                       (depths[0]*depths[0]*vec1.dot(vec1) - depths[1]*depths[1]*vec2.dot(vec2)) / T(2*n);

        // residual for normal alignment
        residuals[1] = N.cross(Cp - C).squared_norm();
        residuals[2] = r*r - (C - P1).squared_norm();
        residuals[3] = r*r - (C - P2).squared_norm();

        return true;
    }
    private:
    osg::Vec2dArray const * const m_proj;
    const Circle3D& m_circle;
    double n;
};

struct CostFunctor_2 {

    CostFunctor_2(const std::vector<Circle3D>& circular_sections) : sections(circular_sections) { }

    template <typename T>
    Vector3D<T> get_normal(size_t index) const {
        return Vector3D<T>(T(sections[index].normal[0]), T(sections[index].normal[1]), T(sections[index].normal[2]));
    }

    template <typename T>
    Vector3D<T> get_center(size_t index) const {
        return Vector3D<T>(T(sections[index].center[0]), T(sections[index].center[1]), T(sections[index].center[2]));
    }

    template <typename T>
    bool operator()(T const* const* parameters, T* residuals) const {

        // residual-0
        residuals[0] = (get_normal<T>(1).cross(get_center<T>(0) - parameters[0][0]*get_center<T>(1))).norm();

        // rest of the residuals:
        for(int i = 2; i < sections.size(); ++i)
            residuals[i-1] = (get_normal<T>(i).cross(parameters[0][i-2]*get_center<T>(i-1) - parameters[0][i-1]*get_center<T>(i))).norm();

        return true;
    }
    private:
    const std::vector<Circle3D>& sections;
};

struct CostFunctor_Depth {

    CostFunctor_Depth(const Circle3D& c0, const Circle3D& c1) : m_circle_0(c0), m_circle_1(c1) { }
    template <typename T>
    bool operator()(const T* const s, T* residual) const {
        Vector3D<T> C0(T(m_circle_0.center[0]), T(m_circle_0.center[1]), T(m_circle_0.center[2]));
        Vector3D<T> n0(T(m_circle_0.normal[0]), T(m_circle_0.normal[1]), T(m_circle_0.normal[2]));
        Vector3D<T> C1(T(m_circle_1.center[0]), T(m_circle_1.center[1]), T(m_circle_1.center[2]));
        // Vector3D<T> n1(T(m_circle_1.normal[0]), T(m_circle_1.normal[1]), T(m_circle_1.normal[2]));
        C1 = C1 * s[0];
        Vector3D<T> diff = C1 - C0;
        // residual[0] = (diff.cross(n0)).squared_norm() + (diff.cross(n1)).squared_norm();
        residual[0] = (diff.cross(n0)).squared_norm();
        return true;
    }
    private:
    const Circle3D& m_circle_0;
    const Circle3D& m_circle_1;
};


class ComponentSolver {
public:
    ComponentSolver(double near) : n(near) {
        options.max_num_iterations = 100;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
    }
    void SolveForSingleCircle(osg::Vec2dArray const * const proj, Circle3D& circle);
    void SolveGeneralizedCylinder(GeneralizedCylinder* gcyl);
    void SolveDepth(const Circle3D& C0, Circle3D& C1);
private:
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    double n;
};

#endif // COMPONENT_SOLVER_HPP
