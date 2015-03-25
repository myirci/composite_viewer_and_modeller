#ifndef COMPONENT_SOLVER_HPP
#define COMPONENT_SOLVER_HPP

#include <osg/Array>
#include <memory>

class Circle3D;
class PersProjParam;

class ComponentSolver {
public:
    ComponentSolver(const std::shared_ptr<PersProjParam>& ppp);
    ~ComponentSolver();
    Circle3D* GetBaseCircles();
    void SolveDepthValues() const;
private:

    Circle3D* base_circles;                 // possible base circles of a generalized cylinder
    osg::ref_ptr<osg::Vec2dArray> m_lfpts;  // local frame end points
    std::shared_ptr<PersProjParam> m_ppp;   // perspective projection parameters

    void calculate_bottom_center_point(osg::Vec3d& cbt, double z1, double z2) const;
    void calculate_top_center_point(const osg::Vec3d& cbt, osg::Vec3d& ctp, double z0, double z3) const;

    void calculate_P0(osg::Vec3d& P0, double z0) const;
    void calculate_P1(osg::Vec3d& P1, double z1) const;
    void calculate_P2(osg::Vec3d& P2, double z2) const;
    void calculate_P3(osg::Vec3d& P3, double z3) const;

    double test_orthogonality(const osg::Vec3d& P0, const osg::Vec3d& P1, const osg::Vec3d& P2, const osg::Vec3d& P3) const;
    double test_orthogonality(const osg::Vec3d& cbt, const osg::Vec3d& ctp, int base) const;

    inline double F(int i, int j) const;
    inline std::pair<double, double> quadratic_root_finder_x(double x, double a1, double a2) const;
    inline std::pair<double, double> quadratic_root_finder_y(double y, double a1, double a2) const;
    inline std::pair<double, double> quadratic_root_finder_z(double z, double a1, double a2) const;
};

#endif // COMPONENT_SOLVER_HPP
