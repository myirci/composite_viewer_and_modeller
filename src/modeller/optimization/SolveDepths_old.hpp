#ifndef SOLVE_DEPTHS_HPP
#define SOLVE_DEPTHS_HPP

#include <ceres/ceres.h>

struct F1 {
    F1(double n, double f, double p) : near(n), far(f), x0(p) { }
    template <typename T>
    bool operator()(const T* const xc0, const T* const z0, T* residual) const {
        residual[0] = xc0[0] * z0[0] + xc0[0] * (T(far) + T(near)) - T(far) * T(x0);
        return true;
    }
    double far;
    double near;
    double x0;
};

struct F2 {
    F2(double n, double f, double p) : near(n), far(f), y0(p) { }
    template <typename T>
    bool operator()(const T* const yc0, const T* const z0, T* residual) const {
        residual[0] = yc0[0] * z0[0] + yc0[0] * (T(far) + T(near)) - T(far) * T(y0);
        return true;
    }
    double far;
    double near;
    double y0;
};

struct F3 {
    F3(double n, double f) : near(n), far(f) { }
    template <typename T>
    bool operator()(const T* const zc0, const T* const z0, T* residual) const {
        residual[0] = zc0[0] * z0[0] + zc0[0] * (T(far) + T(near)) + T(far) * T(near);
        return true;
    }
    double far;
    double near;
};

struct F4 {
    F4(double x, double y, double z, double r) : xcc(x), ycc(y), zcc(z), radius(r) { }
    template <typename T>
    bool operator()(const T* const xc0, const T* const yc0, const T* const zc0, T* residual) const {
        residual[0] = (xc0[0] - xcc) * (xc0[0] - xcc) + (yc0[0] - ycc) * (yc0[0] - ycc) + (zc0[0] - zcc) * (zc0[0] - zcc) - T(radius) * T(radius);
        return true;
    }
    double xcc, ycc, zcc, radius;
};

struct F5 {
    F4(double x, double y, double z,
       double n1, double n2, double n3) :
        xcc(x), ycc(y), zcc(z), nx(n1), ny(n2), nz(n3) { }
    template <typename T>
    bool operator()(const T* const xc0, const T* const yc0, const T* const zc0, T* residual) const {
        residual[0] = (xc0[0] - xcc) * T(nx) + (yc0[0] - ycc) * T(ny) + (zc0[0] - zcc) * T(nz);
        return true;
    }
    double nx, ny, nz;
    double xcc, ycc, zcc;
};

void solve_for_depths_old(double far, double near, double x0, double y0, double xcc, double ycc, double zcc, double radius, double n1, double n2, double n3) {

    double xc0(0.0), yc0(0.0), zc0(0.0), z0(0.0);

    ceres::Problem problem;
    problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<F1, 1, 1, 1>(new F1(near, far, x0)),
                NULL, &xc0, &z0);
    problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<F2, 1, 1, 1>(new F2(near, far, y0)),
                NULL, &yc0, &z0);
    problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<F3, 1, 1, 1>(new F3(near, far)),
                NULL, &zc0, &z0);
    problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<F4, 1, 1, 1, 1>(new F4(xcc, ycc, zcc, radius)),
                NULL, &xc0, &yc0, &zc0);
    problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<F5, 1, 1, 1, 1>(new F5(xcc, ycc, zcc, n1, n2, n3)),
                NULL, &xc0, &yc0, &zc0);

    ceres::Solver::Options options;
    ceres::StringToMinimizerType("trust_region", &options.minimizer_type);
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    // options.parameter_tolerance = 0.000000000000001;

    std::cout << "///////////////////////////////////////////////////"<<std::endl;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << "///////////////////////////////////////////////////"<<std::endl;

    n1.normalize();
    n2.normalize();
}



#endif // SOLVE_DEPTHS_HPP
