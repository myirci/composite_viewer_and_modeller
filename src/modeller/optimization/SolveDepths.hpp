#ifndef SOLVE_DEPTHS_HPP
#define SOLVED_EPTHS_HPP

#include "ceres/ceres.h"
#include <vector>

class SolveDepths {
public:
    SolveDepths(double far, double near) : f(far), n(near) { g = f + n; }

    struct costfunctor_1 {
        costfunctor_1(unsigned short pr1, unsigned short pr2) :
            p1(pr1), p2(pr2) { }
        template <typename T>
         bool operator()(T const* const* parameters, T* residuals) const {
            return true;
        }
        unsigned short p1;
        unsigned short p2;
    };

    static void solve_for_depths() {

    }

private:
    double f;
    double n;
    double g;
    std::vector<double*> parameters;
};
#endif // SOLVE_DEPTHS_HPP
