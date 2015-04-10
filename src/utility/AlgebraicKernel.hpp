#ifndef ALGEBRAIC_KERNEL_HPP
#define ALGEBRAIC_KERNEL_HPP

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpq.h>
#include <vector>

class AlgebraicKernel {
public:
    typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpq>          AK;
    typedef AK::Polynomial_1                                Polynomial_1;
    typedef AK::Algebraic_real_1                            Algebraic_real_1;
    typedef AK::Bound                                       Bound;
    typedef AK::Solve_1                                     Solver_1;
    typedef AK::Approximate_relative_1                      Approximate_relative_1;
    typedef AK::Approximate_absolute_1                      Approximate_absolute_1;
    typedef AK::Multiplicity_type                           Multiplicity_type;

    enum class error_type : unsigned short {
        relative = 0,
        absolute = 1
    };

    AlgebraicKernel();
    void solve(const std::vector<double>& coefficients,
               std::vector<double>& roots,
               const error_type etype = error_type::relative,
               int error_bound = 50);
private:
    AK ak;
    Solver_1 solver;
    Polynomial_1 x;
    std::vector<std::pair<Algebraic_real_1, Multiplicity_type> > m_roots;
    void construct_polynomial(const std::vector<CGAL::Gmpq>& coefficients,
                              Polynomial_1& poly) const;
    inline void initalize_approximator();
    // we do not want clients to copy an algeraic_kernel object.
    AlgebraicKernel(const AlgebraicKernel& rhs);
    AlgebraicKernel& operator=(const AlgebraicKernel& rhs);
};

#endif // ALGEBRAIC_KERNEL_HPP
