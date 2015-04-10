#include "AlgebraicKernel.hpp"

AlgebraicKernel::AlgebraicKernel() :
    ak(),
    x(CGAL::shift(Polynomial_1(1),1)),
    solver(ak.solve_1_object()) {
}

void AlgebraicKernel::solve(const std::vector<double>& coefficients,
                            std::vector<double>& roots,
                            const error_type etype,
                            int error_bound) {

    // coefficients must be ordered from lower to higher degree for this
    // implementation. For client interface the order is from higher
    // to lower.
    std::vector<CGAL::Gmpq> coeffs;
    for(auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        coeffs.push_back(CGAL::Gmpq(*it));
    }
    Polynomial_1 poly;
    construct_polynomial(coeffs, poly);
    solver(poly, std::back_inserter(m_roots));
    if(etype == error_type::relative) {
        Approximate_relative_1 approx_r = ak.approximate_relative_1_object();
        for(auto it = m_roots.begin(); it != m_roots.end(); ++it) {
            for(auto it = m_roots.begin(); it != m_roots.end(); ++it) {
                roots.push_back(((approx_r(it->first, error_bound).first +
                                  approx_r(it->first, error_bound).second)/2.0).to_double());
            }
        }
    }
    else if(etype == error_type::absolute) {
        Approximate_absolute_1 approx_a = ak.approximate_absolute_1_object();
        for(auto it = m_roots.begin(); it != m_roots.end(); ++it) {
            roots.push_back(((approx_a(it->first, error_bound).first +
                              approx_a(it->first, error_bound).second)/2.0).to_double());
        }
    }
    m_roots.clear();
}

void AlgebraicKernel::construct_polynomial(const std::vector<CGAL::Gmpq>& coefficients,
                                            Polynomial_1& poly) const {
    int degree = coefficients.size() - 1;
    if(degree < 0) {
        std::cerr << "ERROR: [algebraic_kernel::construct_polynomial]: "
                  << "Polynomial degree below zero!" << std::endl;
        std::cerr << "ERROR: Degree: " << degree << std::endl;
    }
    while(degree > 0) {
        Polynomial_1 p = x;
        for(int i = 0; i < degree-1; i++) { p *= x; }
        poly += coefficients[degree] * p;
        --degree;
    }
    poly += coefficients[0];
}
