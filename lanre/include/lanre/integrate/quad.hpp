//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QUAD_HPP
#define LANRE_INTEGRATE_QUAD_HPP

#include "lanre/integrate/qag.hpp"
#include "lanre/integrate/qagi.hpp"
#include "lanre/integrate/qagp.hpp"
#include "lanre/integrate/qags.hpp"
//#include "lanre/integrate/qawc.hpp"
#include <utility>

namespace lanre {
namespace integrate {

enum QuadError {
    Success,
    SubdivisionLimitReached,
    RoundoffDetected,
    BadIntegrand,
    DidNotConverge,
    DivergentOrSlowlyConvergent,
    InvalidInput
};

template<class Type>
class Quad {
public:
    template<class Function, class TypeA, class TypeB>
    static auto integrate(
            Function f,
            TypeA a,
            TypeB b,
            double abstol = 1e-8,
            double reltol = 1e-8,
            size_t limit = 500,
            double *abserr = nullptr,
            int *ier = nullptr
    ) -> decltype(f(a) + f(b)) {
        double aa = a, bb = b;
        double result;
        int neval;
        int ierr;
        double abserrr;


        if (b == a) {
            return 0.0;
        } else if (a > b) {
            aa = b;
            bb = a;
        }

        if (std::abs(aa) >= std::numeric_limits<double>::infinity() ||
                std::abs(bb) >= std::numeric_limits<double>::infinity()) {
            int inf;
            double bound;

            if (bb >= std::numeric_limits<double>::infinity() && aa <= -std::numeric_limits<double>::infinity()) {
                // -inf -> inf
                inf = 2;
                bound = 0.0;
            } else if (bb >= std::numeric_limits<double>::infinity() && aa > -std::numeric_limits<double>::infinity()) {
                inf = 1;
                bound = aa;
            } else {
                inf = -1;
                bound = bb;
            }

            result = qagi(f, bound, inf, abstol, reltol, &abserrr, &neval, &ierr, limit);
        } else {
            result = qags(f, aa, bb, abstol, reltol, &abserrr, &neval, &ierr, limit);
        }

        if (ier != nullptr) {
            *ier = ierr;
        }
        if (abserr != nullptr) {
            *abserr = abserrr;
        }

        return result;
    }

    template<class Function, class TypeA, class TypeB>
    static auto integrate(
            Function f,
            TypeA a,
            TypeB b,
            std::vector<double> brkpts,
            double abstol = 1e-8,
            double reltol = 1e-8,
            size_t limit = 500,
            double *abserr = nullptr,
            int *ier = nullptr
    ) -> decltype(f(a) + f(b)) {
        double aa = a, bb = b;
        double result;
        int neval;
        int ierr;
        double abserrr;


        if (b == a) {
            return 0.0;
        } else if (a > b) {
            aa = b;
            bb = a;
        }

        if (std::abs(aa) >= std::numeric_limits<double>::infinity() ||
                std::abs(bb) >= std::numeric_limits<double>::infinity()) {
            throw std::runtime_error("Quad::integrate: Can't use break-pts on infinite interval.");
        } else {
            result = qagp(f, aa, bb, brkpts.size() + 2, brkpts.data(), abstol, reltol, &abserrr, &neval, &ierr, limit);
        }

        if (ier != nullptr) {
            *ier = ierr;
        }
        if (abserr != nullptr) {
            *abserr = abserrr;
        }

        return result;
    }
};


}
}

#endif //LANRE_INTEGRATE_QUAD_HPP
