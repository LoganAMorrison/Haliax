//
// Created by Logan Morrison on 3/30/20.
//

#ifndef LANRE_INTEGRATE_QAG_HPP
#define LANRE_INTEGRATE_QAG_HPP

#include "lanre/integrate/qsrt.hpp"
#include "lanre/integrate/qext.hpp"
#include "lanre/integrate/qk.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <numeric>

namespace lanre {
namespace integrate {

template<class Integrand, size_t Limit>
double qage(
        Integrand f,
        double a,
        double b,
        double epsabs,
        double epsrel,
        int key,
        double *abserr,
        int *neval,
        int *ier,
        std::array<double, Limit> &alist,
        std::array<double, Limit> &blist,
        std::array<double, Limit> &rlist,
        std::array<double, Limit> &elist,
        std::array<int, Limit> &iord,
        int *last
) {
    static const double epmach = std::numeric_limits<double>::epsilon();
    static const double uflow = std::numeric_limits<double>::min();

    // test on validity of parameters
    *ier = 0;
    *neval = 0;
    *last = 0;
    double result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;

    if (epsabs <= 0.0 && epsrel < std::max(50 * epmach, 0.5e-28)) {
        *ier = 6;
        return NAN;
    }

    // first approximation to the integral
    int keyf = key;
    if (key <= 0) keyf = 1;
    if (key < 7) keyf = 6;
    *neval = 0;
    double defabs;
    double resabs;

    if (keyf == 1) {
        result = qk15(f, a, b, abserr, &defabs, &resabs);
    } else if (keyf == 2) {
        result = qk21(f, a, b, abserr, &defabs, &resabs);
    } else if (keyf == 3) {
        result = qk31(f, a, b, abserr, &defabs, &resabs);
    } else if (keyf == 4) {
        result = qk41(f, a, b, abserr, &defabs, &resabs);
    } else if (keyf == 5) {
        result = qk51(f, a, b, abserr, &defabs, &resabs);
    } else {//keyf == 6
        result = qk61(f, a, b, abserr, &defabs, &resabs);
    }

    *last = 1;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 1;

    // test on accuracy.
    double errbnd = std::max(epsabs, epsrel * std::abs(result));
    if (*abserr <= 50.0 * epmach * defabs && *abserr > errbnd) *ier = 2;
    if (Limit == 1) *ier = 1;
    if (!(*ier != 0 || (*abserr <= errbnd && *abserr != resabs) || *abserr == 0.0)) {
        // initialization
        double errmax = *abserr;
        int maxerr = 1;
        double area = result;
        double errsum = *abserr;
        int nrmax = 1;
        int iroff1 = 0;
        int iroff2 = 0;

        // Main loop
        for ((*last) = 2; (*last) < Limit; (*last)++) {
            // bisect the subinterval with the largest error estimate.
            double a1 = alist[maxerr - 1];
            double b1 = 0.5 * (alist[maxerr - 1] + blist[maxerr - 1]);
            double a2 = b1;
            double b2 = blist[maxerr - 1];

            double area1, error1, defab1;
            if (keyf == 1) {
                area1 = qk15(f, a1, b1, &error1, &resabs, &defab1);
            } else if (keyf == 2) {
                area1 = qk21(f, a1, b1, &error1, &resabs, &defab1);
            } else if (keyf == 3) {
                area1 = qk31(f, a1, b1, &error1, &resabs, &defab1);
            } else if (keyf == 4) {
                area1 = qk41(f, a1, b1, &error1, &resabs, &defab1);
            } else if (keyf == 5) {
                area1 = qk51(f, a1, b1, &error1, &resabs, &defab1);
            } else {//keyf == 6
                area1 = qk61(f, a1, b1, &error1, &resabs, &defab1);
            }

            double area2, error2, defab2;
            if (keyf == 1) {
                area2 = qk15(f, a2, b2, &error2, &resabs, &defab2);
            } else if (keyf == 2) {
                area2 = qk21(f, a2, b2, &error2, &resabs, &defab2);
            } else if (keyf == 3) {
                area2 = qk31(f, a2, b2, &error2, &resabs, &defab2);
            } else if (keyf == 4) {
                area2 = qk41(f, a2, b2, &error2, &resabs, &defab2);
            } else if (keyf == 5) {
                area2 = qk51(f, a2, b2, &error2, &resabs, &defab2);
            } else {//keyf == 6
                area2 = qk61(f, a2, b2, &error2, &resabs, &defab2);
            }

            // improve previous approximations to integral
            // and error and test for accuracy.
            neval = neval + 1;
            double area12 = area1 + area2;
            double erro12 = error1 + error2;
            errsum = errsum + erro12 - errmax;
            area = area + area12 - rlist[maxerr - 1];

            if (!(defab1 == error1 || defab2 == error2)) {
                if (std::abs(rlist[maxerr - 1] - area12) <= 1.0e-05 * std::abs(area12) && erro12 < 0.99 * errmax) {
                    iroff1++;
                }
                if (*last > 10 && erro12 > errmax) iroff2++;
            }
            rlist[maxerr - 1] = area1;
            rlist[*last - 1] = area2;
            errbnd = std::max(epsabs, epsrel * std::abs(area));

            if (errsum > errbnd) {
                // test for roundoff error and eventually set error flag.
                if (iroff1 < 6 || iroff2 < 20) *ier = 2;
                // set error flag in the case that the number of subintervals
                // equals limit.
                if (*last == Limit) *ier = 1;
                // set error flag in the case of bad integrand behaviour
                // at a point of the integration range.
                if (std::max(std::abs(a1), std::abs(b2)) <=
                        (0.1e+01 + 0.1e+03 * epmach) * (std::abs(a2) + 0.1e+04 * uflow)) {
                    *ier = 3;
                }
            }
            //  append the newly-created intervals to the list.
            if (error2 < error1) {
                alist[*last - 1] = a2;
                blist[maxerr - 1] = b1;
                blist[*last - 1] = b2;
                elist[maxerr - 1] = error1;
                elist[*last - 1] = error2;
            } else {
                alist[maxerr - 1] = a2;
                alist[*last - 1] = a1;
                blist[*last - 1] = b1;
                rlist[maxerr - 1] = area2;
                rlist[*last - 1] = area1;
                elist[maxerr - 1] = error2;
                elist[*last - 1] = error1;
            }
            // call subroutine dqpsrt to maintain the descending ordering
            // in the list of error estimates and select the subinterval
            // with the largest error estimate (to be bisected next).
            qsrt(*last, &maxerr, &errmax, elist, iord, &nrmax);
            if (*ier != 0 || errsum <= errbnd) break;
        }
        //  compute final result.
        result = std::accumulate(rlist.begin(), rlist.begin() + *last, 0.0);
        *abserr = errsum;
    }

    // 60
    if (keyf != 1) *neval = (10 * keyf + 1) * (2 * (*neval) + 1);
    if (keyf == 1) *neval = 30 * (*neval) + 15;
    return result;
}


template<class Integrand, size_t Limit = 500>
double qag(
        Integrand f,
        double a,
        double b,
        double epsabs,
        double epsrel,
        int key,
        double *abserr,
        int *neval,
        int *ier,
        int *last,
        std::array<double, Limit> &alist,
        std::array<double, Limit> &blist,
        std::array<double, Limit> &rlist,
        std::array<double, Limit> &elist,
        std::array<int, Limit> &iord
) {
    // check validity of lenw.
    *ier = 6;
    *neval = 0;
    *last = 0;
    double result = 0.0;
    *abserr = 0.0;
    if (Limit >= 1) {
        // prepare call for dqage.
        result = qage(f, a, b, epsabs, epsrel, key, abserr, neval,
                      ier, alist, blist, rlist, elist, iord, last);
    }
    return result;
}


template<class Integrand, size_t Limit = 500>
double qag(
        Integrand f,
        double a,
        double b,
        double epsabs,
        double epsrel,
        int key,
        double *abserr
) {
    // check validity of lenw.
    int ier = 6;
    int neval = 0;
    int last = 0;
    double result = 0.0;
    *abserr = 0.0;
    if (Limit >= 1) {
        // prepare call for dqage.
        std::array<double, Limit> alist = {0};
        std::array<double, Limit> blist = {0};
        std::array<double, Limit> rlist = {0};
        std::array<double, Limit> elist = {0};
        std::array<int, Limit> iord = {0};

        result = qage(f, a, b, epsabs, epsrel, key, abserr, &neval,
                      &ier, alist, blist, rlist, elist, iord, &last);
    }
    return result;
}

}

}

#endif //LANRE_INTEGRATE_QAG_HPP
