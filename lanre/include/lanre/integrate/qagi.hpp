//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QAGI_HPP
#define LANRE_INTEGRATE_QAGI_HPP

#include "lanre/integrate/base.hpp"
#include "lanre/integrate/qsrt.hpp"
#include "lanre/integrate/qext.hpp"
#include "lanre/integrate/qk.hpp"
#include <cmath>
#include <array>
#include <numeric>

namespace lanre {
namespace integrate {

template<typename Integrand, size_t Limit>
double qagie(
        Integrand f,
        double bound,
        int inf,
        double epsabs,
        double epsrel,
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
    std::array<double, 3> res3la = {0};
    std::array<double, 52> rlist2 = {0};
    static const double epmach = std::numeric_limits<double>::epsilon();
    static const double uflow = std::numeric_limits<double>::min();
    static const double oflow = std::numeric_limits<double>::max();

    double abseps, area, area1, area12, area2, a1,
            a2, boun, b1, b2, correc = 0.0, defabs, defab1, defab2,
            dres, erlarg = 0.0, erlast,
            errbnd, errmax, error1, error2, erro12, errsum, ertest = 0.0, resabs,
            reseps, result, small = 0.0;
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn,
            ktmin, maxerr, nres, nrmax, numrl2;
    bool extrap, noext;

    // test on validity of parameters
    *ier = 0;
    *neval = 0;
    *last = 0;
    *abserr = 0.0;
    alist[0] = 0.0;
    blist[0] = 0.1e1;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    if (epsabs <= 0.0 && epsrel < std::max(0.5e2 * epmach, 0.5 - 28)) {
        *ier = 6;
        return NAN;
    }

    // first approximation to the integral

    // determine the interval to be mapped onto (0,1).
    // if inf = 2 the integral is computed as i = i1+i2, where
    // i1 = integral of f over (-infinity,0),
    // i2 = integral of f over (0,+infinity).
    boun = bound;
    if (inf == 2) boun = 0.0;
    result = qk15i(f, boun, inf, 0.0, 1.0, abserr, &defabs, &resabs);

    //  test on accuracy
    *last = 1;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 1;
    dres = std::abs(result);
    errbnd = std::max(epsabs, epsrel * dres);
    if (*abserr <= 1.0e2 * epmach * defabs && *abserr > errbnd) *ier = 2;
    if (Limit == 1) *ier = 1;
    if (*ier != 0 || (*abserr <= errbnd && *abserr != resabs) || *abserr == 0.0) {
        goto _130;
    }

    // initialization
    rlist2[0] = result;
    errmax = *abserr;
    maxerr = 1;
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    ktmin = 0;
    numrl2 = 2;
    extrap = false;
    noext = false;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres < (1.0 - 50.0 * epmach) * defabs) ksgn = 1;

    // Main loop
    for ((*last) = 2; (*last) <= Limit; (*last)++) {
        // bisect the subinterval with nrmax-th largest error estimate.
        a1 = alist[maxerr - 1];
        b1 = 0.5 * (alist[maxerr - 1] + blist[maxerr - 1]);
        a2 = b1;
        b2 = blist[maxerr - 1];
        erlast = errmax;
        area1 = qk15i(f, boun, inf, a1, b1, &error1, &resabs, &defab1);
        area2 = qk15i(f, boun, inf, a2, b2, &error2, &resabs, &defab2);

        // improve previous approximations to integral
        // and error and test for accuracy.
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr - 1];

        if (defab1 == error1 || defab2 == error2) goto _15;
        if (std::abs(rlist[maxerr - 1] - area12) > 0.1e-4 * std::abs(area12) || erro12 < 0.99 * errmax) goto _10;
        if (extrap) iroff2 = iroff2 + 1;
        if (!extrap) iroff1 = iroff1 + 1;

        _10:
        if (*last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;

        _15:
        rlist[maxerr - 1] = area1;
        rlist[*last - 1] = area2;
        errbnd = std::max(epsabs, epsrel * std::abs(area));

        // test for roundoff error and eventually set error flag.
        if (iroff1 + iroff2 < 10 || iroff3 < 20) *ier = 2;
        if (iroff2 < 5) ierro = 3;

        // set error flag in the case that the number of
        // subintervals equals limit.
        if (*last == Limit) *ier = 1;

        // set error flag in the case of bad integrand behaviour
        // at some points of the integration range.
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 100.0 * epmach) *
                (std::abs(a2) + 1000.0 * uflow)) {
            *ier = 4;
        }

        //  append the newly-created intervals to the list.
        if (error2 > error1) {
            alist[maxerr - 1] = a2;
            alist[*last - 1] = a1;
            blist[*last - 1] = b1;
            rlist[maxerr - 1] = area2;
            rlist[*last - 1] = area1;
            elist[maxerr - 1] = error2;
            elist[*last - 1] = error1;
        } else {
            alist[*last - 1] = a2;
            blist[maxerr - 1] = b1;
            blist[*last - 1] = b2;
            elist[maxerr - 1] = error1;
            elist[*last - 1] = error2;
        }
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        qsrt(*last, &maxerr, &errmax, elist, iord, &nrmax);
        if (errsum <= errbnd) goto _115;
        if (*ier != 0) goto _100;
        if (*last == 2) goto _80;
        if (noext) continue;
        erlarg = erlarg - erlast;
        if (std::abs(b1 - a1) > small) erlarg = erlarg + erro12;
        if (extrap) goto _40;

        // test whether the interval to be bisected next is the
        // smallest interval.
        if (std::abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) continue;
        extrap = true;
        nrmax = 2;

        _40:
        if (ierro == 3 || erlarg <= ertest) goto _60;
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        id = nrmax;
        jupbnd = *last;
        if (*last > (2 + Limit / 2)) jupbnd = Limit + 3 - *last;
        for (k = id; k <= jupbnd; k++) {
            maxerr = iord[nrmax - 1];
            errmax = elist[maxerr - 1];
            if (std::abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) goto _90;
            nrmax = nrmax + 1;
        }

        // perform extrapolation.
        _60:
        numrl2 = numrl2 + 1;
        rlist2[numrl2 - 1] = area;
        reseps = qext(&numrl2, rlist2, &abseps, res3la, &nres);
        ktmin = ktmin + 1;
        if (ktmin > 5 && *abserr < 0.1e-2 * errsum) *ier = 5;
        if (abseps < *abserr) goto _70;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::abs(reseps));
        if (*abserr <= ertest) goto _100;

        _70:
        // prepare bisection of the smallest interval.
        if (numrl2 == 1) noext = true;
        if (*ier == 5) goto _100;
        maxerr = iord[0];
        errmax = elist[maxerr - 1];
        nrmax = 1;
        extrap = false;
        small *= 0.5;
        erlarg = errsum;
        continue;

        _80:
        small = 0.375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;

        _90:;
    }

    // set final result and error estimate.
    _100:
    if (*abserr == oflow) goto _115;
    if ((*ier + ierro) == 0) goto _110;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if (result != 0.0 && area != 0.0)goto _105;
    if (*abserr > errsum) goto _115;
    if (area == 0.0) goto _130;
    goto _110;

    _105:
    if (*abserr / std::abs(result) > errsum / std::abs(area)) goto _115;

    //  test on divergence
    _110:
    if (ksgn == (-1) && std::max(std::abs(result), std::abs(area)) <= defabs * 0.1e-1) {
        goto _130;
    }
    if (0.1e-1 > (result / area) || (result / area) > 0.1e3 || errsum > std::abs(area)) {
        *ier = 6;
    }

    // compute global integral sum.
    _115:
    result = 0.0;
    for (k = 0; k < (*last); k++) {
        result += rlist[k];
    }
    //result = std::accumulate(rlist.begin(), rlist.begin() + (*last), 0.0);
    *abserr = errsum;

    _130:
    *neval = 30 * (*last) - 15;
    if (inf == 2) (*neval) *= 2;
    if (*ier > 2) (*ier)--;
    return result;
}


template<typename Integrand, size_t Limit = 500>
double qagi(
        Integrand f,
        double bound,
        int inf,
        double epsabs,
        double epsrel,
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
    *ier = 6;
    *neval = 0;
    *last = 0;
    *abserr = 0.0;
    if (Limit >= 1) {
        return qagie(f, bound, inf, epsabs, epsrel, abserr,
                     neval, ier, alist, blist, rlist, elist, iord, last);
    } else {
        return NAN;
    }
}

template<typename Integrand, size_t Limit = 500>
double qagi(
        Integrand f,
        double bound,
        int inf,
        double epsabs,
        double epsrel,
        double *abserr
) {
    int ier = 6;
    int neval = 0;
    int last = 0;
    *abserr = 0.0;
    if (Limit >= 1) {
        std::array<double, Limit> alist = {0};
        std::array<double, Limit> blist = {0};
        std::array<double, Limit> rlist = {0};
        std::array<double, Limit> elist = {0};
        std::array<int, Limit> iord = {0};
        double result = qagie(f, bound, inf, epsabs, epsrel, abserr,
                              &neval, &ier, alist, blist, rlist, elist, iord, &last);
        return result;
    } else {
        return NAN;
    }
}

}
}

#endif //LANRE_INTEGRATE_QAGI_HPP
