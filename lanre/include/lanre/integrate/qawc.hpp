//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QAWC_HPP
#define LANRE_INTEGRATE_QAWC_HPP

#include "lanre/integrate/base.hpp"
#include "lanre/integrate/qmomo.hpp"
#include "lanre/integrate/qext.hpp"


namespace lanre {
namespace integrate {

template<class Function>
double dqawc(Function f, double a, double b, double c, double epsabs,
             double epsrel, double *abserr, int *neval, int *ier) {
    double aa, area, area1, area2, area12, a1, a2, bb, b1, b2;
    double errbnd, errmax, error1, error2, erro12, errsum, result;
    double alist[kLIMIT], blist[kLIMIT], rlist[kLIMIT], elist[kLIMIT];

    int iord[kLIMIT], iroff1, iroff2, k, krule, last, maxerr, nrmax, nev;
    int limit;

    limit = kLIMIT - 1;
    *ier = 6;
    *neval = 0;
    last = 0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((c == a) || (c == b) || ((epsabs < 0.0) && (epsrel < 0.0))) goto _999;

    /*  First approximation to the integral.    */
    aa = a;
    bb = b;
    if (a <= b) goto _10;
    aa = b;
    bb = a;
    _10:
    *ier = 0;
    krule = 1;
    result = qc25c(f, aa, bb, c, abserr, &krule, neval);
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    alist[0] = a;
    blist[0] = b;

    /*  Test on accuracy.   */
    errbnd = fmax(epsabs, epsrel * fabs(result));
    if (limit == 0) *ier = 1;
    if ((*abserr < fmin(1.0e-2 * fabs(result), errbnd)) || (*ier == 1))
        goto _70;

    /*  Initialization. */
    alist[0] = aa;
    blist[0] = bb;
    rlist[0] = result;
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

    /*  Main loop.  */
    for (last = 1; last < limit; last++) {
        /* Bisect the subinterval with nrmax-th largest error estimate.    */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        b2 = blist[maxerr];
        if ((c <= b1) && (c > a1)) b1 = 0.5 * (c + b2);
        if ((c > b1) && (c < b2)) b1 = 0.5 * (a1 + c);
        a2 = b1;
        krule = 2;
        area1 = qc25c(f, a1, b1, c, &error1, &krule, &nev);
        *neval = *neval + nev;
        area2 = qc25c(f, a2, b2, c, &error2, &krule, &nev);
        *neval = *neval + nev;

        /*  Improve previous approximations to integral and error and
         *  test for accuracy.
         */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12 - errmax);
        area += (area12 - rlist[maxerr]);
        if (fabs(rlist[maxerr] - area12) < (1.0e-5 * fabs(area12)) &&
                (erro12 >= 0.99 * errmax) && (krule == 0))
            iroff1++;
        if ((last > 10) && (erro12 > errmax) && (krule == 0))
            iroff2++;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = fmax(epsabs, epsrel * fabs(area));
        if (errsum <= errbnd) goto _15;

        /*  Test for roundoff error and eventually set error flag. */
        if ((iroff1 >= 6) && (iroff2 > 20)) *ier = 2;

        /*  Set error flag in the case that number of interval bisections
         *  exceeds limit.
         */
        if (last == limit) *ier = 1;

        /*  Set error flag in the case of bad integrand behavior at a point
         *  of the integration range.
         */
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 1.0e3 * kEPMACH) *
                (fabs(a2) + 1.0e3 * kUFLOW))
            *ier = 3;
        /* Append newly created intervals to list. */
        _15:
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        } else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

        /*  Call subroutine qsrt to maintain descending ordering in the list. */
        qsrt(limit, last, &maxerr, &errmax, elist, iord, &nrmax);

        /*  Jump out of loop.   */
        if ((*ier != 0) || (errsum <= errbnd)) goto _50;
    }

    /*  Compute final result.   */
    _50:
    result = 0.0;
    for (k = 0; k <= last; k++) {
        result += rlist[k];
    }
    *abserr = errsum;
    _70:
    if (aa == b) result = -result;
    _999:
    return result;
}


template<class Function>
double dqawf(Function f, double a, double omega, int sincos,
             double epsabs, int limlst, int maxp1,
             double *abserr, int *neval, int *ier, double rslst[],
             double erlst[], int ierlst[], double **chebmo) {
    double abseps, correc, cycle, c1, c2, dl, drl;
    double ep, eps, epsa, errsum, fact, p1, psum[52], reseps;
    double res3la[3], result;

    int ktmin, l, ll, momcom, nev, nres, numrl2, lst;

    /* Test on validity of parameters. */
    result = 0.0;
    *abserr = 0.0;
    *neval = 0;
    *ier = 0;
    ll = 0;
    if (((sincos != kCOSINE) && (sincos != kSINE)) || (epsabs <= 0.0) ||
            (limlst < 3))
        *ier = 6;
    if (*ier == 6) goto _999;
    if (omega != 0.0) goto _10;

    /* Integration by DQAGI if omega is zero. */
    if (sincos == kCOSINE)
        result = qagi(f, 0.0, 1, epsabs, 0.0, abserr, neval, ier);
    rslst[0] = result;
    erlst[0] = *abserr;
    ierlst[0] = *ier;
    goto _999;

    /* Initialization. */
    _10:
    res3la[0] = 0.0;    /* res3la must be initialized to 0.0 */
    res3la[1] = 0.0;
    res3la[2] = 0.0;
    l = fabs(omega);
    dl = 2 * l + 1;
    cycle = dl * M_PI / fabs(omega);
    *ier = 0;
    ktmin = 0;
    *neval = 0;
    numrl2 = -1;    /* used as array index. first use is after increment. */
    nres = 0;
    c1 = a;
    c2 = cycle + a;
    p1 = 1.0 - 0.9;
    eps = epsabs;
    if (epsabs > (kUFLOW / p1))
        eps = epsabs * p1;
    ep = eps;
    fact = 1.0;
    correc = 0.0;
    *abserr = 0.0;
    errsum = 0.0;

    /* Main Loop */
    for (lst = 0; lst < limlst; lst++) {

        /* Integrate over current subinterval. */
        /*    dla = lst;  This line is in the original code, but dla is unused. */
        epsa = eps * fact;
        rslst[lst] = dqfour(f, c1, c2, omega, sincos, epsa, 0.0, lst + 1, maxp1,  // lst+1
                            &erlst[lst], &nev, &ierlst[lst], &momcom, chebmo);
        *neval += nev;
        fact *= 0.9;
        errsum += erlst[lst];
        drl = 50.0 * fabs(rslst[lst]);

        /* Test on accuracy with partial sum. */
        if (((errsum + drl) <= epsabs) && (lst >= 5))
            goto _80;
        correc = fmax(correc, erlst[lst]);
        if (ierlst[lst] != 0)
            eps = fmax(ep, correc * p1);
        if (ierlst[lst] != 0)
            *ier = 7;
        if ((*ier == 7) && ((errsum + drl) <= (correc * 10.0))
                && (lst > 4))
            goto _80;
        numrl2++;
        if (lst > 0)
            goto _20;
        psum[0] = rslst[0];
        goto _40;
        _20:
        psum[numrl2] = psum[ll] + rslst[lst];

        if (lst == 1)
            goto _40;

        /* Test on maximum number of subintervals. */
        if (lst == limlst - 1)
            *ier = 8;

        /* Perform new extrapolation. */
        reseps = qext(&numrl2, psum, &abseps, res3la, &nres);

        /* Test whether extrapolated result is influenced by roundoff. */
        ktmin++;
        if ((ktmin >= 15) && (*abserr <= 0.001 * (errsum + drl)))
            *ier = 9;
        if ((abseps > *abserr) && (lst != 2))
            goto _30;
        *abserr = abseps;
        result = reseps;
        ktmin = 0;

        /* If ier is not 0, check whether direct result (partial sum) or
         * extrapolated result yields the best integral approximation.
         */
        if (((*abserr + 10.0 * correc) <= epsabs) || (*abserr <= epsabs) &&
                (10.0 * correc >= epsabs))
            goto _60;
        _30:
        if ((*ier != 0) && (*ier != 7))
            goto _60;
        _40:
        ll = numrl2;
        c1 = c2;
        c2 += cycle;
        _50:;
    }

    /* Set final result and error estimate. */
    _60:
    (*abserr) += (10.0 * correc);
    if (*ier == 0)
        goto _999;
    if ((result != 0.0) && (psum[numrl2] != 0.0))
        goto _70;
    if (*abserr > errsum)
        goto _80;
    if (psum[numrl2] == 0.0)
        goto _999;
    _70:
    if ((*abserr / fabs(result) > (errsum + drl) / fabs(psum[numrl2])))
        goto _80;
    if ((*ier >= 1) && (*ier != 7))
        (*abserr) += drl;
    goto _999;
    _80:
    result = psum[numrl2];
    *abserr = errsum + drl;
    _999:
    return result;

}


template<class Function>
double qaws(Function f, double a, double b, double alfa, double beta,
            int wgtfunc, double epsabs, double epsrel, double *abserr,
            int *neval, int *ier) {
    double alist[kLIMIT], blist[kLIMIT], rlist[kLIMIT], elist[kLIMIT];
    double ri[25], rj[25], rh[25], rg[25];
    double area, area1, area12, area2, a1, a2, b1, b2, centre;
    double errbnd, errmax, error1, erro12, error2, errsum;
    double resas1, resas2, result;

    int iord[kLIMIT], iroff1, iroff2, k, last, limit, maxerr, nev, nrmax;

    limit = kLIMIT - 1;
    /*  Test on validity of parameters. */
    *ier = 6;
    *neval = 0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((b <= a) || ((epsabs < 0.0) && (epsrel < 0.0)) ||
            (alfa <= -1.0) || (beta <= -1.0) || (wgtfunc < 1) ||
            (wgtfunc > 4) || (limit < 1))
        goto _999;
    *ier = 0;

    /*  Compute the modified Chebyshev moments. */
    qmomo(alfa, beta, ri, rj, rg, rh, wgtfunc);

    /*  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b). */
    centre = 0.5 * (a + b);
    area1 = qc25s(f, a, b, a, centre, alfa, beta, ri, rj, rg, rh, &error1,
                  &resas1, wgtfunc, &nev);
    *neval = *neval + nev;
    area2 = qc25s(f, a, b, centre, b, alfa, beta, ri, rj, rg, rh, &error2,
                  &resas2, wgtfunc, &nev);
    *neval = *neval + nev;
    result = area1 + area2;
    *abserr = error1 + error2;

    /* Test on accuracy. */
    errbnd = fmax(epsabs, epsrel * fabs(result));

    /*  Initialization. */
    if (error1 >= error2) {
        alist[0] = a;
        alist[1] = centre;
        blist[0] = centre;
        blist[1] = b;
        rlist[0] = area1;
        rlist[1] = area2;
        elist[0] = error1;
        elist[1] = error2;
    } else {
        alist[0] = centre;
        alist[1] = a;
        blist[0] = b;
        blist[1] = centre;
        rlist[0] = area2;
        rlist[1] = area1;
        elist[0] = error2;
        elist[1] = error1;
    }
    iord[0] = 0;
    iord[1] = 1;
    if (limit == 1) *ier = 1;
    if ((*abserr <= errbnd) || (*ier == 1)) goto _999;
    errmax = elist[0];
    maxerr = 0;
    nrmax = 0;
    area = result;
    errsum = maxerr;
    iroff1 = 0;
    iroff2 = 0;

    /*  Main loop. */
    for (last = 2; last < limit; last++) {
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];

        area1 = qc25s(f, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, &error1,
                      &resas1, wgtfunc, &nev);
        *neval = *neval + nev;
        area2 = qc25s(f, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, &error2,
                      &resas2, wgtfunc, &nev);
        *neval = *neval + nev;

        /*  Improve previous approximation and error test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12 - errmax);
        area += (area12 - rlist[maxerr]);
        if ((a == a1) || (b == b2)) goto _30;
        if ((resas1 == error1) || (resas2 == error2)) goto _30;

        /*  Test for roundoff error. */
        if ((fabs(rlist[maxerr] - area12) < (1.0e-5 * fabs(area12))) &&
                (erro12 >= (0.99 * errmax)))
            iroff1++;
        if ((last > 9) && (erro12 > errmax)) iroff2++;
        _30:
        rlist[maxerr] = area1;
        rlist[last] = area2;

        /*  Test on accuracy. */
        errbnd = fmax(epsabs, epsrel * fabs(area));
        if (errsum <= errbnd) goto _35;

        /*  Set error flag in the case that number of intervals exceeds limit. */
        if (last == limit) *ier = 1;

        /*  Set error flag in the case of roundoff error. */
        if ((iroff1 > 5) || (iroff2 > 19)) *ier = 2;

        /*  Set error flag in case of bad integrand behavior at interior points. */
        if (fmax(fabs(a1), fabs(b2)) <= ((1.0 + 1.0e3 * kEPMACH) *
                (fabs(a2) + 1.0e3 * kUFLOW)))
            *ier = 3;
        /*  Append the newly created intervals to the list. */
        _35:
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        } else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

        /*  Call subroutine qsort to maintain the descending ordering. */
        qsrt(limit, last, &maxerr, &errmax, elist, iord, &nrmax);

        /*  Jump out of loop. */
        if ((*ier != 0) || (errsum <= errbnd)) break;
    }
    result = 0.0;
    for (k = 0; k <= last; k++) {
        result += rlist[k];
    }
    *abserr = errsum;
    _999:
    return result;
}

}
}

#endif //LANRE_INTEGRATE_QAWC_HPP
