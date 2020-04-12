//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_OPTIMIZE_HYBRJ_HPP
#define LANRE_OPTIMIZE_HYBRJ_HPP

#include "lanre/lanre.hpp"
#include "lanre/optimize/qrfac.hpp"
#include "lanre/optimize/qform.hpp"
#include "lanre/optimize/dogleg.hpp"
#include "lanre/optimize/r1updt.hpp"
#include "lanre/optimize/r1mpyq.hpp"
#include <vector>

namespace lanre {
namespace optimize {

template<class Function>
int hybrj_core(
        Function fcn,
        Vector<double> &x,
        Vector<double> &fvec,
        Matrix<double> &fjac,
        int ldfjac,
        double xtol,
        int &maxfev,
        Vector<double> &diag,
        int mode,
        double factor,
        int &nfev,
        int &njev,
        Vector<double> &r,
        Vector<double> &qtf,
        Vector<double> &wa1,
        Vector<double> &wa2,
        Vector<double> &wa3,
        Vector<double> &wa4
) {
    Vector<int> iwa = Vector<int>::Zero(1);
    const double epsmch = std::numeric_limits<double>::epsilon();
    const int n = x.size();
    const int lr = r.size();

    int i, iflag, iter, j, jm1, l, ncfail, ncsuc, nslow1, nslow2, info;
    bool jeval, sing;
    double actred, delta, fnorm, fnorm1, pnorm,
            prered, ratio, sum, temp, xnorm;

    const double p1 = 1.0e-1;
    const double p5 = 5.0e-1;
    const double p001 = 1.0e-3;
    const double p0001 = 1.0e-4;

    info = 0;
    iflag = 0;
    nfev = 0;
    njev = 0;

    if (n <= 0 || ldfjac < n || xtol < 0.0
            || maxfev <= 0 || factor <= 0.0
            || lr < (n * (n + 1)) / 2) {
        goto _300;
    }
    if (mode != 2) {
        goto _20;
    }
    for (j = 0; j < n; j++) {
        if (diag(j) <= 0.0) goto _300;
    }

    _20:
    // evaluate the function at the starting point
    //  and calculate its norm.
    iflag = 1;
    fcn(x, fvec, fjac, iflag);

    nfev = 1;
    if (iflag < 0) goto _300;
    fnorm = fvec.norm();

    // initialize iteration counter and monitors.
    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

    // beginning of the outer loop.
    _30:
    jeval = true;

    // calculate the jacobian matrix.
    iflag = 2;
    fcn(x, fvec, fjac, iflag);
    njev = njev + 1;
    if (iflag < 0) goto _300;

    // compute the qr factorization of the jacobian.
    qrfac(fjac, false, iwa, wa1, wa2);

    // on the first iteration and if mode is 1, scale according
    // to the norms of the columns of the initial jacobian.
    if (iter != 1) goto _70;
    if (mode == 2) goto _50;
    for (j = 0; j < n; j++) {
        diag(j) = wa2(j);
        if (wa2(j) == 0.0) diag(j) = 1.0;
    }


    _50:
    //  on the first iteration, calculate the norm of the scaled x
    // and initialize the step bound delta.
    for (j = 0; j < n; j++) {
        wa3(j) = diag(j) * x(j);
    }
    xnorm = wa3.norm();
    delta = factor * xnorm;
    if (delta == 0.0) delta = factor;

    _70:
    // form (q transpose)*fvec and store in qtf.
    for (i = 0; i < n; i++) {
        qtf(i) = fvec(i);
    }
    for (j = 0; j < n; j++) {
        if (fjac(j, j) == 0.0) continue;
        sum = 0.0;
        for (i = j; i < n; i++) {
            sum = sum + fjac(i, j) * qtf(i);
        }
        temp = -sum / fjac(j, j);
        for (i = j; i < n; i++) {
            qtf(i) = qtf(i) + fjac(i, j) * temp;
        }
    }

    // copy the triangular factor of the qr factorization into r.
    sing = false;
    for (j = 1; j <= n; j++) {
        l = j;
        jm1 = j - 1;
        if (jm1 >= 1) {
            for (i = 1; i <= jm1; i++) {
                r(l - 1) = fjac(i - 1, j - 1);
                l += n - i;
            }
        }
        r(l - 1) = wa1(j - 1);
        if (wa1(j - 1) == 0.0) sing = true;
    }

    // accumulate the orthogonal factor in fjac.
    qform(n, fjac);

    // rescale if necessary.
    if (mode == 2) goto _170;
    for (j = 0; j < n; j++) {
        diag(j) = std::max(diag(j), wa2(j));
    }

    _170:

    // beginning of the inner loop.

    _180:
    // determine the direction p.
    dogleg(r, diag, qtf, delta, wa1, wa2, wa3);

    // store the direction p and x + p. calculate the norm of p.
    for (j = 0; j < n; j++) {
        wa1(j) = -wa1(j);
        wa2(j) = x(j) + wa1(j);
        wa3(j) = diag(j) * wa1(j);
    }
    pnorm = wa3.norm();

    // on the first iteration, adjust the initial step bound.
    if (iter == 1) delta = std::min(delta, pnorm);

    // evaluate the function at x + p and calculate its norm.
    iflag = 1;
    fcn(wa2, wa4, fjac, iflag);
    nfev = nfev + 1;
    if (iflag < 0) goto _300;
    fnorm1 = wa4.norm();

    // compute the scaled actual reduction.
    actred = -1.0;
    if (fnorm1 < fnorm) actred = 1.0 - pow(fnorm1 / fnorm, 2);

    //  compute the scaled predicted reduction.
    l = 0;
    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (j = i; j < n; j++) {
            sum = sum + r(l) * wa1(j);
            l++;
        }
        wa3(i) = qtf(i) + sum;
    }
    temp = wa3.norm();
    prered = 0.0;
    if (temp < fnorm) prered = 1.0 - pow(temp / fnorm, 2);

    // compute the ratio of the actual to the predicted
    // reduction.
    ratio = 0.0;
    if (prered > 0.0) ratio = actred / prered;

    //  update the step bound.
    if (ratio >= p1) goto _230;
    ncsuc = 0;
    ncfail = ncfail + 1;
    delta = p5 * delta;
    goto _240;

    _230:
    ncfail = 0;
    ncsuc = ncsuc + 1;
    if (ratio >= p5 or ncsuc > 1) {
        delta = std::max(delta, pnorm / p5);
    }
    if (std::abs(ratio - 1.0) <= p1) delta = pnorm / p5;

    _240:
    //  test for successful iteration.
    if (ratio < p0001) goto _260;

    // successful iteration. update x, fvec, and their norms.
    for (j = 0; j < n; j++) {
        x(j) = wa2(j);
        wa2(j) = diag(j) * x(j);
        fvec(j) = wa4(j);
    }
    xnorm = wa2.norm();
    fnorm = fnorm1;
    iter = iter + 1;

    _260:
    // determine the progress of the iteration.
    nslow1 = nslow1 + 1;
    if (actred >= p001) nslow1 = 0;
    if (jeval) nslow2++;
    if (actred >= p1) nslow2 = 0;

    // test for convergence.
    if (delta <= xtol * xnorm || fnorm == 0.0) info = 1;
    if (info != 0) goto _300;

    // tests for termination and stringent tolerances.
    if (nfev >= maxfev) info = 2;
    if (p1 * std::max(p1 * delta, pnorm) <= epsmch * xnorm) info = 3;
    if (nslow2 == 5) info = 4;
    if (nslow1 == 10) info = 5;
    if (info != 0) goto _300;

    // criterion for recalculating jacobian.
    if (ncfail == 2) goto _290;

    // calculate the rank one modification to the jacobian
    // and update qtf if necessary.
    for (j = 0; j < n; j++) {
        sum = 0.0;
        for (i = 0; i < n; i++) {
            sum = sum + fjac(i, j) * wa4(i);
        }
        wa2(j) = (sum - wa3(j)) / pnorm;
        wa1(j) = diag(j) * ((diag(j) * wa1(j)) / pnorm);
        if (ratio >= p0001) qtf(j) = sum;
    }

    // compute the qr factorization of the updated jacobian.
    sing = r1updt(r, wa1, wa2, wa3);
    r1mpyq(fjac, wa2, wa3);
    r1mpyq(qtf, wa2, wa3);

    // end of the inner loop.
    jeval = false;
    goto _180;

    _290:
    //  end of the outer loop.
    goto _30;

    _300:
    // termination, either normal or user imposed.
    if (iflag < 0) info = iflag;
    return info;
}

template<class Function>
int hybrj(
        Function fcn,
        Vector<double> &x,
        Vector<double> &fvec,
        double tol
) {
    const int n = x.size();
    double epsfcn;
    double factor;
    int index;
    int info;
    int j;
    int lr;
    int maxfev;
    int ml;
    int mode;
    int mu;
    int nfev;
    int njev;
    int nprint;
    double xtol;
    const int lwa = (n * (3 * n + 13)) / 2;


    info = 0;
    //
    //  Check the input.
    //
    if (n <= 0) {
        return info;
    }
    if (tol <= 0.0) {
        return info;
    }
    if (lwa < (n * (3 * n + 13)) / 2) {
        return info;
    }
    //
    //  Call HYBRD.
    //
    xtol = tol;
    maxfev = 200 * (n + 1);
    ml = n - 1;
    mu = n - 1;
    epsfcn = 0.0;
    mode = 2;
    factor = 100.0;
    nprint = 0;
    nfev = 0;
    njev = 0;
    lr = (n * (n + 1)) / 2;

    Vector<double> diag = Vector<double>::Constant(n, 1.0);
    Matrix<double> fjac = Matrix<double>::Zero(n, n);
    Vector<double> r = Vector<double>::Zero(lr);
    Vector<double> qtf = Vector<double>::Zero(n);
    Vector<double> wa1 = Vector<double>::Zero(n);
    Vector<double> wa2 = Vector<double>::Zero(n);
    Vector<double> wa3 = Vector<double>::Zero(n);
    Vector<double> wa4 = Vector<double>::Zero(n);

    info = hybrj_core(
            fcn,
            x,
            fvec,
            fjac,
            n,
            xtol,
            maxfev,
            diag,
            mode,
            factor,
            nfev,
            njev,
            r,
            qtf,
            wa1,
            wa2,
            wa3,
            wa4
    );

    if (info == 5) {
        info = 4;
    }
    return info;
}

}
}

#endif //LANRE_OPTIMIZE_HYBRJ_HPP
