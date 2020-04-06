//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_FPKNOT_HPP
#define LANRE_INTERPOLATE_FPKNOT_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * locates an additional knot for a spline of degree
 * k and adjusts the corresponding parameters,i.e.
 * @param t The position of the knots.
 * @param n The number of knots.
 * @param fpint The sum of squares of residual right hand sides
 *              for each knot interval.
 * @param nrdata The number of data points inside each knot interval.
 * @param nrint The number of knot intervals.
 * @param istart Indicates that the smallest data point at which the new knot
 *               may be added is x(istart+1)
 */
void fpknot(
        const double *x,
        const int m,
        double *t,
        int &n,
        double *fpint,
        int *nrdata,
        int &nrint,
        const int nest,
        const int &istart
) {
    double an, am, fpmax;
    int ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg = 0, maxpt = 0,
            next, nrx, number = 0;

    k = (n - nrint - 1) / 2;
    // search for knot interval t(number+k) <= x <= t(number+k+1) where
    // fpint(number) is maximal on the condition that nrdata(number)
    // not equals zero.
    fpmax = 0.0;
    jbegin = istart;
    for (j = 1; j <= nrint; j++) {
        jpoint = nrdata[j - 1];
        if (fpmax >= fpint[j - 1] || jpoint == 0) goto _10;
        fpmax = fpint[j - 1];
        number = j;
        maxpt = jpoint;
        maxbeg = jbegin;
        _10:
        jbegin = jbegin + jpoint + 1;
    }
    // let coincide the new knot t(number+k+1) with a data point x(nrx)
    // inside the old knot interval t(number+k) <= x <= t(number+k+1).
    ihalf = maxpt / 2 + 1;
    nrx = maxbeg + ihalf;
    next = number + 1;
    if (next > nrint) goto _40;
    // adjust the different parameters.
    for (j = next; j <= nrint; j++) {
        jj = next + nrint - j;
        fpint[jj] = fpint[jj - 1];
        nrdata[jj] = nrdata[jj - 1];
        jk = jj + k;
        t[jk] = t[jk - 1];
    }
    _40:
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    am = maxpt;
    an = nrdata[number - 1];
    fpint[number - 1] = fpmax * an / am;
    an = nrdata[next - 1];
    fpint[next - 1] = fpmax * an / am;
    jk = next + k;
    t[jk - 1] = x[nrx - 1];
    n = n + 1;
    nrint = nrint + 1;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPKNOT_HPP
