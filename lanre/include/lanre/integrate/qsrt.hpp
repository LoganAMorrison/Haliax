//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QSRT_HPP
#define LANRE_INTEGRATE_QSRT_HPP

#include <vector>

namespace lanre {
namespace integrate {

/**
 *  maintains the descending ordering in the
 *  list of the local error estimated resulting from the
 *  interval subdivision process. at each call two error
 *  estimates are inserted using the sequential search
 *  method, top-down for the largest error estimate and
 *  bottom-up for the smallest error estimate.
 *
 * @paramt limit maximum number of error estimates the list
 *               can contain
 * @param last  number of error estimates currently in the list
 * @param maxerr maxerr points to the nrmax-th largest error
 *               estimate currently in the list
 * @param ermax nrmax-th largest error estimate
 *              ermax = elist(maxerr)
 * @param elist vector of dimension last containing
 *              the error estimates
 * @param iord vector of dimension last, the first k elements
 *             of which contain pointers to the error
 *             estimates, such that
 *             elist(iord(1)),...,  elist(iord(k))
 *             form a decreasing sequence, with
 *             k = last if last<=(limit/2+2), and
 *             k = limit+1-last otherwise
 * @param nrmax maxerr = iord(nrmax)
 */
template<size_t Limit>
void qsrt(
        int last,
        int *maxerr,
        double *ermax,
        const std::array<double, Limit> &elist,
        std::array<int, Limit> &iord,
        int *nrmax
) {
    int i, ibeg, ido, isucc, j, jbnd, jupbn, k;
    double errmax, errmin;


    if (last > 2) goto _10;
    iord[0] = 1;
    iord[1] = 2;
    goto _90;

    // this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    _10:
    errmax = elist[*maxerr - 1];
    if (*nrmax == 1) goto _30;
    ido = *nrmax - 1;
    for (i = 1; i <= ido; i++) {
        isucc = iord[*nrmax - 2];
        if (errmax <= elist[isucc - 1]) goto _30;
        iord[*nrmax - 1] = isucc;
        nrmax = nrmax - 1;
    }

    // compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    _30:
    jupbn = last;
    if (last > (Limit / 2 + 2)) jupbn = Limit + 3 - last;
    errmin = elist[last - 1];

    // insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    jbnd = jupbn - 1;
    ibeg = *nrmax + 1;
    if (ibeg > jbnd) goto _50;
    for (i = ibeg; i <= jbnd; i++) {
        isucc = iord[i - 1];
        if (errmax < elist[isucc - 1]) goto _60;
        iord[i - 2] = isucc;
    }

    _50:
    iord[jbnd - 1] = *maxerr;
    iord[jupbn - 1] = last;
    goto _90;

    // insert errmin by traversing the list bottom-up.
    _60:
    iord[i - 2] = *maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; j++) {
        isucc = iord[k - 1];
        if (errmin < elist[isucc - 1]) goto _80;
        iord[k] = isucc;
        k = k - 1;
    }
    iord[i - 1] = last;
    goto _90;

    _80:
    iord[k] = last;

    // set maxerr and ermax.
    _90:
    *maxerr = iord[*nrmax - 1];
    *ermax = elist[*maxerr - 1];
}

}
}

#endif //LANRE_INTEGRATE_QSRT_HPP
