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
void qsrt(int limit, int last, int *maxerr, double *ermax, double elist[],
          int iord[], int *nrmax) {
    double errmax, errmin;
    int i, ibeg, ido, isucc, j, jbnd, jupbn, k;

    if (last > 1) goto _10;
    iord[0] = 0;
    iord[1] = 1;
    goto _90;
    _10:
    errmax = elist[*maxerr];
    if (*nrmax == 0) goto _30;
    ido = (*nrmax) - 1;
    for (i = 0; i <= ido; i++) {
        isucc = iord[*nrmax - 1];
        if (errmax <= elist[isucc]) goto _30;
        iord[*nrmax] = isucc;
        (*nrmax)--;
    }
    _30:
    jupbn = last;
    if (last > (limit / 2 + 2))
        jupbn = limit + 3 - last;
    errmin = elist[last];
    jbnd = jupbn - 1;
    ibeg = *nrmax + 1;
    if (ibeg > jbnd) goto _50;
    for (i = ibeg; i <= jbnd; i++) {
        isucc = iord[i];
        if (errmax >= elist[isucc]) goto _60;
        iord[i - 1] = isucc;
    }
    _50:
    iord[jbnd] = *maxerr;
    iord[jupbn] = last;
    goto _90;
    _60:
    iord[i - 1] = *maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; j++) {
        isucc = iord[k];
        if (errmin < elist[isucc]) goto _80;
        iord[k + 1] = isucc;
        k--;
    }
    iord[i] = last;
    goto _90;
    _80:
    iord[k + 1] = last;
    _90:
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
}

}
}

#endif //LANRE_INTEGRATE_QSRT_HPP
