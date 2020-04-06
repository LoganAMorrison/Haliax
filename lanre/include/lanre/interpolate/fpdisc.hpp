// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPDISC_HPP
#define LANRE_INTERPOLATE_FPDISC_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

void fpdisc(const double *t, const int n, const int k2, double **b, const int nest) {
    double an, fac, prod;
    int i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint;
    double h[12];

    k1 = k2 - 1;
    k = k1 - 1;
    nk1 = n - k1;
    nrint = nk1 - k;
    an = nrint;
    fac = an / (t[nk1] - t[k1 - 1]);
    for (l = k2; l <= nk1; l++) {
        lmk = l - k1;
        for (j = 1; j <= k1; j++) {
            ik = j + k1;
            lj = l + j;
            lk = lj - k2;
            h[j - 1] = t[l - 1] - t[lk - 1];
            h[ik - 1] = t[l - 1] - t[lj - 1];
        }
        lp = lmk;
        for (j = 1; j <= k2; j++) {
            jk = j;
            prod = h[j - 1];
            for (i = 1; i <= k; i++) {
                jk = jk + 1;
                prod = prod * h[jk - 1] * fac;
            }
            lk = lp + k1;
            b[lmk - 1][j - 1] = (t[lk - 1] - t[lp - 1]) / prod;
            lp = lp + 1;
        }
    }
}


} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPDISC_HPP