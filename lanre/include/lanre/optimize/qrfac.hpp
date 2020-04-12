//
// Created by Logan Morrison on 4/10/20.
//

#ifndef LANRE_OPTIMIZE_QRFAC_HPP
#define LANRE_OPTIMIZE_QRFAC_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

void qrfac(
        Matrix<double> &a,
        bool pivot,
        Vector<int> &ipvt,
        Vector<double> &rdiag,
        Vector<double> &acnorm
) {
    const int m = a.rows();
    const int n = a.cols();

    const double p05 = 0.05;
    Vector<double> wa = Vector<double>::Zero(n);
    const double epsmch = std::numeric_limits<double>::epsilon();

    for (int j = 0; j < n; j++) {
        acnorm[j] = a.col(j).norm();
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        if (pivot) {
            ipvt[j] = j;
        }
    }
    //
    //  Reduce A to R with Householder transformations.
    //
    int minmn = std::min(m, n);

    for (int j = 0; j < minmn; j++) {
        if (pivot) {
            //
            //  Bring the column of largest norm into the pivot position.
            //
            int kmax = j;
            for (int k = j; k < n; k++) {
                if (rdiag[kmax] < rdiag[k]) {
                    kmax = k;
                }
            }
            if (kmax != j) {
                for (int i = 0; i < m; i++) {
                    double temp = a(i, j);
                    a(i, j) = a(i, kmax);
                    a(i, kmax) = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                int k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }
        //
        //  Compute the Householder transformation to reduce the
        //  J-th column of A to a multiple of the J-th unit vector.
        //
        double ajnorm = a.col(j).tail(m - j).norm();

        if (ajnorm != 0.0) {
            if (a(j, j) < 0.0) {
                ajnorm = -ajnorm;
            }
            for (int i = j; i < m; i++) {
                a(i, j) /= ajnorm;
            }
            a(j, j) += 1.0;
            //
            //  Apply the transformation to the remaining columns and update the norms.
            //
            for (int k = j + 1; k < n; k++) {
                double sum = 0.0;
                for (int i = j; i < m; i++) {
                    sum = sum + a(i, j) * a(i, k);
                }
                double temp = sum / a(j, j);
                for (int i = j; i < m; i++) {
                    a(i, k) -= temp * a(i, j);
                }
                if (pivot && rdiag[k] != 0.0) {
                    temp = a(j, k) / rdiag[k];
                    rdiag[k] = rdiag[k] * sqrt(std::max(0.0, 1.0 - temp * temp));
                    if (p05 * (rdiag[k] / wa[k]) * (rdiag[k] / wa[k]) <= epsmch) {

                        rdiag[k] = a.col(j + 1).tail(m - 1 + k).norm();
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}

}
}

#endif //LANRE_OPTIMIZE_QRFAC_HPP
