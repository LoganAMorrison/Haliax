//
// Created by Logan Morrison on 4/10/20.
//

#ifndef LANRE_OPTIMIZE_QFORM_HPP
#define LANRE_OPTIMIZE_QFORM_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

void qform(int n, Matrix<double> &q) {

    const int m = q.rows();
    Vector<double> wa = Vector<double>::Zero(n);
    const int minmn = std::min(m, n);

    for (int j = 1; j < minmn; j++) {
        for (int i = 0; i <= j - 1; i++) {
            q(i, j) = 0.0;
        }
    }
    //
    //  Initialize remaining columns to those of the identity matrix.
    //
    for (int j = n; j < m; j++) {
        for (int i = 0; i < m; i++) {
            q(i, j) = 0.0;
        }
        q(j, j) = 1.0;
    }

    for (int k = minmn - 1; 0 <= k; k--) {
        for (int i = k; i < m; i++) {
            wa[i] = q(i, k);
            q(i, k) = 0.0;
        }
        q(k, k) = 1.0;

        if (wa[k] != 0.0) {
            for (int j = k; j < m; j++) {
                double sum = 0.0;
                for (int i = k; i < m; i++) {
                    sum = sum + q(i, j) * wa[i];
                }
                double temp = sum / wa[k];
                for (int i = k; i < m; i++) {
                    q(i, j) = q(i, j) - temp * wa[i];
                }
            }
        }
    }
}

}
}

#endif //LANRE_OPTIMIZE_QFORM_HPP
