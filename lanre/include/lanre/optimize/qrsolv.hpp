//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_OPTIMIZE_QRSOLV_HPP
#define LANRE_OPTIMIZE_QRSOLV_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

/**
 * Solves a rectangular linear system A*x=b in the least squares sense.
 *
 * Given an M by N matrix A, an N by N diagonal matrix D,
 * and an M-vector B, the problem is to determine an X which
 * solves the system:
 *      A*x = b
 *      D.x = 0
 * in the least squares sense.
 *
 * This function completes the solution of the problem
 * if it is provided with the necessary information from the
 * QR factorization, with column pivoting, of A.  That is, if
 * A*P = Q*R, where P is a permutation matrix, Q has orthogonal
 * columns, and R is an upper triangular matrix with diagonal
 * elements of nonincreasing magnitude, then QRSOLV expects
 * the full upper triangle of R, the permutation matrix p,
 * and the first N components of Q'*B.
 *
 * The system is then equivalent to:
 *
 *      R*Z = Q'*B
 *      P'*D*P*Z = 0
 *
 *  where X = P*Z.  If this system does not have full rank,
 *  then a least squares solution is obtained.  On output QRSOLV
 *  also provides an upper triangular matrix S such that:
 *
 *          P'*(A'*A + D*D)*P = S'*S.
 *
 *   S is computed within QRSOLV and may be of separate interest.
 *
 * @param r On input the full upper triangle must contain the full upper triangle
 *          of the matrix R.
 *          On output the full upper triangle is unaltered, and
 *          the strict lower triangle contains the strict upper triangle
 *          (transposed) of the upper triangular matrix S.
 * @param ipvt defines the permutation matrix P such
 *             that A*P = Q*R.  Column j of p is column ipvt[j-1] of the identity.
 * @param diag the diagonal elements of the matrix D.
 * @param qtb the first N elements of the vector Q'*B.
 * @param x the least squares solution.
 * @param sdiag the diagonal elements of the upper triangular matrix S.
 */
void qrsolv(
        Matrix<double> &r,
        Vector<int> &ipvt,
        Vector<double> &diag,
        Vector<double> &qtb,
        Vector<double> &x,
        Vector<double> &sdiag
) {
    double c;
    double cotan;
    int i;
    int j;
    int k;
    int l;
    int nsing;
    double qtbpj;
    double s;
    double sum2;
    double t;
    double temp;

    const int n = ipvt.size();
    const int ldr = r.rows();

    std::vector<double> wa(n);

    //
    //  Copy R and Q'*B to preserve input and initialize S.
    //
    //  In particular, save the diagonal elements of R in X.
    //
    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            r(i, j) = r(j, i);
        }
    }
    for (j = 0; j < n; j++) {
        x[j] = r(j, j);
    }

    for (j = 0; j < n; j++) {
        wa[j] = qtb[j];
    }
    //
    //  Eliminate the diagonal matrix D using a Givens rotation.
    //
    for (j = 0; j < n; j++) {
        //
        //  Prepare the row of D to be eliminated, locating the
        //  diagonal element using P from the QR factorization.
        //
        l = ipvt[j];

        if (diag[l] != 0.0) {
            sdiag[j] = diag[l];
            for (i = j + 1; i < n; i++) {
                sdiag[i] = 0.0;
            }
            //
            //  The transformations to eliminate the row of D
            //  modify only a single element of Q'*B
            //  beyond the first N, which is initially zero.
            //
            qtbpj = 0.0;

            for (k = j; k < n; k++) {
                //
                //  Determine a Givens rotation which eliminates the
                //  appropriate element in the current row of D.
                //
                if (sdiag[k] != 0.0) {
                    if (fabs(r(k, k)) < fabs(sdiag[k])) {
                        cotan = r(k, k) / sdiag[k];
                        s = 0.5 / sqrt(0.25 + 0.25 * cotan * cotan);
                        c = s * cotan;
                    } else {
                        t = sdiag[k] / r(k, k);
                        c = 0.5 / sqrt(0.25 + 0.25 * t * t);
                        s = c * t;
                    }
                    //
                    //  Compute the modified diagonal element of R and
                    //  the modified element of (Q'*B,0).
                    //
                    r(k, k) = c * r(k, k) + s * sdiag[k];
                    temp = c * wa[k] + s * qtbpj;
                    qtbpj = -s * wa[k] + c * qtbpj;
                    wa[k] = temp;
                    //
                    //  Accumulate the tranformation in the row of S.
                    //
                    for (i = k + 1; i < n; i++) {
                        temp = c * r(i, k) + s * sdiag[i];
                        sdiag[i] = -s * r(i, k) + c * sdiag[i];
                        r(i, k) = temp;
                    }
                }
            }
        }
        //
        //  Store the diagonal element of S and restore
        //  the corresponding diagonal element of R.
        //
        sdiag[j] = r(j, j);
        r(j, j) = x[j];
    }
    //
    //  Solve the triangular system for Z.  If the system is
    //  singular, then obtain a least squares solution.
    //
    nsing = n;

    for (j = 0; j < n; j++) {
        if (sdiag[j] == 0.0 && nsing == n) {
            nsing = j;
        }

        if (nsing < n) {
            wa[j] = 0.0;
        }
    }

    for (j = nsing - 1; 0 <= j; j--) {
        sum2 = 0.0;
        for (i = j + 1; i < nsing; i++) {
            sum2 = sum2 + wa[i] * r(i, j);
        }
        wa[j] = (wa[j] - sum2) / sdiag[j];
    }
    //
    //  Permute the components of Z back to components of X.
    //
    for (j = 0; j < n; j++) {
        l = ipvt[j];
        x[l] = wa[j];
    }
}

}
}

#endif //LANRE_OPTIMIZE_QRSOLV_HPP
