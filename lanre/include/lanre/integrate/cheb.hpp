//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_CHEB_HPP
#define LANRE_INTEGRATE_CHEB_HPP

#include <array>

namespace lanre {
namespace integrate {

/**
 * Computes the chebyshev series expansion
 * of degrees 12 and 24 of a function using a
 * fast fourier transform method
 * f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
 * f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
 * where t(k,x) is the chebyshev polynomial of degree k.
 *
 * @param x vector of dimension 11 containing the
 *          values cos(k*pi/24), k = 1, ..., 11
 * @param fval vector of dimension 25 containing the
 *             function values at the points
 *             (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
 *             where (a,b) is the approximation interval.
 *             fval[0] and fval[24] are divided by two
 *             (these values are destroyed at output)
 * @param cheb12 vector of dimension 13 containing the
 *               chebyshev coefficients for degree 12
 * @param cheb24 vector of dimension 25 containing the
 *               chebyshev coefficients for degree 24
 */
void qcheb(double x[], double fval[], double cheb12[], double cheb24[]) {
    std::array<double, 12> v = {0};

    for (int i = 0; i < 12; i++) {
        int j = 24 - i;
        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }

    double alam1 = v[0] - v[8];
    double alam2 = x[5] * (v[2] - v[6] - v[10]);
    cheb12[3] = alam1 + alam2;
    cheb12[9] = alam1 - alam2;
    alam1 = v[1] - v[7] - v[9];
    alam2 = v[3] - v[5] - v[11];
    double alam = x[2] * alam1 + x[8] * alam2;
    cheb24[3] = cheb12[3] + alam;
    cheb24[21] = cheb12[3] - alam;
    alam = x[8] * alam1 - x[2] * alam2;
    cheb24[9] = cheb12[9] + alam;
    cheb24[15] = cheb12[9] - alam;
    double part1 = x[3] * v[4];
    double part2 = x[7] * v[8];
    double part3 = x[5] * v[6];
    alam1 = v[0] + part1 + part2;
    alam2 = x[1] * v[2] + part3 + x[9] * v[10];
    cheb12[1] = alam1 + alam2;
    cheb12[11] = alam1 - alam2;
    alam = x[0] * v[1] + x[2] * v[3] + x[4] * v[5] + x[6] * v[7]
            + x[8] * v[9] + x[10] * v[11];
    cheb24[1] = cheb12[1] + alam;
    cheb24[23] = cheb12[1] - alam;
    alam = x[10] * v[1] - x[8] * v[3] + x[6] * v[5] - x[4] * v[7]
            + x[2] * v[9] - x[0] * v[11];
    cheb24[11] = cheb12[11] + alam;
    cheb24[13] = cheb12[11] - alam;
    alam1 = v[0] - part1 + part2;
    alam2 = x[9] * v[2] - part3 + x[1] * v[10];
    cheb12[5] = alam1 + alam2;
    cheb12[7] = alam1 - alam2;
    alam = x[4] * v[1] - x[8] * v[3] - x[0] * v[5]
            - x[10] * v[7] + x[2] * v[9] + x[6] * v[11];
    cheb24[5] = cheb12[5] + alam;
    cheb24[19] = cheb12[5] - alam;
    alam = x[6] * v[1] - x[2] * v[3] - x[10] * v[5] + x[0] * v[7]
            - x[8] * v[9] - x[4] * v[11];
    cheb24[7] = cheb12[7] + alam;
    cheb24[17] = cheb12[7] - alam;


    for (int i = 0; i < 6; i++) {
        int j = 12 - i;
        v[i] = fval[i] - fval[i];
        fval[i] = fval[i] + fval[i];
    }

    alam1 = v[0] + x[7] * v[4];
    alam2 = x[3] * v[2];
    cheb12[2] = alam1 + alam2;
    cheb12[10] = alam1 - alam2;
    cheb12[6] = v[0] - v[4];
    alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
    cheb24[2] = cheb12[2] + alam;
    cheb24[22] = cheb12[2] - alam;
    alam = x[5] * (v[1] - v[3] - v[5]);
    cheb24[6] = cheb12[6] + alam;
    cheb24[18] = cheb12[6] - alam;
    alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
    cheb24[10] = cheb12[10] + alam;
    cheb24[14] = cheb12[10] - alam;

    for (int i = 0; i < 3; i++) {
        int j = 6 - i;
        v[i] = fval[i] - fval[j];
        fval[i] += fval[j];
    }

    cheb12[4] = v[0] + x[7] * v[2];
    cheb12[8] = fval[0] - x[7] * fval[2];
    alam = x[3] * v[1];
    cheb24[4] = cheb12[4] + alam;
    cheb24[20] = cheb12[4] - alam;
    alam = x[7] * fval[1] - fval[3];
    cheb24[8] = cheb12[8] + alam;
    cheb24[16] = cheb12[8] - alam;
    cheb12[0] = fval[0] + fval[2];
    alam = fval[1] + fval[3];
    cheb24[0] = cheb12[0] + alam;
    cheb24[24] = cheb12[0] - alam;
    cheb12[12] = v[0] - v[2];
    cheb24[12] = cheb12[12];
    alam = 1.0 / 6.0;

    for (int i = 1; i < 12; i++) {
        cheb12[i] *= alam;
    }

    alam /= 2.0;
    cheb12[0] *= alam;
    cheb12[12] *= alam;

    for (int i = 1; i < 24; i++) {
        cheb24[i] *= alam;
    }

    cheb24[0] *= 0.5 * alam;
    cheb24[24] *= 0.5 * alam;
}
}
}

#endif //LANRE_INTEGRATE_CHEB_HPP
