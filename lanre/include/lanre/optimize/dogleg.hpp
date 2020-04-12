//
// Created by Logan Morrison on 4/10/20.
//

#ifndef LANRE_OPTIMIZE_DOGLEG_HPP
#define LANRE_OPTIMIZE_DOGLEG_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

void dogleg(
        Vector<double> &r,
        const Vector<double> &diag,
        Vector<double> &qtb,
        double delta,
        Vector<double> &x,
        Vector<double> &wa1,
        Vector<double> &wa2
) {
    const int n = x.size();
    static const double epsmch = std::numeric_limits<double>::epsilon();

    int jj = (n * (n + 1)) / 2 + 1;

    for (int k = 1; k <= n; k++) {
        int j = n - k + 1;
        int jp1 = j + 1;
        jj = jj - k;
        int l = jj + 1;
        double sum = 0.0;
        for (int i = jp1; i <= n; i++) {
            sum = sum + r[l - 1] * x[i - 1];
            l = l + 1;
        }
        double temp = r[jj - 1];
        if (temp == 0.0) {
            l = j;
            for (int i = 1; i <= j; i++) {
                temp = std::max(temp, fabs(r[l - 1]));
                l = l + n - i;
            }
            temp = epsmch * temp;
            if (temp == 0.0) {
                temp = epsmch;
            }
        }
        x[j - 1] = (qtb[j - 1] - sum) / temp;
    }
    //
    //  Test whether the Gauss-Newton direction is acceptable.
    //
    for (int j = 0; j < n; j++) {
        wa1[j] = 0.0;
        wa2[j] = diag[j] * x[j];
    }
    double qnorm = wa2.norm();

    if (qnorm <= delta) {
        return;
    }
    //
    //  The Gauss-Newton direction is not acceptable.
    //  Calculate the scaled gradient direction.
    //
    int l = 0;
    for (int j = 0; j < n; j++) {
        double temp = qtb[j];
        for (int i = j; i < n; i++) {
            wa1[i] = wa1[i] + r[l] * temp;
            l++;
        }
        wa1[j] = wa1[j] / diag[j];
    }
    //
    //  Calculate the norm of the scaled gradient and test for
    //  the special case in which the scaled gradient is zero.
    //
    double gnorm = wa1.norm();
    double sgnorm = 0.0;
    double alpha = delta / qnorm;
    //
    //  Calculate the point along the scaled gradient
    //  at which the quadratic is minimized.
    //
    if (gnorm != 0.0) {
        for (int j = 0; j < n; j++) {
            wa1[j] = (wa1[j] / gnorm) / diag[j];
        }
        l = 0;
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = j; i < n; i++) {
                sum = sum + r[l] * wa1[i];
                l = l + 1;
            }
            wa2[j] = sum;
        }
        double temp = wa2.norm();
        sgnorm = (gnorm / temp) / temp;
        alpha = 0.0;
        //
        //  If the scaled gradient direction is not acceptable,
        //  calculate the point along the dogleg at which the quadratic is minimized.
        //
        if (sgnorm < delta) {
            double bnorm = qtb.norm();
            temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
            temp = temp - (delta / qnorm) * (sgnorm / delta) * (sgnorm / delta)
                    + sqrt(pow(temp - (delta / qnorm), 2)
                                   + (1.0 - (delta / qnorm) * (delta / qnorm))
                    * (1.0 - (sgnorm / delta) * (sgnorm / delta)));
            alpha = ((delta / qnorm)
                    * (1.0 - (sgnorm / delta) * (sgnorm / delta))) / temp;
        }
    }
    //
    //  Form appropriate convex combination of the Gauss-Newton
    //  direction and the scaled gradient direction.
    //
    double temp = (1.0 - alpha) * std::min(sgnorm, delta);
    for (int j = 0; j < n; j++) {
        x[j] = temp * wa1[j] + alpha * x[j];
    }
}


}
}

#endif //LANRE_OPTIMIZE_DOGLEG_HPP
