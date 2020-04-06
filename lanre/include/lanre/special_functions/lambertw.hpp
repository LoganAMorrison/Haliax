// Created by Logan Morrison on 3/20/20.
// Implementation of product-log or Lambert-W function
// Taken from:
//      https://github.com/IstvanMezo/LambertW-function
//

#ifndef LANRE_SPECIAL_FUNCTION_LAMBERTW_HPP
#define LANRE_SPECIAL_FUNCTION_LAMBERTW_HPP

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

namespace lanre {
namespace special_functions {

template<class Type>
Type zexpz(Type z) {
    return z * exp(z);
}

//The derivative of z * exp(z) = exp(z) + z * exp(z)>
template<class Type>
Type zexpz_d(Type z) {
    return exp(z) + z * exp(z);
}

//The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
template<class Type>
Type zexpz_dd(Type z) {
    return 2.0 * exp(z) + z * exp(z);
}

//Determine the initial point for the root finding
std::complex<double> init_point(std::complex<double> z, int k) {
    std::complex<double> I{0.0, 1.0};
    std::complex<double> two_pi_k_I{0.0, 2.0 * M_PI * k};
    // initial point coming from the general asymptotic approximation
    std::complex<double> ip{log(z) + two_pi_k_I - log(log(z) + two_pi_k_I)};
    // used when we are close to the branch cut around zero and when k=0,-1
    std::complex<double> p{sqrt(2.0 * (M_E * z + 1.0))};

    //we are close to the branch cut, the initial point must be chosen carefully
    if (abs(z - (-exp(-1.0))) <= 1.0) {
        if (k == 0) {
            ip = -1.0 + p - 1.0 / 3.0 * pow(p, 2) + 11.0 / 72.0 * pow(p, 3);
        }
        if (k == 1 && z.imag() < 0.0) {
            ip = -1.0 - p - 1.0 / 3.0 * pow(p, 2) - 11.0 / 72.0 * pow(p, 3);
        }
        if (k == -1 && z.imag() > 0.0) {
            ip = -1.0 - p - 1.0 / 3.0 * pow(p, 2) - 11.0 / 72.0 * pow(p, 3);
        }
    }

    if (k == 0 && abs(z - 0.5) <= .5) {
        // (1,1) Pade approximant for W(0,a)
        ip = (0.35173371 * (0.1237166 + 7.061302897 * z)) /
                (2. + 0.827184 * (1. + 2. * z));
    }

    if (k == -1 && abs(z - .5) <= .5) {
        // (1,1) Pade approximant for W(-1,a)
        ip = -(((2.2591588985 + 4.22096 * I) * ((-14.073271 - 33.767687754 * I) * z - (12.7127 -
                19.071643 * I) * (1. + 2. * z))) /
                (2. - (17.23103 - 10.629721 * I) * (1. + 2. * z)));
    }

    return ip;
}

std::complex<double> lambertw(std::complex<double> z, int k = 0) {
    //For some particular z and k W(z,k) has simple value:
    if (z.real() == 0.0 && z.imag() == 0.0) {
        return (k == 0) ? 0.0 : -HUGE_VAL;
    }
    if ((z.real() == -exp(-1.0)) && (z.imag() == 0.0) && (k == 0 || k == -1)) {
        return -1.0;
    }
    if ((z.real() == exp(1.0)) && (z.imag() == 0.0) && k == 0) {
        return 1.0;
    }

    //Halley method begins

    // intermediate values in the Halley method
    std::complex<double> w{init_point(z, k)};
    std::complex<double> wprev{w};
    // max number of iterations. This eliminates improbable infinite loops
    const unsigned int maxiter = 30;
    // iteration counter
    unsigned int iter = 0;
    // difference threshold between the last two iteration results (or the iter number of iterations is taken)
    double prec = 1.0e-30;

    do {
        wprev = w;
        w -= 2.0 * ((zexpz(w) - z) * zexpz_d(w)) /
                (2.0 * pow(zexpz_d(w), 2) - (zexpz(w) - z) * zexpz_dd(w));
        iter++;
    } while ((abs(w - wprev) > prec) && iter < maxiter);
    return w;
}

std::complex<double> lambertw(double z, int k = 0) {
    return lambertw(std::complex<double>{z}, k);
}

std::vector<std::complex<double>> lambertw(std::vector<double> zs) {
    std::vector<std::complex<double>> ws(zs.size(), 0.0);
    for (size_t i = 0; i < zs.size(); i++) {
        ws[i] = lambertw(std::complex<double>{zs[i]}, 0);
    }
    return ws;
}

}
}


#endif //LANRE_SPECIAL_FUNCTION_LAMBERTW_HPP
