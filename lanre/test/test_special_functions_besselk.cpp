//
// Created by Logan Morrison on 3/25/20.
//
#include "lanre/special_functions/besselk.hpp"
#include "lanre/autodiff/dual.hpp"
#include "gtest/gtest.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <iomanip>

using namespace lanre::special_functions;
using namespace lanre::autodiff;
using boost::math::cyl_bessel_k;
using boost::math::cyl_bessel_k_prime;

TEST(BesselKTest, BesselK0Test) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    for (int i = 0; i < 100; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double lan = besselk0(x);
        double boost = cyl_bessel_k(0, x);
        double frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }

    for (int i = 0; i < 100; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double lan = besselk0e(x);
        double boost = cyl_bessel_k(0, x) * exp(x);
        double frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }
}

TEST(BesselKTest, BesselK1Test) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    for (int i = 0; i < 100; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double lan = besselk1(x);
        double boost = cyl_bessel_k(1, x);
        double frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }

    for (int i = 0; i < 100; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double lan = besselk1e(x);
        double boost = cyl_bessel_k(1, x) * exp(x);
        double frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }
}

TEST(BesselKTest, BesselKNTest) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    double x, lan, boost, frac_diff;
    for (int n = 2; n <= 10; n++) {
        for (int i = 0; i < 100; i++) {
            x = pow(10.0, logxmin + i * logxstep);

            lan = besselkn(n, x);
            boost = cyl_bessel_k(n, x);
            frac_diff = fabs(boost - lan) / fabs(boost);
            ASSERT_LT(frac_diff, 1e-3);

            lan = besselkne(n, x);
            boost = cyl_bessel_k(n, x) * exp(x);
            frac_diff = fabs(boost - lan) / fabs(boost);
            ASSERT_LT(frac_diff, 1e-3);
        }
    }
}

TEST(BesselKTest, BesselK0DerivTest) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    Dual<double> z;
    double x, lan, boost, frac_diff;
    for (int i = 0; i < 100; i++) {
        z = Dual<double>(pow(10.0, logxmin + i * logxstep), 1.0);
        x = z.val;
        lan = besselk0(z).eps;
        boost = cyl_bessel_k_prime(0, x);
        //std::cout << "lan = " << lan << std::endl;
        //std::cout << "boost = " << boost << std::endl << std::endl;
        frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }
}

TEST(BesselKTest, BesselK1DerivTest) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    Dual<double> z;
    double x, lan, boost, frac_diff;
    for (int i = 0; i < 100; i++) {
        z = Dual<double>(pow(10.0, logxmin + i * logxstep), 1.0);
        x = z.val;
        lan = besselk1(z).eps;
        boost = cyl_bessel_k_prime(1, x);
        //std::cout << "lan = " << lan << std::endl;
        //std::cout << "boost = " << boost << std::endl << std::endl;
        frac_diff = fabs(boost - lan) / fabs(boost);
        ASSERT_LT(frac_diff, 1e-3);
    }
}

TEST(BesselKTest, BesselKNDerivTest) {
    double logxmin = -3;
    double logxmax = 2;
    double logxstep = (logxmax - logxmin) / 99.0;

    Dual<double> z;
    double x, lan, boost, frac_diff;
    for (int n = 2; n <= 10; n++) {
        for (int i = 0; i < 100; i++) {
            z = Dual<double>(pow(10.0, logxmin + i * logxstep), 1.0);
            x = z.val;
            lan = besselkn(n, z).eps;
            boost = cyl_bessel_k_prime(n, x);
            //std::cout << "lan = " << lan << std::endl;
            //std::cout << "boost = " << boost << std::endl << std::endl;
            frac_diff = fabs(boost - lan) / fabs(boost);
            ASSERT_LT(frac_diff, 1e-3);
        }
    }
}
