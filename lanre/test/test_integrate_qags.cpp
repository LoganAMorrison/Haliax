//
// Created by Logan Morrison on 3/30/20.
//

#include "lanre/integrate/qag.hpp"
#include "lanre/integrate/qagi.hpp"
#include "lanre/integrate/qagp.hpp"
#include "lanre/integrate/qags.hpp"
#include "lanre/integrate/qawc.hpp"
#include "lanre/integrate/quad.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace lanre::integrate;

TEST(TestQAG, TestExp) {
    auto f = [](double x) {
        return exp(x);
    };

    double a = 1.0;
    double b = 5.0;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;
    int neval, ier;
    int limit = 500;
    int lenw = limit * 4;
    int last;
    auto iwork = new int[limit];
    auto work = new double[limit];

    std::cout << std::setprecision(30) << std::endl;
    for (int key = 1; key <= 6; key++) {
        result = qag(f, a, b, epsabs, epsrel, key, &abserr, &neval, &ier, &last);
        std::cout << "res = " << result << std::endl;
    }

    delete[] iwork;
    delete[] work;
}

TEST(TestQAGS, TestExp) {
    auto f = [](double x) {
        return exp(x);
    };

    double a = 1.0;
    double b = 5.0;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;
    int neval, ier;

    result = qags(f, a, b, epsabs, epsrel, &abserr, &neval, &ier);

    std::cout << "res = " << result << std::endl;
}

TEST(TestQAGP, TestExp) {
    auto f = [](double x) {
        return exp(x);
    };

    double a = 1.0;
    double b = 5.0;
    int npts2 = 2;
    auto points = new double[npts2];
    double epsabs = 1e-8;
    double epsrel = 1e-8;
    double result, abserr;
    int neval, ier;
    int leniw = 2 * 500 + npts2;

    result = qagp(f, a, b, npts2, points, epsabs, epsrel, &abserr, &neval, &ier);

    std::cout << "res = " << result << std::endl;
}

TEST(TestQAGI, TestExp) {
    auto f = [](double x) {
        return sin(x) * exp(-x * x * x);
    };

    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;
    int neval, ier;
    int last;
    double bound = 0.0;
    int inf = 1;

    result = qagi(f, bound, inf, epsabs, epsrel, &abserr, &neval, &ier);

    std::cout << "res = " << result << std::endl;

}