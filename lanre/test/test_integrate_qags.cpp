//
// Created by Logan Morrison on 3/30/20.
//

#include "lanre/integrate/qag.hpp"
#include "lanre/integrate/qags.hpp"
#include "lanre/integrate/qagi.hpp"
//#include "lanre/integrate/qagp.hpp"
//#include "lanre/integrate/qawc.hpp"
//#include "lanre/integrate/quad.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace lanre::integrate;

TEST(TestQAG, TestExp) {
    auto f = [](double x) {
        return exp(x);
    };

    double exact = exp(1.0) * (exp(4.0) - 1.0);

    double a = 1.0;
    double b = 5.0;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;

    std::cout << std::setprecision(16);
    for (int key = 1; key <= 6; key++) {
        auto start = std::chrono::high_resolution_clock::now();
        result = qag(f, a, b, epsabs, epsrel, key, &abserr);
        auto end = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;
        std::cout << "abserr = " << abserr << std::endl;

        double fracerr = std::abs(exact - result) / exact;
        ASSERT_LT(fracerr, epsrel);
    }
}

TEST(TestQAG, TestSin) {
    auto f = [](double x) {
        return sin(x);
    };

    double exact = 2.0;

    double a = 0.0;
    double b = M_PI;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;

    std::cout << std::setprecision(16);
    for (int key = 1; key <= 6; key++) {
        auto start = std::chrono::high_resolution_clock::now();
        result = qag(f, a, b, epsabs, epsrel, key, &abserr);
        auto end = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;

        double fracerr = std::abs(exact - result) / exact;
        ASSERT_LT(fracerr, epsrel);
    }
}

TEST(TestQAG, TestBesselK) {
    using lanre::special_functions::besselk3;
    auto f = [](double x) {
        return besselk3(x);
    };

    double exact = 2.828628119457949;

    double a = 1.0;
    double b = 10.0;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;

    std::cout << std::setprecision(16);
    for (int key = 1; key <= 6; key++) {
        auto start = std::chrono::high_resolution_clock::now();
        result = qag(f, a, b, epsabs, epsrel, key, &abserr);
        auto end = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;

        double fracerr = std::abs(exact - result) / exact;
        ASSERT_LT(fracerr, epsrel);
    }
}

TEST(TestQAGS, TestExp) {
    auto f = [](double x) {
        return exp(x);
    };

    double exact = exp(1.0) * (exp(4.0) - 1.0);

    double a = 1.0;
    double b = 5.0;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;

    std::cout << std::setprecision(16);
    for (int key = 1; key <= 6; key++) {
        auto start = std::chrono::high_resolution_clock::now();
        result = qags(f, a, b, epsabs, epsrel, &abserr);
        auto end = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;

        double fracerr = std::abs(exact - result) / exact;
        ASSERT_LT(fracerr, epsrel);
    }
}

TEST(TestQAGS, TestSin) {
    auto f = [](double x) {
        return sin(x);
    };

    double exact = 2.0;

    double a = 0.0;
    double b = M_PI;
    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double result, abserr;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = qags(f, a, b, epsabs, epsrel, &abserr);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken * 1e-3 << " (us)" << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, epsrel);
}

TEST(TestQAGS, TestBesselK) {
    using lanre::special_functions::besselk3;
    auto f = [](double x) {
        return besselk3(x);
    };

    double exact = 2.828628119457949;

    double a = 1.0;
    double b = 10.0;
    double epsabs = 1e-3;
    double epsrel = 1e-3;
    double result, abserr;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = qags(f, a, b, epsabs, epsrel, &abserr);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    double fracerr = std::abs(exact - result) / exact;


    std::cout << "res = " << result << ", time = " << time_taken * 1e-3 << " (us)" << std::endl;
    std::cout << "abserr = " << abserr << std::endl;
    std::cout << "fracerr = " << fracerr << std::endl;


    ASSERT_LT(fracerr, epsrel);
}

TEST(TestQAGI, TestExp) {
    auto f = [](double x) {
        return sin(x) / (x * x * x);
    };

    double exact = 0.3785300171241613;

    double bound = 1.0;
    double epsabs = 1e-8;
    double epsrel = 1e-8;
    double result, abserr;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = qagi(f, bound, 1, epsabs, epsrel, &abserr);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;
    std::cout << "abserr = " << abserr << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, epsrel);
}