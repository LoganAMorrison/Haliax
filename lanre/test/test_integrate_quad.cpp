//
// Created by Logan Morrison on 3/30/20.
//

#include "lanre/integrate/quad.hpp"
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
    double abstol = 1e-8;
    double reltol = 1e-5;
    double result, error;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

TEST(TestQAG, TestSin) {
    auto f = [](double x) {
        return sin(x);
    };

    double exact = 2.0;

    double a = 0.0;
    double b = M_PI;
    double abstol = 1e-8;
    double reltol = 1e-5;
    double result, error;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

TEST(TestQAG, TestBesselK) {
    using lanre::special_functions::besselk3;
    auto f = [](double x) {
        return besselk3(x);
    };

    double exact = 2.828628119457949;

    double a = 1.0;
    double b = 10.0;
    double abstol = 1e-10;
    double reltol = 1e-10;
    double result, error;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken << " (ns)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

TEST(TestQAGI, TestBesselK0) {
    using lanre::special_functions::besselk0;
    auto f = [](double x) {
        return besselk0(x);
    };

    double exact = M_PI / 2;

    double a = 0.0;
    double b = std::numeric_limits<double>::infinity();
    double abstol = 1e-5;
    double reltol = 1e-5;
    double result, error;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken * 1e-3 << " (us)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

TEST(TestQAGI, TestGaussian) {
    using lanre::special_functions::besselk0;
    auto f = [](double x) {
        return exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
    };

    double exact = 1.0;

    double a = -std::numeric_limits<double>::infinity();
    double b = std::numeric_limits<double>::infinity();
    double abstol = 1e-5;
    double reltol = 1e-5;
    double result, error;

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken * 1e-3 << " (us)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

TEST(TestQAGP, TestGaussian) {
    using lanre::special_functions::besselk0;
    auto f = [](double x) {
        return pow(x, 3.0) * log(fabs((x * x - 1.0) * (x * x - 2.0)));
    };

    double exact = 52.740748383471444993;

    double a = 0.0;
    double b = 3.0;
    double abstol = 1e-5;
    double reltol = 1e-5;
    double result, error;

    std::vector<double> pts = {1.0, sqrt(2)};

    std::cout << std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    result = Quad<double>::integrate(f, a, b, pts, abstol, reltol, 500, &error);
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "res = " << result << ", time = " << time_taken * 1e-3 << " (us)" << std::endl;
    std::cout << "abserr = " << error << std::endl;

    double fracerr = std::abs(exact - result) / exact;
    ASSERT_LT(fracerr, reltol);
}

