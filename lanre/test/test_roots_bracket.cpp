//
// Created by Logan Morrison on 3/20/20.
//

#include "lanre/roots.hpp"
#include <gtest/gtest.h>

using namespace lanre;

struct Quadratic {
    double a;
    double b;
    double c;

    double operator()(double x) {
        return a * x * x + b * x + c;
    }
};

TEST(RootTest, QuadraticRootTest) {

    Quadratic func{1.0, 1.0, 0.0};

    double root = brent(func, -2.0, -0.5, 1e-5);

    std::cout << "root = " << root << std::endl;
}