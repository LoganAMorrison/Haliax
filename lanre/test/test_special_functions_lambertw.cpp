//
// Created by Logan Morrison on 3/20/20.
//

#include "lanre/special_functions/lambertw.hpp"
#include "lanre/autodiff/dual.hpp"
#include <gtest/gtest.h>
#include <complex>
#include <boost/math/special_functions/lambert_w.hpp>

using lanre::special_functions::lambertw;
using lanre::autodiff::Dual;
using boost::math::lambert_w0;
static const double inf = std::numeric_limits<double>::infinity();

TEST(TestLambertW, TestBranch0) {
    std::pair<double, double> minmax = std::make_pair(1e-5, 1e10);
    size_t num_tests = 100;
    double step = (minmax.second - minmax.first) / double(num_tests - 1);
    for (int i = 0; i < num_tests; i++) {
        double xx = minmax.first + i * step;
        auto lw = lambertw(xx).real();
        double lwboost = lambert_w0(xx);
        std::cout << xx << std::endl;
        std::cout << lw << std::endl;
        std::cout << lwboost << std::endl << std::endl;
        ASSERT_LT(std::abs(1 - lw / lwboost), 1e-3);
    }
}

TEST(TestLambertW, TestBranch0Dual) {
    std::pair<double, double> minmax = std::make_pair(1e-5, 1e10);
    size_t num_tests = 100;
    double step = (minmax.second - minmax.first) / double(num_tests - 1);

    for (int i = 0; i < num_tests; i++) {
        auto xx = Dual<double>(minmax.first + i * step, 1.0);
        auto lw = lambertw(xx).real();
        double lwboost = lambert_w0(xx.val);
        double deriv = lwboost / (xx.val * (lwboost + 1.0));
        std::cout << xx << std::endl;
        std::cout << lw << std::endl;
        std::cout << lwboost << ", " << deriv << std::endl << std::endl;
        //ASSERT_LT(std::abs(1 - lw / lwboost), 1e-3);
    }
}
