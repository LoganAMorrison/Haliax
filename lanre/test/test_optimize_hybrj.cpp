//
// Created by Logan Morrison on 4/11/20.
//

#include <lanre/lanre.hpp>
#include <lanre/optimize/hybrj.hpp>
#include <gtest/gtest.h>

using namespace lanre;
using namespace optimize;


TEST(TestHybrj, TestRosenbrock) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, Matrix<double> &jac, int &iflag) -> void {
        const double xx = x(0);
        const double yy = x(1);

        fvec[0] = 1.0 - xx;
        fvec[1] = 10.0 * pow(yy - xx * xx, 2);

        jac(0, 0) = -1.0;
        jac(0, 1) = 0.0;
        jac(1, 0) = -40.0 * xx * (yy - xx * xx);
        jac(1, 1) = 20.0 * (yy - xx * xx);
    };

    Vector<double> x = Vector<double>::Zero(2);
    Vector<double> fvec = Vector<double>::Zero(2);

    x(0) = 1.2;
    x(1) = 1.0;


    double tol = sqrt(2.22044604926e-16);
    int info = hybrj(f, x, fvec, tol);
    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
    }
}

TEST(TestHybrj, TestPowellSingular) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, Matrix<double> &jac, int &iflag) -> void {
        fvec[0] = x[0] + 10.0 * x[1];
        fvec[1] = sqrt(5.0) * (x[2] - x[3]);
        fvec[2] = pow(x[1] - 2.0 * x[2], 2);
        fvec[3] = sqrt(10.0) * pow(x[0] - x[3], 2);

        jac(0, 0) = 1.0;
        jac(0, 1) = 10.0;
        jac(0, 2) = 0.0;
        jac(0, 3) = 0.0;

        jac(1, 0) = 0.0;
        jac(1, 1) = 0.0;
        jac(1, 2) = sqrt(5);
        jac(1, 3) = -sqrt(5);

        jac(2, 0) = 0.0;
        jac(2, 1) = 2.0 * (x[1] - 2.0 * x[2]);
        jac(2, 2) = -4.0 * (x[1] - 2.0 * x[2]);
        jac(2, 3) = 0.0;

        jac(3, 0) = 2.0 * sqrt(10.0) * (x[0] - x[3]);
        jac(3, 1) = 0.0;
        jac(3, 2) = 0.0;
        jac(3, 3) = -2.0 * sqrt(10.0) * (x[0] - x[3]);
    };

    Vector<double> x(4);
    x << 3.0, -1.0, 0.0, 1.0;
    Vector<double> fvec(x.size());

    std::cout << "Starting vector: " << x << std::endl;

    double tol = sqrt(2.22044604926e-16);
    int info = hybrj(f, x, fvec, tol);

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
    }
}

TEST(TestHybrj, TestPowellBadlyScaled) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, Matrix<double> &jac, int &iflag) -> void {
        fvec[0] = 1.0e4 * x[0] * x[1] - 1.0;
        fvec[1] = exp(-x[0]) + exp(-x[1]) - 1.0001e0;

        jac(0, 0) = 1.0e4 * x[1];
        jac(0, 1) = 1.0e4 * x[1];
        jac(1, 0) = -x[0] * exp(-x[0]);
        jac(1, 1) = -x[1] * exp(-x[1]);

    };

    Vector<double> x(2);
    x << 0.0, 1.0;
    Vector<double> fvec(2);

    double tol = sqrt(2.22044604926e-16);
    int info = hybrj(f, x, fvec, tol);

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
    }
}