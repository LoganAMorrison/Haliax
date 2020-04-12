//
// Created by Logan Morrison on 4/10/20.
//

#include <lanre/lanre.hpp>
#include <lanre/optimize/hybrd.hpp>
#include <vector>
#include <gtest/gtest.h>

using namespace lanre;
using namespace optimize;


TEST(TestHybrd, TestRosenbrock) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        fvec[0] = 1.0 - x[0];
        fvec[1] = 10.0 * pow(x[1] - x[0] * x[0], 2);
    };

    Vector<double> x = Vector<double>::Zero(2);
    Vector<double> fvec = Vector<double>::Zero(2);

    x(0) = 1.2;
    x(1) = 1.0;


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
    }

    double norm = 0.0;
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}

TEST(TestHybrd, TestPowellSingular) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        fvec[0] = x[0] + 10.0 * x[1];
        fvec[1] = sqrt(5.0) * (x[2] - x[3]);
        fvec[2] = pow(x[1] - 2.0 * x[2], 2);
        fvec[3] = sqrt(10.0) * pow(x[0] - x[3], 2);
    };

    Vector<double> x(4);
    x << 3.0, -1.0, 0.0, 1.0;
    Vector<double> fvec(x.size());

    std::cout << "Starting vector: " << x << std::endl;

    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    double norm = 0.0;
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}

TEST(TestHybrd, TestPowellBadlyScaled) {
    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        fvec[0] = 1.0e4 * x[0] * x[1] - 1.0;
        fvec[1] = exp(-x[0]) + exp(-x[1]) - 1.0001e0;
    };

    Vector<double> x(2);
    x << 0.0, 1.0;
    Vector<double> fvec(2);

    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    double norm = 0.0;
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}

TEST(TestHybrd, TestWoodFunction) {
    const double c3 = 2.0e2;
    const double c4 = 2.02e1;
    const double c5 = 1.98e1;
    const double c6 = 1.8e2;

    auto f = [c3, c4, c5, c6](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        double temp1 = x[1] - x[0] * x[0];
        double temp2 = x[3] - x[2] * x[2];
        fvec[0] = -c3 * x[0] * temp1 - (1.0 - x[0]);
        fvec[1] = c3 * temp1 + c4 * (x[1] - 1.0) + c5 * (x[3] - 1.0);
        fvec[2] = -c6 * x[2] * temp2 - (1.0 - x[2]);
        fvec[3] = c6 * temp2 + c4 * (x[3] - 1.0) + c5 * (x[1] - 1.0);
    };

    Vector<double> x(4);
    x << -3.0, -1.0, -3.0, -1.0;
    Vector<double> fvec(4);


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}

TEST(TestHybrd, TestHelicalValley) {
    const double c7 = 2.5e-1;
    const double c8 = 5.0e-1;

    auto f = [c7, c8](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        double tpi = 8.0 * atan(1.0);
        double temp1 = std::copysign(c7, x[1]);
        if (x[0] > 0.0) {
            temp1 = atan(x[1] / x[0]) / tpi;
        }
        if (x[0] < 0.0) {
            temp1 = atan(x[1] / x[0]) / tpi + c8;
        }
        double temp2 = sqrt(x[0] * x[0] + x[1] * x[1]);
        fvec[0] = 10.0 * (x[2] - 10.0 * temp1);
        fvec[1] = 10.0 * (temp2 - 1.0);
        fvec[2] = x[2];
    };

    Vector<double> x(3);
    x << -1.0, 0.0, 0.0;
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}

TEST(TestHybrd, TestWatson) {
    const double c9 = 2.9e1;
    auto f = [c9](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {
        const int n = x.size();
        int k, i, j;
        double ti, sum1, sum2, temp, temp1, temp2;
        for (k = 0; k < n; k++) {
            fvec[k] = 0.0;
        }
        for (i = 1; i <= 29; i++) {
            ti = double(i) / c9;
            sum1 = 0.0;
            temp = 1.0;
            for (j = 2; j <= n; j++) {
                sum1 += double(j - 1) * temp * x[j - 1];
                temp = ti * temp;
            }
            sum2 = 0.0;
            temp = 1.0;
            for (j = 1; j <= n; j++) {
                sum2 += temp * x[j - 1];
                temp = ti * temp;
            }
            temp1 = sum1 - sum2 * sum2 - 1.0;
            temp2 = 2.0 * ti * sum2;
            temp = 1.0 / ti;
            for (k = 1; k <= n; k++) {
                fvec[k - 1] += temp * (double(k - 1) - temp2) * temp1;
                temp = ti * temp;
            }
        }
        temp = x[1] - x[0] * x[0] - 1.0;
        fvec[0] += x[0] * (1.0 - 2.0 * temp);
        fvec[1] += temp;
    };
    int n = 4;
    Vector<double> x = Vector<double>::Zero(n);
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
    ASSERT_LT(norm, tol);
    ASSERT_EQ(info, 1);
}


TEST(TestHybrd, TestChebyQuad) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };
    int n = 10;
    Vector<double> x = Vector<double>::Zero(n);
    Vector<double> fvec(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        x[i] = double(i) / double(n + 1);
    }


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestBrownAlmostLinear) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };

    int n = 10;
    Vector<double> x = Vector<double>::Constant(n, 0.5);
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestDiscreteBoundaryValueAndIntegralEquation) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };
    int n = 10;
    Vector<double> x = Vector<double>::Zero(n);
    Vector<double> fvec(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        double ti = double(i) / double(n + 1);
        x[i] = ti * (ti - 1.0);
    }


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestTrig) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };
    int n = 10;
    Vector<double> x = Vector<double>::Constant(n, 1.0 / double(n));
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestVariableDimensioned) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };
    int n = 10;
    Vector<double> x = Vector<double>::Zero(n);
    Vector<double> fvec(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        x[i] = 1.0 - double(i) / double(n);
    }


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestBroydenTridiagonalBanded) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };
    int n = 10;
    Vector<double> x = Vector<double>::Constant(n, -1.0);
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}

TEST(TestHybrd, TestMultipleInitInput) {

    auto f = [](Vector<double> &x, Vector<double> &fvec, int &iflag) -> void {

    };

    int n = 10;
    Vector<double> x = Vector<double>::Constant(n, 0.0);
    Vector<double> fvec(x.size());


    double tol = sqrt(2.22044604926e-16);
    int info = hybrid(f, x, fvec, tol);
    double norm = 0.0;

    std::cout << "info = " << info << std::endl;

    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    for (size_t i = 0; i < fvec.size(); i++) {
        std::cout << "f[" << i << "] = " << fvec[i] << std::endl;
        norm += fvec[i] * fvec[i];
    }
    std::cout << "norm(f) = " << sqrt(norm) << std::endl;
}






