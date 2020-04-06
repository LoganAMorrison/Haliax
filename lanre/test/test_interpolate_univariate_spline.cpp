//
// Created by Logan Morrison on 4/3/20.
//

#include "lanre/interpolate/univariate_spline.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace lanre::interpolate;

/**
 * Test the UnivariateSpline class on sin(x) on the interval from (0,pi)
 */
TEST(TestUnivariateSpline, TestSin) {
    double frac_diff;

    const size_t num_pts = 50;
    const std::pair<double, double> interval{0.0, M_PI};

    const double step = (interval.second - interval.first) / double(num_pts - 1);
    std::vector<double> abscissas(num_pts);
    std::vector<double> ordinates(num_pts);
    std::vector<double> weights(num_pts);

    // Fill the abscissas with 50 values ranging from (0,pi) and the ordinates
    // with the sine evaluated at the corresponding abscissa. Set all weights
    // to unity
    for (size_t i = 0; i < num_pts; i++) {
        abscissas[i] = interval.first + i * step;
        ordinates[i] = sin(abscissas[i]);
        weights[i] = 1.0;
    }

    // Boundary box is set to endpoints
    std::pair<double, double> bbox = std::make_pair(abscissas.front(), abscissas.back());
    // Make a cubic spline
    int k = 3;
    // No smoothing
    double s = 0.0;
    // Set extrapolation to result in zero.
    int ext = 0;

    // Construct spline
    auto spline = UnivariateSpline(abscissas, ordinates, weights, bbox, k, s, ext);
    // Evaluate spline at abscissas
    auto ss = spline(abscissas);

    // Check that the spline evaluated at abscissas is close (within 1%) to ordinates
    for (size_t i = 0; i < num_pts; i++) {
        double spval = spline(abscissas[i]);
        if (ordinates[i] == 0.0) {
            ASSERT_LT(fabs(spval), 1e-3);
        } else {
            frac_diff = fabs(spval - ordinates[i]) / fabs(ordinates[i]);
            ASSERT_LT(frac_diff, 1e-2);
        }
    }

    // Construct new points which don't coincide with abscissas
    const size_t num_new_pts = 199;
    const double step2 = (interval.second - interval.first) / double(num_new_pts - 1);

    // Check that the spline evaluated at new points is close (within 1%) to sin(x)
    for (size_t i = 0; i < num_new_pts; i++) {
        double xx = interval.first + i * step2;
        double spval = spline(xx);
        double sinx = sin(xx);
        if (sinx == 0.0) {
            ASSERT_LT(fabs(spval), 1e-3);
        } else {
            frac_diff = fabs(spval - sinx) / fabs(sinx);
            ASSERT_LT(frac_diff, 1e-2);
        }
    }

    // Check the integral of sin(x) over (0, pi)
    double integral = spline.integrate(0.0, M_PI);
    frac_diff = fabs(integral - 2.0) / 2.0;
    ASSERT_LT(frac_diff, 1e-3);

    // Check the derivatives of spline
    for (size_t i = 0; i < num_new_pts; i++) {
        double xx = interval.first + i * step2;
        double cosx = cos(xx);
        double sinx = sin(xx);

        // Check first derivative: d/dx sin = cos
        double dspval = spline.derivative(xx, 1);
        if (fabs(cosx) <= 10 * std::numeric_limits<double>::epsilon()) {
            ASSERT_LT(fabs(dspval), 1e-3);
        } else {
            frac_diff = fabs(dspval - (cosx)) / fabs(cosx);
            ASSERT_LT(frac_diff, 1e-2);
        }
        // Check second derivative: d^2/dx^2 sin = -sin
        double d2spval = spline.derivative(xx, 2);
        if (abs(sinx) <= 10 * std::numeric_limits<double>::epsilon()) {
            ASSERT_LT(fabs(d2spval), 1e-3);
        } else {
            frac_diff = fabs(d2spval - (-sinx)) / fabs(sinx);
            ASSERT_LT(frac_diff, 1e-2);
        }
    }
}

/**
 * Test the UnivariateSpline class on exp(x) on the interval from (0,4)
 */
TEST(TestUnivariateSpline, TestExp) {
    double frac_diff;

    const size_t num_pts = 50;
    const std::pair<double, double> interval{0.0, 4.0};

    const double step = (interval.second - interval.first) / double(num_pts - 1);
    std::vector<double> abscissas(num_pts);
    std::vector<double> ordinates(num_pts);
    std::vector<double> weights(num_pts);

    // Fill the abscissas with 50 values ranging from (0,pi) and the ordinates
    // with the sine evaluated at the corresponding abscissa. Set all weights
    // to unity
    for (size_t i = 0; i < num_pts; i++) {
        abscissas[i] = interval.first + i * step;
        ordinates[i] = exp(abscissas[i]);
        weights[i] = 1.0;
    }

    // Boundary box is set to endpoints
    std::pair<double, double> bbox = std::make_pair(abscissas.front(), abscissas.back());
    // Make a cubic spline
    int k = 3;
    // No smoothing
    double s = 0.0;
    // Set extrapolation to result in zero.
    int ext = 0;

    // Construct spline
    auto spline = UnivariateSpline(abscissas, ordinates, weights, bbox, k, s, ext);
    // Evaluate spline at abscissas
    auto ss = spline(abscissas);

    // Check that the spline evaluated at abscissas is close (within 1%) to ordinates
    for (size_t i = 0; i < num_pts; i++) {
        double spval = spline(abscissas[i]);
        if (ordinates[i] == 0.0) {
            ASSERT_LT(fabs(spval), 1e-3);
        } else {
            frac_diff = fabs(spval - ordinates[i]) / fabs(ordinates[i]);
            ASSERT_LT(frac_diff, 1e-2);
        }
    }

    // Construct new points which don't coincide with abscissas
    const size_t num_new_pts = 199;
    const double step2 = (interval.second - interval.first) / double(num_new_pts - 1);

    // Check that the spline evaluated at new points is close (within 1%) to sin(x)
    for (size_t i = 0; i < num_new_pts; i++) {
        double xx = interval.first + i * step2;
        double spval = spline(xx);
        double expx = exp(xx);

        frac_diff = fabs(spval - expx) / fabs(expx);
        ASSERT_LT(frac_diff, 1e-2);
    }

    // Check the integral of sin(x) over (0, pi)
    double integral = spline.integrate(interval.first, interval.second);
    double exact = exp(4.0) - 1.0;
    frac_diff = fabs(integral - exact) / exact;
    ASSERT_LT(frac_diff, 1e-3);

    // Check the derivatives of spline
    for (size_t i = 0; i < num_new_pts; i++) {
        double xx = interval.first + i * step2;
        double expx = exp(xx);

        // Check first derivative
        double dspval = spline.derivative(xx, 1);
        frac_diff = fabs(dspval - (expx)) / fabs(expx);
        ASSERT_LT(frac_diff, 1e-2);
        // Check second derivative
        double d2spval = spline.derivative(xx, 2);
        frac_diff = fabs(dspval - (expx)) / fabs(expx);
        ASSERT_LT(frac_diff, 1e-2);
    }
}