// Created by Logan Morrison on 4/3/20.
// Mimics the SciPy class scipy.interpolate.UnivariateSpline

#ifndef LANRE_INTERPOLATE_INTERP1D_HPP
#define LANRE_INTERPOLATE_INTERP1D_HPP

#include "lanre/interpolate/curfit.hpp"
#include "lanre/interpolate/spalde.hpp"
#include "lanre/interpolate/splder.hpp"
#include "lanre/interpolate/splev.hpp"
#include "lanre/interpolate/splint.hpp"
#include "lanre/interpolate/sproot.hpp"
#include "lanre/autodiff/dual.hpp"

#include <cmath>
#include <utility>
#include <vector>
#include <iostream>

namespace lanre {
namespace interpolate {

class UnivariateSpline {
private:
    /* 1-D array of independent input data. Must be increasing;
     * must be strictly increasing if `s` is 0.
     */
    std::vector<double> m_x;

    /* 1-D array of dependent input data, of the same length as `x`. */
    std::vector<double> m_y;

    /* Weights for spline fitting.  Must be positive.  If None (default),
     * weights are all equal.
     */
    std::vector<double> m_w;

    /* 2-sequence specifying the boundary of the approximation interval. If
     * None (default), ``bbox=[x[0], x[-1]]``.
     */
    std::pair<double, double> m_bbox;

    /*Degree of the smoothing spline.  Must be <= 5.
     * Default is k=3, a cubic spline.
     */
    int m_k = 3;

    /* Positive smoothing factor used to choose the number of knots.  Number
     * of knots will be increased until the smoothing condition is satisfied::
     * sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s
     * If None (default), ``s = len(w)`` which should be a good value if
     * ``1/w[i]`` is an estimate of the standard deviation of ``y[i]``.
     * If 0, spline will interpolate through all data points.
     */
    double m_s;

    /* Controls the extrapolation mode for elements
     * not in the interval defined by the knot sequence.
     */
    int m_ext = 0;

    int m_n; // Number of knots
    int m_m; // size of x, y, w
    int m_nest; // size of t, c, fpint, z
    int m_k1; // width of a and q
    int m_k2; // width of b and g
    std::vector<double> m_t;
    std::vector<double> m_c;
    double *m_fpint;
    double *m_z;
    double **m_a;
    double **m_b;
    double **m_g;
    double **m_q;
    int *m_iwrk;
public:
    UnivariateSpline(
            std::vector<double> x,
            std::vector<double> y,
            std::vector<double> w,
            std::pair<double, double> bbox,
            int k,
            double s,
            int ext
    ) : m_x(std::move(x)), m_y(std::move(y)), m_w(std::move(w)),
        m_bbox(std::move(bbox)), m_k(k), m_s(s), m_ext(ext) {
        if (x.size() != y.size() || x.size() != w.size() || y.size() != w.size()) {
            throw std::runtime_error("UnivariateSpline: x, y and w must be of same size");
        }
        m_m = m_x.size();

        m_nest = m_m + m_k + 1;
        m_k1 = m_k + 1;
        m_k2 = m_k1 + 1;

        m_t = std::vector<double>(m_nest, 0.0);
        m_c = std::vector<double>(m_nest, 0.0);
        m_fpint = new double[m_nest];
        m_z = new double[m_nest];
        m_a = new double *[m_nest];
        m_b = new double *[m_nest];
        m_g = new double *[m_nest];
        m_q = new double *[m_m];
        m_iwrk = new int[m_nest];

        for (int i = 0; i < m_nest; i++) {
            m_a[i] = new double[m_k1];
            m_b[i] = new double[m_k2];
            m_g[i] = new double[m_k2];
        }
        for (int i = 0; i < m_m; i++) {
            m_q[i] = new double[m_k1];
        }

        // fit the spline
        int ier;
        int iopt = 0;
        double fp;

        dierckx::curfit(
                iopt,
                m_m,
                m_x.data(),
                m_y.data(),
                m_w.data(),
                m_bbox.first,
                m_bbox.second,
                m_k,
                m_s,
                m_nest,
                m_n,
                m_t.data(),
                m_c.data(),
                fp,
                m_fpint,
                m_z,
                m_a,
                m_b,
                m_g,
                m_q,
                m_iwrk,
                ier
        );

        if (ier > 0) {
            if (ier == 1) {
                throw std::runtime_error(
                        "UnivariateSpline: the required storage space exceeds the available\n"
                        "storage space, as specified by the parameter nest.\n"
                        "probably causes : nest too small. if nest is already\n"
                        "large (say nest > m/2), it may also indicate that s is\n"
                        "too small\n"
                        "the approximation returned is the weighted least-squares\n"
                        "spline according to the knots t(1),t(2),...,t(n). (n=nest)\n"
                        "the parameter fp gives the corresponding weighted sum of\n"
                        "squared residuals (fp>s).");
            } else if (ier == 2) {
                throw std::runtime_error(
                        "UnivariateSpline: a theoretically impossible result was found during\n"
                        "the iteration process for finding a smoothing spline with\n"
                        "fp = s. probably causes : s too small.\n"
                        "there is an approximation returned but the corresponding\n"
                        "weighted sum of squared residuals does not satisfy the\n"
                        "condition abs(fp-s)/s < tol.");
            } else if (ier == 3) {
                throw std::runtime_error(
                        "UnivariateSpline: error. the maximal number of iterations maxit (set to 20\n"
                        "by the program) allowed for finding a smoothing spline\n"
                        "with fp=s has been reached. probably causes : s too small\n"
                        "there is an approximation returned but the corresponding\n"
                        "weighted sum of squared residuals does not satisfy the\n"
                        "condition abs(fp-s)/s < tol.");
            } else {
                throw std::runtime_error(
                        "UnivariateSpline:  error. on entry, the input data are controlled on validity\n"
                        "the following restrictions must be satisfied.\n"
                        "-1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m\n"
                        "xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)\n"
                        "if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)\n"
                        "\txb<t(k+2)<t(k+3)<...<t(n-k-1)<xe\n"
                        "\tthe schoenberg-whitney conditions, i.e. there\n"
                        "\tmust be a subset of data points xx(j) such that\n"
                        "\tt(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1\n"
                        "if iopt>=0: s>=0\n"
                        "if s=0 : nest >= m+k+1\n"
                        "if one of these conditions is found to be violated,control\n"
                        "is immediately repassed to the calling program. in that\n"
                        "case there is no approximation returned.");
            }
        }
    }

    ~UnivariateSpline() {

        delete[] m_fpint;
        delete[] m_z;
        delete[] m_iwrk;
        for (int i = 0; i < m_nest; i++) {
            delete[] m_a[i];
            delete[] m_b[i];
            delete[] m_g[i];
        }
        for (int i = 0; i < m_m; i++) {
            delete[] m_q[i];
        }
    }

    /**
     * Evaluate the spline at a given x
     * @param x Location where spline is to be evaluated.
     * @return Value of spline at x.
     */
    double operator()(double x) {
        int ier = 0;
        double y = dierckx::splev(m_t.data(), m_n, m_c.data(), m_k, x, m_ext, ier);
        if (ier != 0) {
            if (ier == 0) {
                throw std::runtime_error(
                        "UnivariateSpline::operator(): argument out of bounds and ext == 2"
                );
            } else {
                throw std::runtime_error(
                        "UnivariateSpline::operator(): invalid input data"
                        " restrictions:\n"
                        "\tm >= 1\n"
                        "\tt(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1."
                );
            }
        }
        return y;
    }

    /**
     * Evaluate the spline at a given x
     * @param x Location where spline is to be evaluated.
     * @return Value of spline at x.
     */
    lanre::autodiff::Dual<double> operator()(const lanre::autodiff::Dual<double> &x) {
        return lanre::autodiff::Dual<double>{
                operator()(x.val), x.eps * derivative(x.val)
        };
    }

    /**
     * Evaluate the spline at a set of x values.
     * @param x Locations where spline is to be evaluated.
     * @return Vector of values of the spline at each x.
     */
    std::vector<double> operator()(std::vector<double> &x) {
        int ier = 0;
        std::vector<double> y(x.size());
        dierckx::splev(m_t.data(), m_n, m_c.data(), m_k, x.data(), y.data(), x.size(), m_ext, ier);
        if (ier != 0) {
            if (ier == 0) {
                throw std::runtime_error(
                        "UnivariateSpline::operator(): argument out of bounds and ext == 2"
                );
            } else {
                throw std::runtime_error(
                        "UnivariateSpline::operator(): invalid input data"
                        " restrictions:\n"
                        "\tm >= 1\n"
                        "\tt(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1."
                );
            }
        }
        return y;
    }

    /**
     * Integrate the spline from a to b.
     * @param a Lower bound of integration.
     * @param b Upper bound of integration.
     * @return Integral of spline from a to b.
     */
    double integrate(double a, double b) {
        std::vector<double> wrk(m_n, 0.0);
        return dierckx::splint(m_t.data(), m_n, m_c.data(), m_k, a, b, wrk.data());
    }

    /**
     * Compute the nu-th derivative of the spline at a x.
     * @param x Value to evaluate nu-th derivative of spline with.
     * @param nu Order of the derivative to compute.
     * @return nu-th derivative of the spline at the x.
     */
    double derivative(double x, int nu = 1) {
        if (nu > m_k || nu < 0) {
            throw std::runtime_error("UnivariateSpline::derivative : nu must be 0 <= nu <= k");
        }
        int ier;
        std::vector<double> wrk(m_n, 0.0);
        std::vector<double> xx(1, x);
        std::vector<double> y(1, 0.0);

        dierckx::splder(m_t.data(), m_n, m_c.data(), m_k, nu, xx.data(), y.data(),
                        1, m_ext, wrk.data(), ier);
        return y[0];
    }

    /**
     * Compute the nu-th derivative of the spline at a set of x values.
     * @param x Values to evaluate nu-th derivative of spline with.
     * @param nu Order of the derivative to compute.
     * @return nu-th derivative of the spline at the x values.
     */
    std::vector<double> derivative(const std::vector<double> &x, int nu = 1) {
        if (nu > m_k || nu < 0) {
            throw std::runtime_error("UnivariateSpline::derivative : nu must be 0 <= nu <= k");
        }
        int ier;
        std::vector<double> wrk(m_n, 0.0);
        std::vector<double> y(x.size(), 0.0);

        dierckx::splder(m_t.data(), m_n, m_c.data(), m_k, nu, x.data(), y.data(),
                        x.size(), m_ext, wrk.data(), ier);
        return y;
    }

    /**
     * Compute all the derivatives of the spline at 'x'. The number of
     * derivatives is equal to the order of the spline, 'k'.
     *
     * @param x Value to compute derivatives at.
     * @return Derivatives of spline at 'x'.
     */
    std::vector<double> derivatives(double x) {
        int ier;
        int k1 = m_k + 1;
        std::vector<double> d(k1, 0.0);
        dierckx::spalde(m_t.data(), m_n, m_c.data(), k1, x, d.data(), ier);
        std::cout << "ier = " << ier << std::endl;
        return d;
    }


};


}
}

#endif //LANRE_INTERPOLATE_INTERP1D_HPP
