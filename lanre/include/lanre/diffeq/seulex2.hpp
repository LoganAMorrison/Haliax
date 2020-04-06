//
// Created by Logan Morrison on 3/18/20.
//

/*
 * Numerical solution of a stiff (or differential algebraic) system of
 * first order ordinary differential equations M u'(t) = f(t, u(t)).
 * This is an extrapolation-algorithm, based on the linearly implicit
 * euler method (with step size control and order selection).
 *
 * Adapted from E. Hairer and G. Wanner
 */

#ifndef LANRE_DIFFEQ_SEULEX2_HPP
#define LANRE_DIFFEQ_SEULEX2_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/algorithm.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/solution.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "lanre/diffeq/decsol.hpp"
#include <boost/math/special_functions/pow.hpp>

namespace lanre {
namespace diffeq {

struct SeulexCache2 : public ODEAlgorithmCache {
    bool reject = false;
    bool last = false;
    bool atov = false;
    bool caljac = false;
    bool calhes = false;
    int nm1 = 0;
    int job = 0;
    // CCOSUE
    double xoldd = 0.0;
    double hhh = 0.0;
    int nrd = 0;
    int nnrd = 0;
    int kright = 0;
    int kc = 0;
    int kopt = 0;

    // LINAL
    int mle = 0;
    int mue = 0;
    int mbjac = 0;
    int mbb = 0;
    int mdiag = 0;
    int mdiff = 0;
    int mbdiag = 0;
    int lrde = 0;
    int ldjac = 0;
    int ldmas = 0;
    int lde = 0;

    int wkfcn = 0;
    int wkjac = 0;
    int wkrow = 0;
    int wkdec = 0;
    int wksol = 0;
    int ipt = 0;

    double err = 0.0;
    double errold = 0.0;
    double dtopt = 0.0;
    double dtmaxn = 0.0;
    double theta = 0.0;

    // Vectors
    Vector<double> yh{};
    Vector<double> dy{};
    Vector<double> fx{};
    Vector<double> yhh{};
    Vector<double> dyh{};
    Vector<double> del{};
    Vector<double> wh{};
    Vector<double> scal{};
    Vector<double> hh{};
    Vector<double> w{};
    Vector<double> a{};
    Vector<double> dens{};
    Vector<double> facul{};
    Vector<double> dfdt{};
    // Matrices
    Matrix<double> dfdu{};
    Matrix<double> mass_matrix{};
    Matrix<double> T{};
    Matrix<double> fsafe{};
    Matrix<double> e{};
    // Integer vectors
    Vector<int> ip{};
    Vector<int> nj{};
    Vector<int> iphes{};
    Vector<int> icomp{};
};

struct Seulex2 : public ODEAlgorithm {
    bool autnms = false; // Is problem autonomous?
    bool implct = false; // Is problem implicit?
    bool jband = false; // Is jacobian banded?

    int mujac = -1; // Upper bandwidth of jacobian
    int mljac = -1; // Lower bandwidth of jacobian
    int mumas = -1; // Upper bandwidth of mass matrix
    int mlmas = -1; // Lower bandwidth of mass matrix
    // IWORK parameters
    bool hess = false; // Transform to Hessenberg form?
    int km = 0; // maximum number of columns in extrapolation
    int nsequ = 0; // Choice of step size sequence
    int lambda = 0;
    int nrdens = 0; // Number of components for which dense output is required
    int m1 = 0;
    int m2 = 0;
    // WORK parameters
    double thet = 0.0;
    double fac1 = 0.0;
    double fac2 = 0.0;
    double fac3 = 0.0;
    double fac4 = 0.0;
    double safe1 = 0.0;
    double safe2 = 0.0;
};

static void init(ODEIntegrator &integrator, Seulex2 &alg, SeulexCache2 &cache) {
    // alg.km: Maximum number of columns in the extrapolation
    if (alg.km == 0) {
        alg.km = 12;
    } else {
        if (alg.km <= 2) {
            throw std::runtime_error("Seulex: Invalid input, km = " + std::to_string(alg.km));
        }
    }
    // alg.nsequ: Choice of step size sequence
    if (alg.nsequ == 0) {
        alg.nsequ = 2;
    }
    if (alg.nsequ <= 0 || alg.nsequ >= 5) {
        throw std::runtime_error("Seulex: Invalid input, nsequ = " + std::to_string(alg.nsequ));
    }
    // alg.lambda: Parameter for dense output
    if (alg.lambda < 0 || alg.lambda >= 2) {
        throw std::runtime_error("Seulex: Invalid input, lambda = " + std::to_string(alg.lambda));
    }
    // alg.nrdens: Number of dense output components
    if (alg.nrdens < 0 || alg.nrdens > integrator.n) {
        throw std::runtime_error("Seulex: Invalid input, nrdens = " + std::to_string(alg.nrdens));
    }
    // alg.m1, alg.m2: Parameters for second order equations
    cache.nm1 = integrator.n - alg.m1;
    if (alg.m1 == 0) {
        alg.m2 = integrator.n;
    }
    if (alg.m2 == 0) {
        alg.m2 = alg.m1;
    }
    if (alg.m1 < 0 || alg.m2 < 0 || alg.m1 + alg.m2 > integrator.n) {
        throw std::runtime_error(
                "Seulex: Invalid input m1,m2 = " +
                        std::to_string(alg.m1) +
                        ", " +
                        std::to_string(alg.m2)
        );
    }
    // alg.thet: Decides whether the jacobian should be recomputed
    if (alg.thet == 0.0) {
        alg.thet = std::min(1e-4, integrator.opts.reltol);
    }
    // alg.fac1, alg.fac2: Parameters for step size selection
    if (alg.fac1 == 0.0) {
        alg.fac1 = 0.1;
    }
    if (alg.fac2 == 0.0) {
        alg.fac2 = 4.0;
    }
    // alg.fac3, alg.fac4: Parameters for the order selection
    if (alg.fac3 == 0.0) {
        alg.fac3 = 0.7;
    }
    if (alg.fac4 == 0.0) {
        alg.fac4 = 0.9;
    }
    // alg.safe1, alg.safe2: Safety factors for step size prediction
    if (alg.safe1 == 0.0) {
        alg.safe1 = 0.6;
    }
    if (alg.safe2 == 0.0) {
        alg.safe2 = 0.93;
    }
    // cache.wkfcn, wkjac, wjdec, wksol: estimated work for func evals
    if (cache.wkfcn == 0.0) {
        cache.wkfcn = 1.0;
    }
    if (cache.wkjac == 0.0) {
        cache.wkjac = 5.0;
    }
    if (cache.wkdec == 0.0) {
        cache.wkdec = 1.0;
    }
    if (cache.wksol == 0.0) {
        cache.wksol = 1.0;
    }
    cache.wkrow = cache.wkfcn + cache.wksol;

    // Check if the tolerances are okay
    if (integrator.opts.abstol <= 0.0 || integrator.opts.reltol <= 1.0 * std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("Seulex: Tolerances are too small.");
    }


    // Computation of arrays

    if (alg.mujac < 0) {
        alg.mujac = integrator.n;
    }
    if (alg.mljac < 0) {
        alg.mljac = integrator.n;
    }
    if (alg.mumas < 0) {
        alg.mumas = integrator.n;
    }
    if (alg.mlmas < 0) {
        alg.mlmas = integrator.n;
    }

    // Computation of row-dimensions of the 2-d arrays for jac and e
    if (alg.jband) {
        cache.ldjac = alg.mljac + alg.mujac + 1;
        cache.lde = alg.mljac + cache.ldjac;
    } else {
        alg.mljac = cache.nm1;
        alg.mujac = cache.nm1;
        cache.ldjac = cache.nm1;
        cache.lde = cache.nm1;
    }
    // Mass matrix
    if (alg.implct) {
        if (alg.mlmas != cache.nm1) {
            cache.ldmas = alg.mlmas + alg.mumas + 1;
            if (alg.jband) {
                cache.job = 4;
            } else {
                cache.job = 3;
            }
        } else {
            cache.ldmas = cache.nm1;
            cache.job = 5;
        }
        // Bandwith of mass matrix not larger that bandwith of jacobian
        if (alg.mlmas > alg.mljac || alg.mumas > alg.mujac) {
            throw std::runtime_error("Seulex: Bandwidth of mass matrix not larger than bandwith of jacobian.");
        }
    } else {
        cache.ldmas = 0;
        if (alg.jband) {
            cache.job = 2;
        } else {
            cache.job = 1;
            if (integrator.n > 2 && alg.hess) {
                cache.job = 7;
            }
        }
    }
    cache.ldmas = std::max(1, cache.ldmas);

    // Hessenberg option only for explicit equ. with full jacobian
    if ((alg.implct || alg.jband) && cache.job == 7) {
        throw std::runtime_error(
                "Hessenberg option only for explicit equations with full jacobian."
        );
    }

    cache.yh.resize(integrator.n);
    cache.dy.resize(integrator.n);
    cache.dfdt.resize(integrator.n);
    cache.yhh.resize(integrator.n);
    cache.dyh.resize(integrator.n);
    cache.del.resize(integrator.n);
    cache.wh.resize(integrator.n);
    cache.scal.resize(integrator.n);

    cache.hh.resize(alg.km);
    cache.w.resize(alg.km);
    cache.a.resize(alg.km);
    cache.facul.resize(alg.km);
    cache.dens.resize((alg.km + 2) * cache.nrd);

    cache.ip.resize(cache.nm1);
    cache.iphes.resize(integrator.n);
    cache.nj.resize(alg.km);
    cache.icomp.resize(cache.nrd);

    cache.dfdu.resize(cache.ldjac, integrator.n);
    cache.mass_matrix.resize(cache.ldmas, cache.nm1);
    cache.T.resize(alg.km, cache.nm1);
    cache.fsafe.resize((alg.km + 1) * alg.km / 2, cache.nrd);

    cache.e.resize(cache.lde, cache.nm1);
}

/**
 * Perform a decomposition
 * @param integrator
 * @param alg
 * @param cache
 * @param fac
 * @return
 */
static int decomp_real(ODEIntegrator &integrator, Seulex2 &alg, SeulexCache2 &cache, double fac) {
    int mm, ier = 0;
    double sum;

    switch (cache.job) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.e(i, j) = -cache.dfdu(i, j);
                }
                cache.e(j, j) += fac;
            }
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (2):
            // mass = identity, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j);
                }
                cache.e(cache.mdiag, j) += fac;
            }
            ier = decb(integrator.n, cache.e, cache.mle, cache.mue, cache.ip);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++)
                    cache.e(i, j) = -cache.dfdu(i, j);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e(i, j) += fac * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e(i + cache.mdiff, j) += fac * cache.mass_matrix(i, j);
            }
            ier = decb(integrator.n, cache.e, cache.mle, cache.mue, cache.ip);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++)
                for (int i = 0; i < integrator.n; i++)
                    cache.e(i, j) = cache.mass_matrix(i, j) * fac - cache.dfdu(i, j);
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // This option is not provided
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.job << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            if (cache.calhes) elmhes(integrator.n, 0, integrator.n, cache.dfdu, cache.iphes);
            cache.calhes = false;
            for (int j = 0; j < integrator.n - 1; j++) cache.e(j + 1, j) = -cache.dfdu(j + 1, j);
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i <= j; i++) cache.e(i, j) = -cache.dfdu(i, j);
                cache.e(j, j) += fac;
            }
            ier = dech(integrator.n, cache.e, 1, cache.ip);
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.nm1; i++) {
                    cache.e(i, j) = -cache.dfdu(i, j + alg.m1);
                }
                cache.e(j, j) += fac;
            }
            break;

        case (12):
            // mass = identity, Jacobian a banded matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                cache.e(cache.mdiag, j) += fac;
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.nm1; i++)
                    cache.e(i, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e(i, j) += fac * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            break;

        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e(i + cache.mdiff, j) += fac * cache.mass_matrix(i, j);
            }
            break;

        case (15):
            // mass is a full matrix, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++)
                for (int i = 0; i < cache.nm1; i++)
                    cache.e(i, j) = cache.mass_matrix(i, j) * fac - cache.dfdu(i, j + alg.m1);
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.job << std::endl;
            ier = -1;
            return (ier);
    }

    switch (cache.job) {
        case (1):
        case (2):
        case (3):
        case (4):
        case (5):
        case (7):
            break;

        case (11):
        case (13):
        case (15):
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                for (int i = 0; i < cache.nm1; i++) {
                    sum = 0.0;
                    for (int k = 0; k < mm; k++)
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / fac;
                    cache.e(i, j) -= sum;
                }
            }
            ier = dec(cache.nm1, cache.e, cache.ip);
            break;

        case (12):
        case (14):
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    sum = 0.0;
                    for (int k = 0; k < mm; k++)
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / fac;
                    cache.e(i + cache.mle, j) -= sum;
                }
            }
            ier = decb(cache.nm1, cache.e, cache.mle, cache.mue, cache.ip);
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.job << std::endl;
            ier = -1;
            return ier;
    }

    return ier;

}

/**
 * Perform a linear solve
 * @param integrator
 * @param alg
 * @param cache
 * @param fac
 */
static void linear_solve(ODEIntegrator &integrator, Seulex2 &alg, SeulexCache2 &cache, double fac) {
    if (cache.job == 1) {
        sol(integrator.n, cache.e, cache.del, cache.ip);
        return;
    } else if (cache.job == 11) {
        int mm = alg.m1 / alg.m2;
        for (int j = 1; j <= alg.m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * alg.m2;
                sum = (cache.del(jkm - 1) + sum) / fac;
                for (int i = 1; i <= cache.nm1; i++) {
                    int im1 = i + alg.m1;
                    cache.del(im1 - 1) = cache.del(im1 - 1) + cache.dfdu(i - 1, jkm - 1) * sum;
                }
            }
        }
        Vector<double> del_seg = cache.del.segment(alg.m1, cache.del.size() - 1);
        sol(cache.nm1, cache.e, del_seg, cache.ip);
        cache.del.segment(alg.m1, cache.del.size() - 1) = del_seg;
        for (int i = alg.m1; i >= 1; i--) {
            cache.del(i - 1) = (cache.del(i - 1) + cache.del(alg.m2 + i - 1)) / fac;
        }
        return;
    } else if (cache.job == 2) {
        solb(integrator.n, cache.e, cache.mle, cache.mue, cache.del, cache.ip);
        return;
    } else if (cache.job == 12) {
        int mm = alg.m1 / alg.m2;
        for (int j = 1; j <= alg.m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * alg.m2;
                sum = (cache.del(jkm - 1) + sum) / fac;
                for (int i = std::max(1, j - alg.mujac); i <= std::min(cache.nm1, j + alg.mljac); i++) {
                    int im1 = i + alg.m1;
                    cache.del(im1 - 1) += cache.dfdu(i + alg.mujac + 1 - j - 1, jkm - 1) * sum;
                }
            }
        }
        Vector<double> del_seg = cache.del.segment(alg.m1, cache.del.size() - 1);
        solb(cache.nm1, cache.e, cache.mle, cache.mue, del_seg, cache.ip);
        cache.del.segment(alg.m1, cache.del.size() - 1) = del_seg;
        for (int i = alg.m1; i >= 1; i--) {
            cache.del(i - 1) = (cache.del(i - 1) + cache.del(alg.m2 + i - 1)) / fac;
        }
        return;
    } else if (cache.job == 7) {
        for (int mmm = integrator.n - 2; mmm >= 1; mmm--) {
            int mp = integrator.n - mmm;
            int mp1 = mp - 1;
            int i = cache.iphes(mp - 1);

            if (i == mp) {
                continue;
            }
            double zsafe = cache.del(i - 1);
            cache.del(mp - 1) = cache.del(i - 1);
            cache.del(i - 1) = zsafe;

            for (i = mp + 1; i <= integrator.n; i++) {
                cache.del(i - 1) -= cache.dfdu(i - 1, mp1 - 1) * cache.del(mp - 1);
            }
        }
        solh(integrator.n, cache.e, 1, cache.del, cache.ip);
        for (int mmm = 1; mmm <= integrator.n - 2; mmm++) {
            int mp = integrator.n - mmm;
            int mp1 = mp - 1;
            for (int i = mp + 1; i <= integrator.n; i++) {
                cache.del(i - 1) += cache.dfdu(i - 1, mp1 - 1) * cache.del(mp - 1);
            }
            int i = cache.iphes(mp - 1);
            if (i == mp) {
                continue;
            }
            double zsafe = cache.del(mp - 1);
            cache.del(mp - 1) = cache.del(i - 1);
            cache.del(i - 1) = zsafe;
        }
        return;
    } else {
        return;
    }
}

/**
 * Compute the j-th line of the extrapolation table and provide an estimate
 * for the optimal step size.
 * @param integrator
 * @param alg
 * @param cache
 */
static void seul(ODEIntegrator &integrator, Seulex2 &alg, SeulexCache2 &cache, int jj) {

    // Comput the matrix e and its decomposition
    double hj = integrator.dt / double(cache.nj(jj - 1));
    double hji = 1.0 / hj;
    int ier = decomp_real(integrator, alg, cache, hji);

    if (ier != 0) {
        cache.atov = true;
        integrator.dt += 0.5;
        cache.reject = true;
    }
    integrator.num_decompositions++;

    // Starting procedure
    if (!alg.autnms) {
        integrator.f->dudt(cache.dy, integrator.u, integrator.t + hj);
    }
    for (int i = 1; i <= integrator.n; i++) {
        cache.yh(i - 1) = integrator.u(i - 1);
        cache.del(i - 1) = cache.dy(i - 1);
    }

    linear_solve(integrator, alg, cache, hji);
    integrator.num_linear_solves++;
    int m = cache.nj(jj - 1);
    if (integrator.opts.dense && m == jj) {
        cache.ipt += 1;
        for (int i = 1; i <= cache.nrd; i++) {
            cache.fsafe(cache.ipt - 1, i - 1) = cache.del(cache.icomp(i - 1));
        }
    }

    // Semi-implicit euler method
    if (m > 1) {
        for (int mm = 1; mm <= m - 1; mm++) {
            for (int i = 1; i <= integrator.n; i++) {
                cache.yh(i - 1) += cache.del(i - 1);
            }
            if (alg.autnms) {
                integrator.f->dudt(cache.dyh, cache.yh, integrator.t + hj);
            } else {
                integrator.f->dudt(cache.dyh, cache.yh, integrator.t + hj * (mm + 1));
            }
            integrator.num_function_evaluations++;

            if (mm == 1 && jj <= 2) {
                // Stability check
                double del1 = 0.0;
                for (int i = 1; i <= integrator.n; i++) {
                    del1 += boost::math::pow<2>(cache.del(i - 1) / cache.scal(i - 1));
                }
                del1 = std::sqrt(del1);
                if (alg.implct) {
                    for (int i = 1; i <= cache.nm1; i++) {
                        cache.wh(i - 1) = cache.del(i + alg.m1 - 1);
                    }
                    if (alg.mlmas == cache.nm1) {
                        for (int i = 1; i <= cache.nm1; i++) {
                            double sum = 0.0;
                            for (int j = 1; j <= cache.nm1; j++) {
                                sum += cache.mass_matrix(i - 1, j - 1) * cache.wh(j - 1);
                            }
                            cache.del(i + alg.m1 - 1) = sum;
                        }
                    } else {
                        for (int i = 1; i <= cache.nm1; i++) {
                            double sum = 0.0;
                            for (int j = std::max(1, i - alg.mlmas); j <= std::min(cache.nm1, i + alg.mumas); j++) {
                                sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j - 1) * cache.wh(j - 1);
                            }
                            cache.del(i + alg.m1 - 1) = sum;
                        }
                    }
                }

                if (!alg.autnms) {
                    integrator.f->dudt(cache.wh, cache.yh, integrator.t + hj);
                    integrator.num_function_evaluations++;
                    for (int i = 1; i <= integrator.n; i++) {
                        cache.del(i - 1) = cache.wh(i - 1) - cache.del(i - 1) * hji;
                    }
                } else {
                    for (int i = 1; i <= integrator.n; i++) {
                        cache.del(i - 1) = cache.dyh(i - 1) - cache.del(i - 1) * hji;
                    }
                }

                linear_solve(integrator, alg, cache, hji);
                integrator.num_linear_solves++;

                double del2 = 0.0;
                for (int i = 1; i <= integrator.n; i++) {
                    del2 += boost::math::pow<2>(cache.del(i - 1) / cache.scal(i - 1));
                }
                del2 = std::sqrt(del2);
                cache.theta = del2 / std::max(1.0, del1);
                if (cache.theta > 1.0) {
                    cache.atov = true;
                    integrator.dt += 0.5;
                    cache.reject = true;
                }
            }

            linear_solve(integrator, alg, cache, hji);
            integrator.num_linear_solves++;

            for (int i = 1; i <= integrator.n; i++) {
                cache.del(i - 1) = cache.dyh(i - 1);
            }

            if (integrator.opts.dense) {
                cache.ipt += 1;
                for (int i = 1; i <= cache.nrd; i++) {
                    cache.fsafe(cache.ipt - 1, i - 1) = cache.del(cache.icomp(i - 1));
                }
            }
        }
    }

    for (int i = 1; i <= integrator.n; i++) {
        cache.T(jj - 1, i - 1) = cache.yh(i - 1) + cache.del(i - 1);
    }

    // Polynomial extrapolation
    if (jj == 1) {
        return;
    }
    for (int l = jj; l >= 2; l--) {
        double fac = (double(cache.nj(jj - 1)) / double(cache.nj(l - 1 - 1))) - 1.0;
        for (int i = 1; i <= integrator.n; i++) {
            cache.T(l - 1 - 1, i - 1) = (
                    cache.T(l - 1, i - 1) +
                            (cache.T(l - 1, i - 1) - cache.T(l - 1 - 1, i - 1)) / fac);
        }
    }
    cache.err = 0.0;
    for (int i = 1; i <= integrator.n; i++) {
        cache.err += boost::math::pow<2>(
                std::min(std::abs((cache.T(1 - 1, i - 1) - cache.T(2 - 1, i - 1))) / cache.scal(i - 1), 1.15));
    }
    if (cache.err >= 1.30) {
        cache.atov = true;
        integrator.dt += 0.5;
        cache.reject = true;
    }
    cache.err = std::sqrt(cache.err / double(integrator.n));
    if (jj > 2 && cache.err >= cache.errold) {
        cache.atov = true;
        integrator.dt += 0.5;
        cache.reject = true;
    }
    cache.errold = std::max(4.0 * cache.err, 1.0);

    // Compute the optimal step sizes
    double expo = 1.0 / double(jj);
    double facmin = std::pow(alg.fac1, expo);
    double fac = std::min(alg.fac2 / facmin, std::max(facmin, std::pow(cache.err / alg.safe1, expo) / alg.safe2));
    fac = 1.0 / fac;
    cache.hh(jj - 1) = std::min(std::abs(integrator.dt) * fac, cache.dtmaxn);
    cache.w(jj - 1) = cache.a(jj - 1) / cache.hh(jj - 1);
    return;
}

/**
 * Core integrator for Seulex
 * @param integrator
 * @param alg
 * @param cache
 * @return
 */
static ODESolution seulex_core(ODEIntegrator &integrator, Seulex2 &alg, SeulexCache2 &cache) {
    ODESolution solution{};

    solution.ts.push_back(integrator.t);
    solution.us.push_back(integrator.u);

    init(integrator, alg, cache);
    // Compute the coefficients for dense output
    if (integrator.opts.dense) {
        cache.nnrd = cache.nrd;
        cache.facul(0) = 1.0;
        for (int i = 1; i <= alg.km - 1; i++) {
            cache.facul(i) = i * cache.facul(i - 1);
        }
    }

    // Compute the mass matrix for the implicit case
    if (alg.implct) {
        cache.mass_matrix = integrator.f->mass_matrix;
    }

    // Initializations
    cache.lrde = (alg.km + 2) * cache.nrd;
    cache.mle = alg.mljac;
    cache.mue = alg.mujac;
    cache.mbjac = alg.mljac + alg.mujac + 1;
    cache.mbb = alg.mumas + alg.mlmas + 1;
    cache.mdiag = cache.mle + cache.mue + 1;
    cache.mdiff = cache.mle + cache.mue - alg.mumas;
    cache.mbdiag = alg.mumas + 1;
    if (alg.m1 > 0) {
        cache.job += 10;
    }

    // Define the step size sequence
    if (alg.nsequ == 1) {
        cache.nj(1 - 1) = 1;
        cache.nj(2 - 1) = 2;
        cache.nj(3 - 1) = 3;
        for (int i = 4; i <= alg.km; i++) {
            cache.nj(i - 1) = 2 * cache.nj(i - 2 - 1);
        }
    }
    if (alg.nsequ == 2) {
        cache.nj(1 - 1) = 2;
        cache.nj(2 - 1) = 3;
        for (int i = 3; i <= alg.km; i++) {
            cache.nj(i - 1) = 2 * cache.nj(i - 2 - 1);
        }
    }
    for (int i = 1; i <= alg.km; i++) {
        if (alg.nsequ == 3) {
            cache.nj(i - 1) = i;
        } else if (alg.nsequ == 4) {
            cache.nj(i - 1) = i + 1;
        }
    }
    cache.a(0) = cache.wkjac + cache.nj(0) * cache.wkrow + cache.wkdec;
    for (int i = 2; i <= alg.km; i++) {
        cache.a(i - 1) = cache.a(i - 1 - 1) + (cache.nj(i - 1) - 1) * cache.wkrow + cache.wkdec;
    }

    int k = std::max(2,
                     std::min(alg.km - 2, int(-log10(integrator.opts.reltol + integrator.opts.abstol) * 0.6 + 1.5)));
    cache.dtmaxn = std::min(std::abs(integrator.opts.dtmax), std::abs(integrator.t_final - integrator.t));
    integrator.dt = std::max(std::abs(integrator.dt), 1e-6);
    cache.theta = 2 * std::abs(alg.thet);

    if (integrator.opts.dense) {
        // do something...
    }

    cache.err = 0.0;
    cache.w(0) = 1.3;
    for (int i = 1; i <= integrator.n; i++) {
        cache.scal(i - 1) = integrator.opts.abstol + integrator.opts.reltol * std::abs(integrator.u(i - 1));
    }

    cache.caljac = false;
    cache.reject = false;
    cache.last = false;

    goto10:
    if (cache.reject) {
        cache.theta = 2.0 * std::abs(alg.thet);
    }
    cache.atov = false;
    // Is tfinal reached in the next step?
    if (0.1 * std::abs(integrator.t_final - integrator.t) <= std::abs(integrator.t) *
            std::numeric_limits<double>::epsilon()) {
        goto goto110;
    }
    cache.dtopt = integrator.dt;
    integrator.dt = integrator.tdir * std::min(
            std::abs(integrator.dt),
            std::min(
                    std::abs(integrator.t_final - integrator.t),
                    cache.dtmaxn)
    );
    if ((integrator.t + 1.01 * integrator.dt - integrator.t_final) * integrator.tdir > 0.0) {
        integrator.dt = integrator.t_final - integrator.t;
        cache.last = true;
    }
    if (alg.autnms) {
        integrator.f->dudt(cache.dy, integrator.u, integrator.t);
        integrator.num_function_evaluations++;
    }
    if (cache.theta > alg.thet && !cache.caljac) {
        integrator.num_jacobian_evaluations++;
        integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
        cache.caljac = true;
        cache.calhes = true;
    }

    // The first and last step
    if (integrator.num_steps == 0 || cache.last) {
        cache.ipt = 0;
        integrator.num_steps++;
        for (int j = 1; j <= k; j++) {
            cache.kc = j;
            seul(integrator, alg, cache, j);
            if (cache.atov) {
                goto goto10;
            }
            if (j > 1 && cache.err <= 1.0) {
                goto goto60;
            }
        }
        goto goto55;
    }

    // Basic integration step
    goto30:
    cache.ipt = 0;
    integrator.num_steps++;
    if (integrator.num_steps >= integrator.opts.max_num_steps) {
        goto goto120;
    }
    cache.kc = k - 1;
    for (int j = 1; j <= cache.kc; j++) {
        seul(integrator, alg, cache, j);
        if (cache.atov) {
            goto goto10;
        }
    }

    // Convergence monitor
    if (k == 2 || cache.reject) {
        goto goto50;
    }
    if (cache.err <= 1.0) {
        goto goto60;
    }
    if (cache.err > double(cache.nj(k + 1 - 1) * cache.nj(k - 1)) * 4.0) {
        goto goto100;
    }

    goto50:
    seul(integrator, alg, cache, k);
    if (cache.atov) {
        goto goto10;
    }
    cache.kc = k;
    if (cache.err <= 1.0) {
        goto goto60;
    }

    // Hope for convergence in line k + 1
    goto55:
    if (cache.err > double(cache.nj(k + 1 - 1)) * 2.0) {
        goto goto100;
    }
    cache.kc = k + 1;
    seul(integrator, alg, cache, cache.kc);
    if (cache.atov) {
        goto goto10;
    }
    if (cache.err > 1.0) {
        goto goto100;
    }

    // Step is accepted
    goto60:
    integrator.t += integrator.dt;
    if (integrator.opts.dense) {
        cache.kright = cache.kc;
        for (int i = 1; i <= cache.nrd; i++) {
            cache.dens(i - 1) = integrator.u(cache.icomp(i - 1));
        }
    }

    for (int i = 1; i <= integrator.n; i++) {
        double t1i = cache.T(1 - 1, i - 1);
        cache.scal(i - 1) = integrator.opts.abstol + integrator.opts.reltol * t1i;
        integrator.u(i - 1) = t1i;
    }
    integrator.num_accept++;
    cache.caljac = false;

    if (integrator.opts.dense) {
        cache.xoldd = integrator.tprev;
        cache.hhh = integrator.dt;
        for (int i = 1; i <= cache.nrd; i++) {
            cache.dens(cache.nrd + i - 1) = integrator.u(cache.icomp(i - 1));
        }
        for (int klr = 1; klr <= cache.kright - 1; klr++) {
            // Compute the differences
            if (klr >= 2) {
                for (int kk = klr; kk <= cache.kc; kk++) {
                    int lbeg = ((kk + 1) * kk) / 2;
                    int lend = lbeg - kk + 2;
                    for (int l = lbeg; l >= lend; l--) {
                        for (int i = 1; i <= cache.nrd; i++) {
                            cache.fsafe(l - 1, i - 1) -= cache.fsafe(l - 1 - 1, i - 1);
                        }
                    }
                }
            }
            // Compute the derivatives at the right end
            for (int kk = klr + alg.lambda; kk <= cache.kc; kk++) {
                auto facnj = double(cache.nj(kk - 1));
                facnj = std::pow(facnj, klr) / cache.facul(klr + 1 - 1);
                cache.ipt = ((kk + 1) * kk) / 2;
                for (int i = 1; i <= cache.nrd; i++) {
                    int krn = (kk - alg.lambda + 1) * cache.nrd;
                    cache.dens(krn + i - 1) = cache.fsafe(cache.ipt - 1, i - 1) * facnj;
                }
            }
            for (int j = klr + alg.lambda + 1; j <= cache.kc; j++) {
                auto dblenj = double(cache.nj(j - 1));
                for (int l = j; j >= klr + alg.lambda + 1; l--) {
                    double factor = dblenj / double(cache.nj(l - 1 - 1)) - 1.0;
                    for (int i = 1; i <= cache.nrd; i++) {
                        int krn = (l - alg.lambda + 1) * cache.nrd + i;
                        cache.dens(krn - cache.nrd - 1) = (
                                cache.dens(krn - 1) +
                                        (cache.dens(krn - 1) - cache.dens(krn - cache.nrd - 1)) / factor);
                    }
                }
            }
        }
        // Copute the coefficients of the interpolation polynomial
        for (int in = 1; in <= cache.nrd; in++) {
            for (int j = 1; j <= cache.kright; j++) {
                int ii = cache.nrd * j + in;
                cache.dens(ii - 1) -= cache.dens(ii - cache.nrd - 1);
            }
        }
    }

    if (integrator.opts.dense) {
        // Do something...
    }
    solution.ts.push_back(integrator.t);
    solution.us.push_back(integrator.u);

    // Compute the optimal order
    if (cache.kc == 2) {
        cache.kopt = std::min(3, alg.km - 1);
        if (cache.reject) {
            cache.kopt = 2;
        }
        goto goto80;
    }
    if (cache.kc <= k) {
        cache.kopt = cache.kc;
        if (cache.w(cache.kc - 1 - 1) < cache.w(cache.kc - 1) * alg.fac3) {
            cache.kopt = cache.kc - 1;
        }
        if (cache.w(cache.kc - 1) < cache.w(cache.kc - 2) * alg.fac4) {
            cache.kopt = std::min(cache.kc + 1, alg.km - 1);
        }
    } else {
        cache.kopt = cache.kc - 1;
        if (cache.kc > 3 && cache.w(cache.kc - 3) < cache.w(cache.kc - 2) * alg.fac3) {
            cache.kopt = cache.kc - 2;
        }
        if (cache.w(cache.kc - 1) < cache.w(cache.kopt - 1) * alg.fac4) {
            cache.kopt = std::min(cache.kc, alg.km - 1);
        }
    }

    // After a rejected step
    goto80:
    if (cache.reject) {
        k = std::min(cache.kopt, cache.kc);
        integrator.dt = integrator.tdir * std::min(std::abs(integrator.dt), std::abs(cache.hh(k - 1)));
        cache.reject = false;
        goto goto10;
    }
    if (cache.kopt <= cache.kc) {
        integrator.dt = cache.hh(cache.kopt - 1);
    } else {
        if (cache.kc < k && cache.w(cache.kc - 1) < cache.w(cache.kc - 2) * alg.fac4) {
            integrator.dt = cache.hh(cache.kc - 1) * cache.a(cache.kopt) / cache.a(cache.kc - 1);
        } else {
            integrator.dt = cache.hh(cache.kc - 1) * cache.a(cache.kopt - 1) / cache.a(cache.kc - 1);
        }
    }
    k = cache.kopt;
    integrator.dt *= integrator.tdir;
    goto goto10;

    // Step is rejected
    goto100:
    k = std::min(k, std::min(cache.kc, alg.km - 1));
    if (k > 2 && cache.w(k - 2) < cache.w(k - 1) * alg.fac3) {
        k -= 1;
    }
    integrator.num_reject++;
    integrator.dt = integrator.tdir * cache.hh(k - 1);
    cache.last = false;
    cache.reject = true;
    if (cache.caljac) {
        goto goto30;
    }
    goto goto10;

    // Solution exit
    goto110:
    integrator.dt = cache.dtopt;
    return solution;

    // Fail exit
    goto120:
    std::cout << "Seulex: exit at t = " << integrator.t
              << ", dt = " << integrator.dt << std::endl;
    return solution;
}


ODESolution solve(ODEProblem &prob, Seulex2 &alg, ODEIntegratorOptions &opts) {
    ODEIntegrator integrator{prob, opts};
    SeulexCache2 cache{};
    return seulex_core(integrator, alg, cache);
}

ODESolution solve(ODEProblem &prob, Seulex2 &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}

}
}

#endif //LANRE_DIFFEQ_SEULEX2_HPP
