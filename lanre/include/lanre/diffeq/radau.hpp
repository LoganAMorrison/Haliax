//
// Created by Logan Morrison on 3/16/20.
//

#ifndef LANRE_DIFFEQ_RADAU_HPP
#define LANRE_DIFFEQ_RADAU_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/algorithm.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/solution.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "lanre/diffeq/decsol.hpp"
#include <exception>
#include <iostream>

namespace lanre {
namespace diffeq {


struct RadauToleranceTooSmall : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Tolerances are too small.";
    }
};

struct RadauLinearAlgebraError : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Error encountered in linear algebra routine.";
    }
};

struct RadauMaxStepExceeded : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Maximum number of steps exceeded.";
    }
};

struct RadauStepSizeUnderflow : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Step size too small.";
    }
};

struct RadauSingularMatrixEncountered : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Encountered repeatedly singular matrix.";
    }
};

struct RadauInvalidJobOption : public std::exception {
    [[nodiscard]] const char *what() const noexcept override {
        return "Radau: Job option not implemented.";
    }
};

/*
 * Cache parameters:
 *      Flags:
 *          cache.implct, cache.jband
 *          cache.caljac, cache.calhes
 *          cache.first, cache.loop
 *      Scalars:
 *          ints:
 *              cache.mle, cache.mue
 *              cache.mdiff, cache.mbjac, cache.mbb, cache.mdiag, cache.mbdiag
 *              cache.ldmas, cache.ldjac, cache.lde1
 *              cache.newt, cache.ier
 *          floats:
 *              cache.err
 *              m_fac1, cache.theta
 *              m_dtfac, m_dtopt
 *              cache.faccon, cache.thqold, cache.cfac
 *              m_dtacc, m_erracc
 *              cache.alphn, cache.betan
 *      Vectors:
 *          cache.u0
 *          cache.z1, cache.z2, cache.z3
 *          cache.z1_seg, cache.z2_seg, cache.z3_seg
 *          cache.f1, cache.f2, cache.f3
 *          cache.cont, cache.cont_seg
 *          cache.scal
 *          cache.ip1, cache.ip2, cache.iphes
 *      Matrices:
 *          cache.dfdu, cache.e1, cache.e2r, cache.e2i
 */

/*
 * Input parameters:
 *  ijac, alg.mljac, alg.mujac
 *  
 *  imas, alg.mumas, alg.mlmas
 *  
 *  alg.hess
 *  max_allowed_steps
 *  alg.nit
 *  alg.startn
 *  alg.nind1, alg.nind2, alg.nind3
 *  alg.pred
 *  alg.m1, alg.m2
 *  alg.thet
 *  alg.fnewt
 *  alg.quot1, alg.quot2
 *  m_dtmax
 *  alg.facl, alg.facr
 *  alg.alph
 *  alg.beta
 *  
 */

struct Radau5 {
    const double t11 = 9.1232394870892942792e-02;
    const double t12 = -0.14125529502095420843;
    const double t13 = -3.0029194105147424492e-02;
    const double t21 = 0.24171793270710701896;
    const double t22 = 0.20412935229379993199;
    const double t23 = 0.38294211275726193779;
    const double t31 = 0.96604818261509293619;
    const double ti11 = 4.3255798900631553510;
    const double ti12 = 0.33919925181580986954;
    const double ti13 = 0.54177053993587487119;
    const double ti21 = -4.1787185915519047273;
    const double ti22 = -0.32768282076106238708;
    const double ti23 = 0.47662355450055045196;
    const double ti31 = -0.50287263494578687595;
    const double ti32 = 2.5719269498556054292;
    const double ti33 = -0.59603920482822492497;

    // alg.c1 = (4 - sqrt(6)) / 10
    const double c1 = 0.15505102572168219018;
    // alg.c2 = (4 + sqrt(6)) / 10
    const double c2 = 0.64494897427831780982;
    // alg.c1m1 = alg.c1 - 1 = -(6 + sqrt(6)) / 10
    const double c1m1 = -0.84494897427831780982;
    // alg.c2m1 = alg.c2 - 1 = (sqrt(6) - 6) / 10
    const double c2m1 = -0.35505102572168219018;
    // alg.c1mc2 = alg.c1 - alg.c2 = - sqrt(6) / 5
    const double c1mc2 = -0.48989794855663561964;
    // u1 = 3 - 3^(1/3) + 3^(2/3)
    const double u1 = 3.6378342527444957322;
    // alph = ( 6 + 3^(1/3) - 3^(2/3) ) / 2
    const double alph = 2.6810828736277521339;
    // beta = sqrt[3(6 + 3^(4/3) + 3^(2/3) )] / 2
    const double beta = 3.0504301992474105694;
    // cno = (6 - 3^(1/3) + 3^(2/3) ) / 60
    const double cno = 0.060630570879074928870;

    double safe = 0.9;
    int nit = 0, ijob = 0, m1 = 0, m2 = 0, nm1 = 0, npred = 0, imas = 0;
    double thet = 0;
    double quot1 = 0, quot2 = 0;
    int nind1 = 0, nind2 = 0, nind3 = 0;
    bool startn = false, pred = false, hess = false;
    double facr = 0, facl = 0;
    double fnewt = 0;
    int mujac = -1, mljac = -1;
    int mumas = -1, mlmas = -1;

    Radau5() = default;
};

struct Radau5Cache {
    int ijob;
    int mle, mue;
    int mdiff, mbjac, mbb, mdiag, mbdiag;
    bool implct, jband;
    int ldmas, hess, ldjac, lde1;

    int ier = 0, newt;

    double alphn = 0, betan = 0;
    double err = 0, fac1 = 0, theta;
    double dtfac = 0, faccon, thqold, cfac, dtacc, erracc, dtopt;
    bool caljac = true, calhes = true, first = true, loop = false;
    bool reject = false, last = false;

    int num_singular;

    Vector<double> u0;
    Vector<double> z1, z2, z3;
    Vector<double> z1_seg, z2_seg, z3_seg;
    Vector<double> f1, f2, f3;
    Vector<double> cont, cont_seg;
    Vector<double> scal;
    Matrix<double> dfdu, e1, e2r, e2i;
    Matrix<double> mass_matrix;
    Vector<int> ip1, ip2, iphes;

    int numSingular = 0;

    Radau5Cache() = default;
};


static void init(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    // Check and change the tolerances
    if ((integrator.opts.abstol <= 0.0) ||
            (integrator.opts.reltol <= 10.0 * std::numeric_limits<double>::epsilon())) {
        std::cout << " tolerances are too small" << std::endl;
        throw -1;
    } else {
        double quot = integrator.opts.abstol / integrator.opts.reltol;
        integrator.opts.reltol = 0.1 * pow(integrator.opts.reltol, 2.0 / 3.0);
        integrator.opts.abstol = integrator.opts.reltol * quot;
    }

    // initial step length
    if (fabs(integrator.dt) < 10.0 * std::numeric_limits<double>::epsilon()) {
        integrator.dt = 1.0e-6;
    }

    // max step length

    // facl, facr--parameters for step size selection
    if (alg.facl == 0.0) {
        alg.facl = 5.0;
    }
    if (alg.facr == 0.0) {
        alg.facr = 1.0 / 8.0;
    }
    if ((alg.facl < 1.0) || (alg.facr > 1.0)) {
        std::cout << " curious input facl, facr = " << alg.facl << ", " << alg.facr << std::endl;
        throw -1;
    }

    // nit--maximal number of Newton iterations
    if (alg.nit == 0) {
        alg.nit = 7;
    }
    if (alg.nit <= 0) {
        std::cout << " curious input nit = " << alg.nit << std::endl;
        throw -1;
    }

    // startn--switch for starting values of Newton iterations
    alg.startn = alg.startn != 0;

    // parameters for differential-algebraic components
    if (alg.nind1 == 0) {
        alg.nind1 = integrator.n;
    }
    if (alg.nind1 + alg.nind2 + alg.nind3 != integrator.n) {
        std::cout << " curious input for nind1, nind2, nind3 = "
                  << alg.nind1 << alg.nind2 << alg.nind3 << std::endl;
        throw -1;
    }

    // pred--step size control
    alg.pred = alg.npred <= 1;

    // parameters for second order equations
    alg.nm1 = integrator.n - alg.m1;
    if (alg.m1 == 0) {
        alg.m2 = integrator.n;
    }
    if (alg.m2 == 0) {
        alg.m2 = alg.m1;
    }
    if ((alg.m1 < 0) || (alg.m2 < 0) || (alg.m1 + alg.m2 > integrator.n)) {
        std::cout << " curious input for m_m1, m_m2 = " << alg.m1 << ", " << alg.m2 << std::endl;
        throw -1;
    }

    // m_fnewt--stopping criterion for Newton's method, usually chosen < 1
    if (alg.fnewt == 0.0) {
        alg.fnewt = std::max(10.0 * std::numeric_limits<double>::epsilon() / integrator.opts.reltol,
                             std::min(0.03, sqrt(integrator.opts.reltol)));
    }
    if (alg.fnewt <= std::numeric_limits<double>::epsilon() / integrator.opts.reltol) {
        std::cout << " curious input for m_fnewt = " << alg.fnewt << std::endl;
        throw -1;
    }

    // quot1 and quot2--if quot1 < m_dtNew/hold < quot2, step size = const
    if (alg.quot1 == 0.0) {
        alg.quot1 = 1.0;
    }
    if (alg.quot2 == 0.0) {
        alg.quot2 = 1.2;
    }
    if ((alg.quot1 > 1.0) || (alg.quot2 < 1.0)) {
        std::cout << " curious input for quot1, quot2 = " << alg.quot1 << ", " << alg.quot2 << std::endl;
        throw -1;
    }

    // thet--decides whether the Jacobian should be recomputed
    if (alg.thet == 0.0) {
        alg.thet = 0.001;
    }
    if (alg.thet >= 1.0) {
        std::cout << " curious input for thet = " << alg.thet << std::endl;
        throw -1;
    }

    if (alg.mljac < 0) {
        alg.mljac = integrator.n;
    }
    if (alg.mujac < 0) {
        alg.mujac = integrator.n;
    }
    if (alg.mlmas < 0) {
        alg.mlmas = integrator.n;
    }
    if (alg.mumas < 0) {
        alg.mumas = integrator.n;
    }

    // implicit, banded or not?
    cache.implct = (alg.imas != 0);
    cache.jband = (alg.mljac < alg.nm1);

    // Computation of the row-dimensions of the 2-D arrays
    // Jacobian and matrices e1, e2
    if (cache.jband) {
        cache.ldjac = alg.mljac + alg.mujac + 1;
        cache.lde1 = alg.mljac + cache.ldjac;
    } else {
        alg.mljac = alg.nm1;
        alg.mujac = alg.nm1;
        cache.ldjac = alg.nm1;
        cache.lde1 = alg.nm1;
    }


    // mass matrix
    if (cache.implct) {
        if (alg.mlmas != alg.nm1) {
            cache.ldmas = alg.mlmas + alg.mumas + 1;
            if (cache.jband) {
                cache.ijob = 4;
            } else {
                cache.ijob = 3;
            }
        } else {
            alg.mumas = alg.nm1;
            cache.ldmas = alg.nm1;
            cache.ijob = 5;
        }
        // bandwith of "mas" not smaller than bandwith of "jac"
        if ((alg.mlmas > alg.mljac) || (alg.mumas > alg.mujac)) {
            std::cout << "bandwith of 'mas' not smaller than bandwith of 'jac'" << std::endl;
            throw -1;
        }
    } else {
        cache.ldmas = 0;
        if (cache.jband) {
            cache.ijob = 2;
        } else {
            cache.ijob = 1;
            if ((integrator.n > 2) && (alg.hess != 0)) {
                cache.ijob = 7;
            }
        }
    }
    cache.ldmas = std::max(1, cache.ldmas);

    // Hessenberg option only for explicit equations with full Jacobian
    if ((cache.implct || cache.jband) && (cache.ijob == 7)) {
        std::cout << " Hessenberg option only for explicit equations with " <<
                  "full Jacobian" << std::endl;
        throw -1;
    }

    // for second-order equations increase ijob by 10
    if (alg.m1 > 0) {
        cache.ijob += 10;
    }

    // Define constants used in linear algebra routines
    cache.mle = alg.mljac;
    cache.mue = alg.mujac;
    cache.mbjac = alg.mljac + alg.mujac + 1;
    cache.mbb = alg.mlmas + alg.mumas + 1;
    cache.mdiag = cache.mle + cache.mue;
    cache.mdiff = cache.mle + cache.mue - alg.mumas;
    cache.mbdiag = alg.mumas + 1;

    cache.z1.resize(integrator.n);
    cache.z2.resize(integrator.n);
    cache.z3.resize(integrator.n);
    cache.u0.resize(integrator.n);
    cache.scal.resize(integrator.n);
    cache.f1.resize(integrator.n);
    cache.f2.resize(integrator.n);
    cache.f3.resize(integrator.n);
    cache.cont.resize(4 * integrator.n);
    cache.ip1.resize(alg.nm1);
    cache.ip2.resize(alg.nm1);
    cache.iphes.resize(alg.nm1);

    cache.dfdu.resize(cache.ldjac, integrator.n);
    cache.e1.resize(cache.lde1, alg.nm1);
    cache.e2r.resize(cache.lde1, alg.nm1);
    cache.e2i.resize(cache.lde1, alg.nm1);

    cache.cfac = alg.safe * (1 + 2 * alg.nit);

    if (cache.implct) {
        cache.mass_matrix.resize(cache.ldmas, integrator.n);
        for (int i = 0; i < cache.ldmas; i++) {
            for (int j = 0; j < integrator.n; j++) {
                cache.mass_matrix(i, j) = integrator.f->mass_matrix(i, j);
            }
        }
    }

}


static int decomp_real(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    int mm, ier = 0;
    double sum;

    switch (cache.ijob) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.e1(i, j) = -cache.dfdu(i, j);
                }
                cache.e1(j, j) += cache.fac1;
            }
            ier = dec(integrator.n, cache.e1, cache.ip1);
            break;

        case (2):
            // mass = identity, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e1(i + cache.mle, j) = -cache.dfdu(i, j);
                }
                cache.e1(cache.mdiag, j) += cache.fac1;
            }
            ier = decb(integrator.n, cache.e1, cache.mle, cache.mue, cache.ip1);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++)
                    cache.e1(i, j) = -cache.dfdu(i, j);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e1(i, j) += cache.fac1 * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            ier = dec(integrator.n, cache.e1, cache.ip1);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e1(i + cache.mle, j) = -cache.dfdu(i, j);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e1(i + cache.mdiff, j) += cache.fac1 * cache.mass_matrix(i, j);
            }
            ier = decb(integrator.n, cache.e1, cache.mle, cache.mue, cache.ip1);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++)
                for (int i = 0; i < integrator.n; i++)
                    cache.e1(i, j) = cache.mass_matrix(i, j) * cache.fac1 - cache.dfdu(i, j);
            ier = dec(integrator.n, cache.e1, cache.ip1);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // This option is not provided
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            if (cache.calhes) elmhes(integrator.n, 0, integrator.n, cache.dfdu, cache.iphes);
            cache.calhes = false;
            for (int j = 0; j < integrator.n - 1; j++) cache.e1(j + 1, j) = -cache.dfdu(j + 1, j);
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i <= j; i++) cache.e1(i, j) = -cache.dfdu(i, j);
                cache.e1(j, j) += cache.fac1;
            }
            ier = dech(integrator.n, cache.e1, 1, cache.ip1);
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < alg.nm1; i++) {
                    cache.e1(i, j) = -cache.dfdu(i, j + alg.m1);
                }
                cache.e1(j, j) += cache.fac1;
            }
            break;

        case (12):
            // mass = identity, Jacobian a banded matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e1(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                cache.e1(cache.mdiag, j) += cache.fac1;
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < alg.nm1; i++)
                    cache.e1(i, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e1(i, j) += cache.fac1 * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            break;

        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e1(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e1(i + cache.mdiff, j) += cache.fac1 * cache.mass_matrix(i, j);
            }
            break;

        case (15):
            // mass is a full matrix, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++)
                for (int i = 0; i < alg.nm1; i++)
                    cache.e1(i, j) = cache.mass_matrix(i, j) * cache.fac1 - cache.dfdu(i, j + alg.m1);
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return (ier);
    }

    switch (cache.ijob) {
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
                for (int i = 0; i < alg.nm1; i++) {
                    sum = 0.0;
                    for (int k = 0; k < mm; k++)
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / cache.fac1;
                    cache.e1(i, j) -= sum;
                }
            }
            ier = dec(alg.nm1, cache.e1, cache.ip1);
            break;

        case (12):
        case (14):
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    sum = 0.0;
                    for (int k = 0; k < mm; k++)
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / cache.fac1;
                    cache.e1(i + cache.mle, j) -= sum;
                }
            }
            ier = decb(alg.nm1, cache.e1, cache.mle, cache.mue, cache.ip1);
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;
    }

    return ier;

} // decomp_real


static int decomp_complex(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {

    int mm, ier = 0;
    double bb, ffma, abno, alp, bet, sumr, sumi, sums;

    switch (cache.ijob) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.e2r(i, j) = -cache.dfdu(i, j); // cache.e2r\[(.*?)\]\[(.*?)\] -> cache.e2r($1,$2)
                    cache.e2i(i, j) = 0.0;
                }
                cache.e2r(j, j) += cache.alphn;
                cache.e2i(j, j) = cache.betan;
            }
            ier = decc(integrator.n, cache.e2r, cache.e2i, cache.ip2);
            break;

        case (2):
            // mass = identiy, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e2r(i + cache.mle, j) = -cache.dfdu(i, j);
                    cache.e2i(i + cache.mle, j) = 0.0;
                }
                cache.e2r(cache.mdiag, j) += cache.alphn;
                cache.e2i(cache.mdiag, j) = cache.betan;
            }
            ier = decbc(integrator.n, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.ip2);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.e2r(i, j) = -cache.dfdu(i, j);
                    cache.e2i(i, j) = 0.0;
                }
            }
            for (int j = 0; j < integrator.n; j++) {
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++) {
                    bb = cache.mass_matrix(i - j + cache.mbdiag - 1, j);
                    cache.e2r(i, j) += cache.alphn * bb;
                    cache.e2i(i, j) = cache.betan * bb;
                }
            }
            ier = decc(integrator.n, cache.e2r, cache.e2i, cache.ip2);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e2r(i + cache.mle, j) = -cache.dfdu(i, j);
                    cache.e2i(i + cache.mle, j) = 0.0;
                }
                for (int i = std::max(0, alg.mumas - j); i < std::min(cache.mbb, alg.mumas - j + integrator.n); i++) {
                    bb = cache.mass_matrix(i, j);
                    cache.e2r(i + cache.mdiff, j) += cache.alphn * bb;
                    cache.e2i(i + cache.mdiff, j) = cache.betan * bb;
                }
            }
            ier = decbc(integrator.n, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.ip2);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    bb = cache.mass_matrix(i, j);
                    cache.e2r(i, j) = cache.alphn * bb - cache.dfdu(i, j);
                    cache.e2i(i, j) = cache.betan * bb;
                }
            }
            ier = decc(integrator.n, cache.e2r, cache.e2i, cache.ip2);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // This option is not provided
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            for (int j = 0; j < integrator.n - 1; j++) {
                cache.e2r(j + 1, j) = -cache.dfdu(j + 1, j);
                cache.e2i(j + 1, j) = 0.0;
            }
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i <= j; i++) {
                    cache.e2i(i, j) = 0.0;
                    cache.e2r(i, j) = -cache.dfdu(i, j);
                }
                cache.e2r(j, j) += cache.alphn;
                cache.e2i(j, j) = cache.betan;
            }
            ier = dechc(integrator.n, cache.e2r, cache.e2i, 1, cache.ip2);
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < alg.nm1; i++) {
                    cache.e2r(i, j) = -cache.dfdu(i, j + alg.m1);
                    cache.e2i(i, j) = 0.0;
                }
                cache.e2r(j, j) += cache.alphn;
                cache.e2i(j, j) = cache.betan;
            }
            break;

        case (12):
            // mass = identity, Jacobian a banded matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e2r(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                    cache.e2i(i + cache.mle, j) = 0.0;
                }
                cache.e2r(cache.mdiag, j) += cache.alphn;
                cache.e2i(cache.mdiag, j) += cache.betan;
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < alg.nm1; i++) {
                    cache.e2r(i, j) = -cache.dfdu(i, j + alg.m1);
                    cache.e2i(i, j) = 0.0;
                }
                for (int i = std::max(0, j - alg.mumas); i < std::min(alg.nm1, j + alg.mlmas + 1); i++) {
                    ffma = cache.mass_matrix(i - j + cache.mbdiag - 1, j);
                    cache.e2r(j, j) += cache.alphn * ffma;
                    cache.e2i(j, j) += cache.betan * ffma;
                }
            }
            break;

        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e2r(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                    cache.e2i(i + cache.mle, j) = 0.0;
                }
                for (int i = 0; i < cache.mbb; i++) {
                    ffma = cache.mass_matrix(i, j);
                    cache.e2r(i + cache.mdiff, j) += cache.alphn * ffma;
                    cache.e2i(i + cache.mdiff, j) += cache.betan * ffma;
                }
            }
            break;

        case (15):
            // mass is a full matrix, Jacobian a full matrix, second order
            for (int j = 0; j < alg.nm1; j++) {
                for (int i = 0; i < alg.nm1; i++) {
                    cache.e2r(i, j) = cache.alphn * cache.mass_matrix(i, j) - cache.dfdu(i, j + alg.m1);
                    cache.e2i(i, j) = cache.betan * cache.mass_matrix(i, j);
                }
            }
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

    }

    switch (cache.ijob) {
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
            abno = cache.alphn * cache.alphn + cache.betan * cache.betan;
            alp = cache.alphn / abno;
            bet = cache.betan / abno;
            for (int j = 0; j < alg.m2; j++) {
                for (int i = 0; i < alg.nm1; i++) {
                    sumr = sumi = 0.0;
                    for (int k = 0; k < mm; k++) {
                        sums = sumr + cache.dfdu(i, j + k * alg.m2);
                        sumr = sums * alp + sumi * bet;
                        sumi = sumi * alp - sums * bet;
                    }
                    cache.e2r(i, j) -= sumr;
                    cache.e2i(i, j) -= sumi;
                }
            }
            ier = decc(alg.nm1, cache.e2r, cache.e2i, cache.ip2);
            break;

        case (12):
        case (14):
            mm = alg.m1 / alg.m2;
            abno = cache.alphn * cache.alphn + cache.betan * cache.betan;
            alp = cache.alphn / abno;
            bet = cache.betan / abno;
            for (int j = 0; j < alg.m2; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    sumr = sumi = 0.0;
                    for (int k = 0; k < mm; k++) {
                        sums = sumr + cache.dfdu(i, j + k * alg.m2);
                        sumr = sums * alp + sumi * bet;
                        sumi = sumi * alp - sums * bet;
                    }
                    cache.e2r(i + cache.mle, j) -= sumr;
                    cache.e2i(i + cache.mle, j) -= sumi;
                }
            }
            ier = decbc(alg.nm1, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.ip2);
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;
    }

    return ier;

} // decomp_complex


static int linear_solve(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    int mm, mp, mp1, ii, jkm, mpi, ier = 0;
    double abno = 0, bb, e1imp, s1, s2, s3, sum1, sum2, sum3, sumh;
    double ffja, z2i, z3i, zsafe;

    switch (cache.ijob) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            sol(integrator.n, cache.e1, cache.z1, cache.ip1);
            solc(integrator.n, cache.e2r, cache.e2i, cache.z2, cache.z3, cache.ip2);
            break;

        case (2):
            // mass = identity, Jacobian a banded matrix
            for (int i = 0; i < integrator.n; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            solb(integrator.n, cache.e1, cache.mle, cache.mue, cache.z1, cache.ip1);
            solbc(integrator.n, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.z2, cache.z3, cache.ip2);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++) {
                s1 = s2 = s3 = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(integrator.n, i + alg.mumas + 1); j++) {
                    bb = cache.mass_matrix(i - j + cache.mbdiag - 1, j);
                    s1 -= bb * cache.f1[j];
                    s2 -= bb * cache.f2[j];
                    s3 -= bb * cache.f3[j];
                }
                cache.z1[i] += s1 * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            sol(integrator.n, cache.e1, cache.z1, cache.ip1);
            solc(integrator.n, cache.e2r, cache.e2i, cache.z2, cache.z3, cache.ip2);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int i = 0; i < integrator.n; i++) {
                s1 = s2 = s3 = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(integrator.n, i + alg.mumas + 1); j++) {
                    bb = cache.mass_matrix(i - j + cache.mbdiag - 1, j);
                    s1 -= bb * cache.f1[j];
                    s2 -= bb * cache.f2[j];
                    s3 -= bb * cache.f3[j];
                }
                cache.z1[i] += s1 * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            solb(integrator.n, cache.e1, cache.mle, cache.mue, cache.z1, cache.ip1);
            solbc(integrator.n, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.z2, cache.z3, cache.ip2);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++) {
                s1 = s2 = s3 = 0.0;
                for (int j = 0; j < integrator.n; j++) {
                    bb = cache.mass_matrix(i, j);
                    s1 -= bb * cache.f1[j];
                    s2 -= bb * cache.f2[j];
                    s3 -= bb * cache.f3[j];
                }
                cache.z1[i] += s1 * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            sol(integrator.n, cache.e1, cache.z1, cache.ip1);
            solc(integrator.n, cache.e2r, cache.e2i, cache.z2, cache.z3, cache.ip2);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // This option is not provided
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            for (int i = 0; i < integrator.n; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            for (int mm1 = integrator.n - 3; mm1 >= 0; mm1--) {
                mp = integrator.n - mm1 - 2;
                mp1 = mp - 1;
                ii = cache.iphes[mp];
                if (ii != mp) {
                    zsafe = cache.z1[mp];
                    cache.z1[mp] = cache.z1[ii];
                    cache.z1[ii] = zsafe;
                    zsafe = cache.z2[mp];
                    cache.z2[mp] = cache.z2[ii];
                    cache.z2[ii] = zsafe;
                    zsafe = cache.z3[mp];
                    cache.z3[mp] = cache.z3[ii];
                    cache.z3[ii] = zsafe;
                }
                for (int i = mp + 1; i < integrator.n; i++) {
                    e1imp = cache.dfdu(i, mp1);
                    cache.z1[i] -= e1imp * cache.z1[mp];
                    cache.z2[i] -= e1imp * cache.z2[mp];
                    cache.z3[i] -= e1imp * cache.z3[mp];
                }
            }
            solh(integrator.n, cache.e1, 1, cache.z1, cache.ip1);
            solhc(integrator.n, cache.e2r, cache.e2i, 1, cache.z2, cache.z3, cache.ip2);
            for (int mm1 = 0; mm1 < integrator.n - 2; mm1++) {
                mp = integrator.n - mm1 - 2;
                mp1 = mp - 1;
                for (int i = mp; i < integrator.n; i++) {
                    e1imp = cache.dfdu(i, mp1);
                    cache.z1[i] += e1imp * cache.z1[mp];
                    cache.z2[i] += e1imp * cache.z2[mp];
                    cache.z3[i] += e1imp * cache.z3[mp];
                }
                ii = cache.iphes[mp];
                if (ii != mp) {
                    zsafe = cache.z1[mp];
                    cache.z1[mp] = cache.z1[ii];
                    cache.z1[ii] = zsafe;
                    zsafe = cache.z2[mp];
                    cache.z2[mp] = cache.z2[ii];
                    cache.z2[ii] = zsafe;
                    zsafe = cache.z3[mp];
                    cache.z3[mp] = cache.z3[ii];
                    cache.z3[ii] = zsafe;
                }
            }
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
        case (12):
            // ---  b = identity, Jacobian a banded matrix, second order
            for (int i = 0; i < integrator.n; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int i = 0; i < alg.m1; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            for (int i = 0; i < alg.nm1; i++) {
                s1 = s2 = s3 = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(alg.nm1, i + alg.mumas + 1); j++) {
                    bb = cache.mass_matrix(i - j + cache.mbdiag - 1, j);
                    s1 -= bb * cache.f1[j + alg.m1];
                    s2 -= bb * cache.f2[j + alg.m1];
                    s3 -= bb * cache.f3[j + alg.m1];
                }
                cache.z1[i + alg.m1] += s1 * cache.fac1;
                cache.z2[i + alg.m1] = cache.z2[i + alg.m1] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i + alg.m1] = cache.z3[i + alg.m1] + s3 * cache.alphn + s2 * cache.betan;
            }
            break;

        case (15):
            // mass is a full matrix, Jacobian a full matrix, second order
            for (int i = 0; i < alg.m1; i++) {
                s2 = -cache.f2[i];
                s3 = -cache.f3[i];
                cache.z1[i] -= cache.f1[i] * cache.fac1;
                cache.z2[i] = cache.z2[i] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i] = cache.z3[i] + s3 * cache.alphn + s2 * cache.betan;
            }
            for (int i = 0; i < alg.nm1; i++) {
                s1 = s2 = s3 = 0.0;
                for (int j = 0; j < alg.nm1; j++) {
                    bb = cache.mass_matrix(i, j);
                    s1 -= bb * cache.f1[j + alg.m1];
                    s2 -= bb * cache.f2[j + alg.m1];
                    s3 -= bb * cache.f3[j + alg.m1];
                }
                cache.z1[i + alg.m1] += s1 * cache.fac1;
                cache.z2[i + alg.m1] = cache.z2[i + alg.m1] + s2 * cache.alphn - s3 * cache.betan;
                cache.z3[i + alg.m1] = cache.z3[i + alg.m1] + s3 * cache.alphn + s2 * cache.betan;
            }
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;
    }

    switch (cache.ijob) {
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
            abno = cache.alphn * cache.alphn + cache.betan * cache.betan;
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                sum1 = sum2 = sum3 = 0.0;
                for (int k = mm - 1; k >= 0; k--) {
                    jkm = j + k * alg.m2;
                    sum1 = (cache.z1[jkm] + sum1) / cache.fac1;
                    sumh = (cache.z2[jkm] + sum2) / abno;
                    sum3 = (cache.z3[jkm] + sum3) / abno;
                    sum2 = sumh * cache.alphn + sum3 * cache.betan;
                    sum3 = sum3 * cache.alphn - sumh * cache.betan;
                    for (int i = 0; i < alg.nm1; i++) {
                        cache.z1[i + alg.m1] += cache.dfdu(i, jkm) * sum1;
                        cache.z2[i + alg.m1] += cache.dfdu(i, jkm) * sum2;
                        cache.z3[i + alg.m1] += cache.dfdu(i, jkm) * sum3;
                    }
                }
            }
            cache.z1_seg = cache.z1.segment(alg.m1, integrator.n);
            sol(alg.nm1, cache.e1, cache.z1_seg, cache.ip1);
            cache.z1.segment(alg.m1, integrator.n) = cache.z1_seg;

            cache.z2_seg = cache.z2.segment(alg.m1, integrator.n);
            cache.z3_seg = cache.z3.segment(alg.m1, integrator.n);
            solc(alg.nm1, cache.e2r, cache.e2i, cache.z2_seg, cache.z3_seg, cache.ip2);
            cache.z2.segment(alg.m1, integrator.n) = cache.z2_seg;
            cache.z3.segment(alg.m1, integrator.n) = cache.z3_seg;
            break;
        case (12):
        case (14):
            abno = cache.alphn * cache.alphn + cache.betan * cache.betan;
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                sum1 = sum2 = sum3 = 0.0;
                for (int k = mm - 1; k >= 0; k--) {
                    jkm = j + k * alg.m2;
                    sum1 = (cache.z1[jkm] + sum1) / cache.fac1;
                    sumh = (cache.z2[jkm] + sum2) / abno;
                    sum3 = (cache.z3[jkm] + sum3) / abno;
                    sum2 = sumh * cache.alphn + sum3 * cache.betan;
                    sum3 = sum3 * cache.alphn - sumh * cache.betan;
                    for (int i = std::max(0, j - alg.mujac); i < std::min(alg.nm1, j + alg.mljac + 1); i++) {
                        ffja = cache.dfdu(i + alg.mujac - j, jkm);
                        cache.z1[i + alg.m1] += ffja * sum1;
                        cache.z2[i + alg.m1] += ffja * sum2;
                        cache.z3[i + alg.m1] += ffja * sum3;
                    }
                }
            }

            cache.z1_seg = cache.z1.segment(alg.m1, integrator.n);
            solb(alg.nm1, cache.e1, cache.mle, cache.mue, cache.z1_seg, cache.ip1);
            cache.z1.segment(alg.m1, integrator.n) = cache.z1_seg;

            cache.z2_seg = cache.z2.segment(alg.m1, integrator.n);
            cache.z3_seg = cache.z3.segment(alg.m1, integrator.n);
            solbc(alg.nm1, cache.e2r, cache.e2i, cache.mle, cache.mue, cache.z2_seg, cache.z3_seg, cache.ip2);
            cache.z2.segment(alg.m1, integrator.n) = cache.z2_seg;
            cache.z3.segment(alg.m1, integrator.n) = cache.z3_seg;
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;
    }

    switch (cache.ijob) {
        case (1):
        case (2):
        case (3):
        case (4):
        case (5):
        case (7):
            break;

        case (11):
        case (12):
        case (13):
        case (14):
        case (15):
            for (int i = alg.m1 - 1; i >= 0; i--) {
                mpi = alg.m2 + i;
                cache.z1[i] = (cache.z1[i] + cache.z1[mpi]) / cache.fac1;
                z2i = cache.z2[i] + cache.z2[mpi];
                z3i = cache.z3[i] + cache.z3[mpi];
                cache.z3[i] = (z3i * cache.alphn - z2i * cache.betan) / abno;
                cache.z2[i] = (z2i * cache.alphn + z3i * cache.betan) / abno;
            }
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

    }

    return ier;

} // linear_solve


static int perform_decompositions(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    // compute the matrices cache.e1 and e2 and their decompositions
    cache.fac1 = alg.u1 / integrator.dt;
    cache.alphn = alg.alph / integrator.dt;
    cache.betan = alg.beta / integrator.dt;

    cache.ier = decomp_real(integrator, alg, cache);

    if (cache.ier != 0) {
        if (cache.ier == -1) {
            throw RadauLinearAlgebraError();
        }
        cache.num_singular++;
        if (cache.num_singular >= 5) {
            throw RadauSingularMatrixEncountered();
        }
        integrator.dt *= 0.5;
        cache.dtfac = 0.5;
        cache.reject = true;
        cache.last = false;
        if (!cache.caljac) {
            integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
        }
        return cache.ier;
    }

    cache.ier = decomp_complex(integrator, alg, cache);

    if (cache.ier != 0) {
        if (cache.ier == -1) throw RadauLinearAlgebraError();
        cache.num_singular++;
        if (cache.num_singular >= 5) {
            throw RadauSingularMatrixEncountered();
        }
        integrator.dt *= 0.5;
        cache.dtfac = 0.5;
        cache.reject = true;
        cache.last = false;
        if (!cache.caljac) integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
        return cache.ier;
    }
    integrator.num_decompositions++;
    return cache.ier;
}


static void perpare_newton_iteration(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    if (cache.first || alg.startn) {
        for (int i = 0; i < integrator.n; i++)
            cache.z1[i] = cache.z2[i] = cache.z3[i] = cache.f1[i] = cache.f2[i] = cache.f3[i] = 0.0;
    } else {
        double c3q = integrator.dt / integrator.dtprev;
        double c1q = alg.c1 * c3q;
        double c2q = alg.c2 * c3q;
        double ak1, ak2, ak3;

        for (int i = 0; i < integrator.n; i++) {
            ak1 = cache.cont[i + integrator.n];
            ak2 = cache.cont[i + 2 * integrator.n];
            ak3 = cache.cont[i + 3 * integrator.n];
            cache.z1[i] = c1q * (ak1 + (c1q - alg.c2m1) * (ak2 + (c1q - alg.c1m1) * ak3));
            cache.z2[i] = c2q * (ak1 + (c2q - alg.c2m1) * (ak2 + (c2q - alg.c1m1) * ak3));
            cache.z3[i] = c3q * (ak1 + (c3q - alg.c2m1) * (ak2 + (c3q - alg.c1m1) * ak3));
            cache.f1[i] = alg.ti11 * cache.z1[i] + alg.ti12 * cache.z2[i] + alg.ti13 * cache.z3[i];
            cache.f2[i] = alg.ti21 * cache.z1[i] + alg.ti22 * cache.z2[i] + alg.ti23 * cache.z3[i];
            cache.f3[i] = alg.ti31 * cache.z1[i] + alg.ti32 * cache.z2[i] + alg.ti33 * cache.z3[i];
        }
    }
}


static void perform_newton_iteration(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    cache.newt = 0;
    cache.faccon = pow(std::max(cache.faccon, std::numeric_limits<double>::epsilon()), 0.8);
    cache.theta = fabs(alg.thet);
    double dyno, dynold;

    while (true) {
        if (cache.newt >= alg.nit) {
            if (cache.ier != 0) {
                cache.num_singular++;
                if (cache.num_singular >= 5) {
                    throw RadauSingularMatrixEncountered();
                }
            }
            integrator.dt *= 0.5;
            cache.dtfac = 0.5;
            cache.reject = true;
            cache.last = false;
            if (!cache.caljac) {
                integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
            }
            cache.loop = true;
            break;
        }
        // compute the right-hand side
        for (int i = 0; i < integrator.n; i++) cache.cont[i] = integrator.u[i] + cache.z1[i];
        integrator.f->dudt(cache.z1, cache.cont, integrator.t + alg.c1 * integrator.dt);

        for (int i = 0; i < integrator.n; i++) cache.cont[i] = integrator.u[i] + cache.z2[i];
        integrator.f->dudt(cache.z2, cache.cont, integrator.t + alg.c2 * integrator.dt);

        for (int i = 0; i < integrator.n; i++) cache.cont[i] = integrator.u[i] + cache.z3[i];
        integrator.f->dudt(cache.z3, cache.cont, integrator.t + integrator.dt);

        integrator.num_function_evaluations += 3;

        // solve the linear systems
        for (int i = 0; i < integrator.n; i++) {
            double a1 = cache.z1[i];
            double a2 = cache.z2[i];
            double a3 = cache.z3[i];
            cache.z1[i] = alg.ti11 * a1 + alg.ti12 * a2 + alg.ti13 * a3;
            cache.z2[i] = alg.ti21 * a1 + alg.ti22 * a2 + alg.ti23 * a3;
            cache.z3[i] = alg.ti31 * a1 + alg.ti32 * a2 + alg.ti33 * a3;
        }
        cache.ier = linear_solve(integrator, alg, cache);
        if (cache.ier == -1) throw RadauLinearAlgebraError();
        integrator.num_linear_solves++;
        cache.newt++;
        dyno = 0.0;
        double denom;
        for (int i = 0; i < integrator.n; i++) {
            denom = cache.scal[i];
            dyno = dyno + pow(cache.z1[i] / denom, 2) + pow(cache.z2[i] / denom, 2) +
                    pow(cache.z3[i] / denom, 2);
        }
        dyno = sqrt(dyno / double(3 * integrator.n));
        // bad convergence or number of iterations to large
        if ((cache.newt > 1) && (cache.newt < alg.nit)) {
            double thq = dyno / dynold;
            if (cache.newt == 2) cache.theta = thq;
            else cache.theta = sqrt(thq * cache.thqold);
            cache.thqold = thq;
            if (cache.theta < 0.99) {
                cache.faccon = cache.theta / (1.0 - cache.theta);
                double dyth = cache.faccon * dyno * pow(cache.theta, alg.nit - 1 - cache.newt) / alg.fnewt;
                if (dyth >= 1.0) {
                    double qnewt = std::max(1.0e-4, std::min(20.0, dyth));
                    cache.dtfac = 0.8 * pow(qnewt, -1.0 / (4.0 + alg.nit - 1 - cache.newt));
                    integrator.dt *= cache.dtfac;
                    cache.reject = true;
                    cache.last = false;
                    if (cache.caljac) {
                        integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
                    }
                    cache.loop = true;
                    break;
                }
            } else {
                if (cache.ier != 0) {
                    cache.num_singular++;
                    if (cache.num_singular >= 5) {
                        throw RadauSingularMatrixEncountered();
                    }
                }
                integrator.dt *= 0.5;
                cache.dtfac = 0.5;
                cache.reject = true;
                cache.last = false;
                if (!cache.caljac) {
                    integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
                }
                cache.loop = true;
                break;
            }
        }
        dynold = std::max(dyno, std::numeric_limits<double>::epsilon());
        for (int i = 0; i < integrator.n; i++) {
            cache.f1[i] = cache.f1[i] + cache.z1[i];
            cache.f2[i] = cache.f2[i] + cache.z2[i];
            cache.f3[i] = cache.f3[i] + cache.z3[i];
            cache.z1[i] = alg.t11 * cache.f1[i] + alg.t12 * cache.f2[i] + alg.t13 * cache.f3[i];
            cache.z2[i] = alg.t21 * cache.f1[i] + alg.t22 * cache.f2[i] + alg.t23 * cache.f3[i];
            cache.z3[i] = alg.t31 * cache.f1[i] + cache.f2[i];
        }
        if (cache.faccon * dyno <= alg.fnewt) break;
    }
}


static int error(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    int mm = 0, ii, mp, ier = 0;
    double sum, zsafe;

    double hee1 = -(13.0 + 7.0 * sqrt(6.0)) / (3.0 * integrator.dt);
    double hee2 = (-13.0 + 7.0 * sqrt(6.0)) / (3.0 * integrator.dt);
    double hee3 = -1.0 / (3.0 * integrator.dt);

    switch (cache.ijob) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            sol(integrator.n, cache.e1, cache.cont, cache.ip1);
            break;

        case (2):
            // mass = identity, Jacobian a banded matrix
            for (int i = 0; i < integrator.n; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            solb(integrator.n, cache.e1, cache.mle, cache.mue, cache.cont, cache.ip1);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < integrator.n; i++) {
                sum = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(integrator.n, i + alg.mumas + 1); j++)
                    sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j) * cache.f1[j];
                cache.f2[i] = sum;
                cache.cont[i] = sum + cache.u0[i];
            }
            sol(integrator.n, cache.e1, cache.cont, cache.ip1);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int i = 0; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < integrator.n; i++) {
                sum = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(integrator.n, i + alg.mumas + 1); j++)
                    sum = sum + cache.mass_matrix(i - j + cache.mbdiag - 1, j) * cache.f1[j];
                cache.f2[i] = sum;
                cache.cont[i] = sum + cache.u0[i];
            }
            solb(integrator.n, cache.e1, cache.mle, cache.mue, cache.cont, cache.ip1);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int i = 0; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < integrator.n; i++) {
                sum = 0.0;
                for (int j = 0; j < integrator.n; j++)
                    sum += cache.mass_matrix(j, i) * cache.f1[j];
                cache.f2[i] = sum;
                cache.cont[i] = sum + cache.u0[i];
            }
            sol(integrator.n, cache.e1, cache.cont, cache.ip1);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // this option is not provided
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            for (int i = 0; i < integrator.n; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            for (int mm1 = integrator.n - 3; mm1 >= 0; mm1--) {
                mp = integrator.n - mm1 - 2;
                ii = cache.iphes[mp];
                if (ii != mp) {
                    zsafe = cache.cont[mp];
                    cache.cont[mp] = cache.cont[ii];
                    cache.cont[ii] = zsafe;
                }
                for (int i = mp; i < integrator.n; i++)
                    cache.cont[i] -= cache.dfdu(i, mp - 1) * cache.cont[mp];
            }
            solh(integrator.n, cache.e1, 1, cache.cont, cache.ip1);
            for (int mm1 = 0; mm1 < integrator.n - 2; mm1++) {
                mp = integrator.n - mm1 - 2;
                for (int i = mp; i < integrator.n; i++)
                    cache.cont[i] += cache.dfdu(i, mp - 1) * cache.cont[mp];
                ii = cache.iphes[mp];
                if (ii != mp) {
                    zsafe = cache.cont[mp];
                    cache.cont[mp] = cache.cont[ii];
                    cache.cont[ii] = zsafe;
                }
            }
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
            for (int i = 0; i < integrator.n; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            break;

        case (12):
            // mass = identity, Jacobian a banded matrix, second order
            for (int i = 0; i < integrator.n; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int i = 0; i < alg.m1; i++) {
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            for (int i = alg.m1; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < alg.nm1; i++) {
                sum = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(alg.nm1, i + alg.mumas + 1); j++)
                    sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j) * cache.f1[j + alg.m1];
                cache.f2[i + alg.m1] = sum;
                cache.cont[i + alg.m1] = sum + cache.u0[i + alg.m1];
            }
            break;

        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int i = 0; i < alg.m1; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            for (int i = alg.m1; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < alg.nm1; i++) {
                sum = 0.0;
                for (int j = std::max(0, i - alg.mlmas); j < std::min(alg.nm1, i + alg.mumas + 1); j++)
                    sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j) * cache.f1[j + alg.m1];
                cache.f2[i + alg.m1] = sum;
                cache.cont[i + alg.m1] = sum + cache.u0[i + alg.m1];
            }
            break;

        case (15):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int i = 0; i < alg.m1; i++) {
                cache.f2[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
                cache.cont[i] = cache.f2[i] + cache.u0[i];
            }
            for (int i = alg.m1; i < integrator.n; i++)
                cache.f1[i] = hee1 * cache.z1[i] + hee2 * cache.z2[i] + hee3 * cache.z3[i];
            for (int i = 0; i < alg.nm1; i++) {
                sum = 0.0;
                for (int j = 0; j < alg.nm1; j++)
                    sum += cache.mass_matrix(j, i) * cache.f1[j + alg.m1];
                cache.f2[i + alg.m1] = sum;
                cache.cont[i + alg.m1] = sum + cache.u0[i + alg.m1];
            }
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

    }

    switch (cache.ijob) {
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
                sum = 0.0;
                for (int k = mm - 1; k >= 0; k--) {
                    sum = (cache.cont[j + k * alg.m2] + sum) / cache.fac1;
                    for (int i = 0; i < alg.nm1; i++)
                        cache.cont[i + alg.m1] += cache.dfdu(i, j + k * alg.m2) * sum;
                }
            }
            cache.cont_seg = cache.cont.segment(alg.m1, cache.cont.size());
            sol(alg.nm1, cache.e1, cache.cont_seg, cache.ip1);
            cache.cont.segment(alg.m1, cache.cont.size()) = cache.cont_seg;

            for (int i = alg.m1 - 1; i >= 0; i--)
                cache.cont[i] = (cache.cont[i] + cache.cont[alg.m2 + i]) / cache.fac1;
            break;

        case (12):
        case (14):
            mm = alg.m1 / alg.m2;
            for (int j = 0; j < alg.m2; j++) {
                sum = 0.0;
                for (int k = mm - 1; k >= 0; k--) {
                    sum = (cache.cont[j + k * alg.m2] + sum) / cache.fac1;
                    for (int i = std::max(0, j - alg.mujac); i < std::min(alg.nm1, j + alg.mljac); i++)
                        cache.cont[i + alg.m1] += cache.dfdu(i + alg.mujac - j, j + k * alg.m2) * sum;
                }
            }
            cache.cont_seg = cache.cont.segment(alg.m1, cache.cont.size());
            solb(alg.nm1, cache.e1, cache.mle, cache.mue, cache.cont_seg, cache.ip1);
            cache.cont.segment(alg.m1, cache.cont.size()) = cache.cont_seg;

            for (int i = alg.m1 - 1; i >= 0; i--)
                cache.cont[i] = (cache.cont[i] + cache.cont[alg.m2 + i]) / cache.fac1;
            break;
        default:
            std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
            ier = -1;
            return ier;

    }

    cache.err = 0.0;
    for (int i = 0; i < integrator.n; i++)
        cache.err += pow(cache.cont[i] / cache.scal[i], 2);
    cache.err = std::max(sqrt(cache.err / integrator.n), 1.0e-10);

    if (cache.err < 1.0) return ier;

    if (cache.first || cache.reject) {
        for (int i = 0; i < integrator.n; i++) cache.cont[i] = integrator.u[i] + cache.cont[i];
        integrator.f->dudt(cache.f1, cache.cont, integrator.t);
        integrator.num_function_evaluations++;
        for (int i = 0; i < integrator.n; i++) cache.cont[i] = cache.f1[i] + cache.f2[i];

        switch (cache.ijob) {
            case (1):
            case (3):
            case (5):
                // full matrix option
                sol(integrator.n, cache.e1, cache.cont, cache.ip1);
                break;

            case (2):
            case (4):
                // banded matrix option
                solb(integrator.n, cache.e1, cache.mle, cache.mue, cache.cont, cache.ip1);
                break;

            case (7):
                //Hessenberg matrix option
                // mass = identity, Jacobian a full matrix, Hessenberg-option
                for (int mm1 = integrator.n - 3; mm1 >= 0; mm1--) {
                    mp = integrator.n - mm1 - 2;
                    ii = cache.iphes[mp];
                    if (ii != mp) {
                        zsafe = cache.cont[mp];
                        cache.cont[mp] = cache.cont[ii];
                        cache.cont[ii] = zsafe;
                    }
                    for (int i = mp; i < integrator.n; i++)
                        cache.cont[i] -= cache.dfdu(i, mp - 1) * cache.cont[mp];
                }
                solh(integrator.n, cache.e1, 1, cache.cont, cache.ip1);
                for (int mm1 = 0; mm1 < integrator.n - 2; mm1++) {
                    mp = integrator.n - mm1 - 2;
                    for (int i = mp; i < integrator.n; i++)
                        cache.cont[i] += cache.dfdu(i, mp - 1) * cache.cont[mp];
                    ii = cache.iphes[mp];
                    if (ii != mp) {
                        zsafe = cache.cont[mp];
                        cache.cont[mp] = cache.cont[ii];
                        cache.cont[ii] = zsafe;
                    }
                }
                break;

            case (11):
            case (13):
            case (15):
                // Full matrix option, second order
                for (int j = 0; j < alg.m2; j++) {
                    sum = 0.0;
                    for (int k = mm - 1; k >= 0; k--) {
                        sum = (cache.cont[j + k * alg.m2] + sum) / cache.fac1;
                        for (int i = 0; i < alg.nm1; i++)
                            cache.cont[i + alg.m1] += cache.dfdu(i, j + k * alg.m2) * sum;
                    }
                }

                cache.cont_seg = cache.cont.segment(alg.m1, cache.cont.size());
                sol(alg.nm1, cache.e1, cache.cont_seg, cache.ip1);
                cache.cont.segment(alg.m1, cache.cont.size()) = cache.cont_seg;

                for (int i = alg.m1 - 1; i >= 0; i--)
                    cache.cont[i] = (cache.cont[i] + cache.cont[alg.m2 + i]) / cache.fac1;
                break;

            case (12):
            case (14):
                // Banded matrix option, second order
                for (int j = 0; j < alg.m2; j++) {
                    sum = 0.0;
                    for (int k = mm - 1; k >= 0; k--) {
                        sum = (cache.cont[j + k * alg.m2] + sum) / cache.fac1;
                        for (int i = std::max(0, j - alg.mujac); i < std::min(alg.nm1, j + alg.mljac); i++)
                            cache.cont[i + alg.m1] += cache.dfdu(i + alg.mujac - j, j + k * alg.m2) * sum;
                    }
                }
                cache.cont_seg = cache.cont.segment(alg.m1, cache.cont.size());
                solb(alg.nm1, cache.e1, cache.mle, cache.mue, cache.cont_seg, cache.ip1);
                cache.cont.segment(alg.m1, cache.cont.size()) = cache.cont_seg;

                for (int i = alg.m1 - 1; i >= 0; i--)
                    cache.cont[i] = (cache.cont[i] + cache.cont[alg.m2 + i]) / cache.fac1;
                break;
            default:
                std::cout << "Not a value cache.ijob. \tcache.ijob = " << cache.ijob << std::endl;
                ier = -1;
                return ier;
        }

        cache.err = 0.0;
        for (int i = 0; i < integrator.n; i++)
            cache.err += pow(cache.cont[i] / cache.scal[i], 2);
        cache.err = std::max(sqrt(cache.err / integrator.n), 1.0e-10);
    }

    return ier;

} // error


static double
dense_output(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache, int t_i, double t_t, double /*t_dt*/) {
    double s = (t_t - integrator.t) / integrator.dtprev;

    return (cache.cont[t_i] +
            s * (cache.cont[t_i + integrator.n] +
                    (s - alg.c2m1) * (cache.cont[t_i + 2 * integrator.n] +
                            (s - alg.c1m1) * cache.cont[t_i + 3 * integrator.n])));
}


static ODESolution core_integrator(ODEIntegrator &integrator, Radau5 &alg, Radau5Cache &cache) {
    init(integrator, alg, cache);
    ODESolution solution{};

    integrator.dt = std::min(fabs(integrator.dt), integrator.opts.dtmax);
    integrator.dt = copysign(integrator.dt, integrator.tdir);
    integrator.dtprev = integrator.dt;

    cache.last = false;

    if ((integrator.t + integrator.dt * 1.0001 - integrator.t_final) * integrator.tdir >= 0.0) {
        integrator.dt = integrator.t_final - integrator.t;
        cache.last = true;
    }

    cache.dtopt = integrator.dt;
    cache.faccon = 1.0;

    if (integrator.opts.dense) {
        for (int i = 0; i < integrator.n; i++) {
            cache.cont[i] = integrator.u[i];
        }
    }
    solution.ts.push_back(integrator.t);
    solution.us.push_back(integrator.u);

    for (int i = 0; i < integrator.n; i++) {
        cache.scal[i] = integrator.opts.abstol + integrator.opts.reltol * fabs(integrator.u[i]);
    }

    integrator.f->dudt(cache.u0, integrator.u, integrator.t);
    integrator.num_function_evaluations++;

    cache.dtfac = integrator.dt;
    cache.ier = 0;

    // basic integration step
    integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);

    do {
        cache.loop = false;

        do {
            cache.ier = perform_decompositions(integrator, alg, cache);
        } while (cache.ier != 0);

        while (true) {
            integrator.num_steps++;
            if (integrator.num_steps >= integrator.opts.max_num_steps) {
                std::cout << " exit of RADAU5 at integrator.t = " << integrator.t << std::endl;
                std::cout << " more than integrator.opts.max_num_steps = " << integrator.opts.max_num_steps
                          << " steps are needed" << std::endl;
                solution.retcode = Retcode::MaxIters;
                return solution;
            }

            if (0.1 * fabs(integrator.dt) <= fabs(integrator.t) * std::numeric_limits<double>::epsilon()) {
                std::cout << " exit of RADAU5 at t = " << integrator.t << std::endl;
                std::cout << " step size too small, dt = " << integrator.dt << std::endl;
                solution.retcode = Retcode::DtLessThanMin;
                return solution;
            }

            // check the index of the problem
            if (alg.nind2 != 0) { // is index 2
                for (int i = alg.nind1; i < alg.nind1 + alg.nind2; i++)
                    cache.scal[i] = cache.scal[i] / cache.dtfac;
            }

            if (alg.nind3 != 0) { // is index 3
                for (int i = alg.nind1 + alg.nind2; i < alg.nind1 + alg.nind2 + alg.nind3; i++)
                    cache.scal[i] = cache.scal[i] / (cache.dtfac * cache.dtfac);
            }

            perpare_newton_iteration(integrator, alg, cache);
            perform_newton_iteration(integrator, alg, cache);

            if (cache.loop) {
                break;
            }

            // error estimation
            cache.err = 0.0;
            cache.ier = error(integrator, alg, cache);
            if (cache.ier == -1) {
                //throw RadauLinearAlgebraError();
                solution.retcode = LinearAlgError;
                return solution;
            }

            // computation of m_dtNew -- require 0.2 <= m_dtNew/integrator.dt <= 8.
            double fac = std::min(alg.safe, cache.cfac / (cache.newt + 2 * alg.nit));
            double quot = std::max(alg.facr, std::min(alg.facl, pow(cache.err, 0.25) / fac));
            double m_dtNew = integrator.dt / quot;

            //  is the error small enough ?
            if (cache.err < 1.0) {
                // step is accepted
                cache.first = false;
                integrator.num_accept++;
                if (alg.pred) {
                    // Use the predictive controller of Gustafsson
                    if (integrator.num_accept > 1) {
                        double facgus =
                                (cache.dtacc / (integrator.dt)) * pow(cache.err * cache.err / cache.erracc, 0.25) /
                                        alg.safe;
                        facgus = std::max(alg.facr, std::min(alg.facl, facgus));
                        quot = std::max(quot, facgus);
                        m_dtNew = integrator.dt / quot;
                    }
                    cache.dtacc = integrator.dt;
                    cache.erracc = std::max(1.0e-2, cache.err);
                }
                integrator.tprev = integrator.t;
                integrator.dtprev = integrator.dt;
                integrator.t += integrator.dt;
                double ak, acont3;
                for (int i = 0; i < integrator.n; i++) {
                    integrator.u[i] = integrator.u[i] + cache.z3[i];
                    cache.cont[i + integrator.n] = (cache.z2[i] - cache.z3[i]) / alg.c2m1;
                    ak = (cache.z1[i] - cache.z2[i]) / alg.c1mc2;
                    acont3 = cache.z1[i] / alg.c1;
                    acont3 = (ak - acont3) / alg.c2;
                    cache.cont[i + 2 * integrator.n] = (ak - cache.cont[i + integrator.n]) / alg.c1m1;
                    cache.cont[i + 3 * integrator.n] = cache.cont[i + 2 * integrator.n] - acont3;
                }

                for (int i = 0; i < integrator.n; i++) {
                    cache.scal[i] = integrator.opts.abstol + integrator.opts.reltol * fabs(integrator.u[i]);
                }
                if (integrator.opts.dense) {
                    for (int i = 0; i < integrator.n; i++) {
                        cache.cont[i] = integrator.u[i];
                    }
                }
                solution.ts.push_back(integrator.t);
                solution.us.push_back(integrator.u);

                cache.caljac = false;
                if (cache.last) {
                    integrator.dt = cache.dtopt;
                    solution.retcode = Success;
                    return solution;
                }
                integrator.f->dudt(cache.u0, integrator.u, integrator.t);
                integrator.num_function_evaluations++;
                m_dtNew = integrator.tdir * std::min(fabs(m_dtNew), integrator.opts.dtmax);
                cache.dtopt = std::min(integrator.dt, m_dtNew);
                if (cache.reject) m_dtNew = integrator.tdir * std::min(fabs(m_dtNew), fabs(integrator.dt));
                cache.reject = false;
                if ((integrator.t + m_dtNew / alg.quot1 - integrator.t_final) * integrator.tdir >= 0.0) {
                    integrator.dt = integrator.t_final - integrator.t;
                    cache.last = true;
                } else {
                    double qt = m_dtNew / (integrator.dt);
                    cache.dtfac = integrator.dt;
                    if ((cache.theta <= alg.thet) && (qt >= alg.quot1) && (qt <= alg.quot2)) {
                        continue;
                    }
                    integrator.dt = m_dtNew;
                }
                cache.dtfac = integrator.dt;
                if (cache.theta > alg.thet) {
                    integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
                }
                cache.loop = true;
            } else {
                // step is rejected
                cache.reject = true;
                cache.last = false;
                if (cache.first) {
                    integrator.dt *= 0.1;
                    cache.dtfac = 0.1;
                } else {
                    cache.dtfac = m_dtNew / integrator.dt;
                    integrator.dt = m_dtNew;
                }
                if (integrator.num_accept >= 1) integrator.num_reject++;
                if (!cache.caljac) integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
                cache.loop = true;
            }
            break;
        }
    } while (cache.loop);
    solution.retcode = Success;
    return solution;
}


ODESolution solve(ODEProblem &prob, Radau5 &alg, ODEIntegratorOptions &opts) {
    ODEIntegrator integrator{prob, opts};
    Radau5Cache cache{};
    return core_integrator(integrator, alg, cache);
}

ODESolution solve(ODEProblem &prob, Radau5 &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}

}
}

#endif //LANRE_DIFFEQ_RADAU_HPP
