//
// Created by Logan Morrison on 3/16/20.
//

#ifndef LANRE_DIFFEQ_RODAS_HPP
#define LANRE_DIFFEQ_RODAS_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/algorithm.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/solution.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "lanre/diffeq/decsol.hpp"

namespace lanre {
namespace diffeq {

struct Rodas : public ODEAlgorithm {
    bool hess = false;
    bool implicit = false, autonomous = false;
    int npred = 0;
    int m1 = 0, m2 = 0;
    double fac1 = 5.0, fac2 = 1.0 / 6.0;
    double safe = 0.9;

    int mujac = -1, mljac = -1;
    int mumas = -1, mlmas = -1;

    double c2 = .386, c3 = .21, c4 = .63;
    double bet2p = .0317, bet3p = .0635, bet4p = .3438;
    double d1 = .25, d2 = -.1043, d3 = .1035, d4 = -.03620000000000023;
    double a21 = 1.544;
    double a31 = .9466785280815826, a32 = .2557011698983284;
    double a41 = 3.314825187068521, a42 = 2.896124015972201, a43 = .9986419139977817;
    double a51 = 1.221224509226641, a52 = 6.019134481288629, a53 = 12.53708332932087;
    double a54 = -.687886036105895;
    double c21 = -5.6688, c31 = -2.430093356833875, c32 = -.2063599157091915;
    double c41 = -.1073529058151375, c42 = -9.594562251023355, c43 = -20.47028614809616;
    double c51 = 7.496443313967647, c52 = -10.24680431464352, c53 = -33.99990352819905;
    double c54 = 11.7089089320616;
    double c61 = 8.083246795921522, c62 = -7.981132988064893, c63 = -31.52159432874371;
    double c64 = 16.31930543123136, c65 = -6.058818238834054;
    double gamma = .25;
    double d21 = 10.12623508344586, d22 = -7.487995877610167, d23 = -34.80091861555747;
    double d24 = -7.992771707568823, d25 = 1.025137723295662;
    double d31 = -.6762803392801253, d32 = 6.087714651680015, d33 = 16.43084320892478;
    double d34 = 24.76722511418386, d35 = -6.594389125716872;

    Rodas() = default;

    Rodas(int method);
};

struct RodasCache : public ODEAlgorithmCache {
    double dtopt, dtacc, erracc;
    double fac;
    double hd1, hd2, hd3, hd4;
    bool calhes;

    int job;
    bool pred;
    int nm1;
    int mle, mue;
    int mdiff, mbjac, mbb, mdiag, mbdiag;
    bool jband;
    int ldmas, ldjac, lde;

    Vector<double> ak1, ak2, ak3, ak4, ak5, ak6;
    Vector<double> du, du1, unew;
    Vector<double> dfdt;
    Matrix<double> dfdu, e, mass_matrix;
    Vector<int> ip, iphes;
    Vector<double> cont;

    // Counting variables
    int num_singular;

    bool last;
    bool reject;
};


Rodas::Rodas(int method) {
    if (method == 1) {
        return;
    } else if (method == 2) {
        c2 = .3507221;
        c3 = .2557041;
        c4 = .681779;
        bet2p = .0317;
        bet3p = .0047369;
        bet4p = .3438;
        d1 = .25;
        d2 = -.06902209999999998;
        d3 = -9.671999999999459e-4;
        d4 = -.08797900000000025;
        a21 = 1.4028884;
        a31 = .6581212688557198;
        a32 = -1.320936088384301;
        a41 = 7.131197445744498;
        a42 = 16.02964143958207;
        a43 = -5.561572550509766;
        a51 = 22.73885722420363;
        a52 = 67.38147284535289;
        a53 = -31.2187749303856;
        a54 = .7285641833203814;
        c21 = -5.1043536;
        c31 = -2.899967805418783;
        c32 = 4.040399359702244;
        c41 = -32.64449927841361;
        c42 = -99.35311008728094;
        c43 = 49.99119122405989;
        c51 = -76.46023087151691;
        c52 = -278.5942120829058;
        c53 = 153.9294840910643;
        c54 = 10.97101866258358;
        c61 = -76.29701586804983;
        c62 = -294.2795630511232;
        c63 = 162.0029695867566;
        c64 = 23.6516690309527;
        c65 = -7.652977706771382;
        gamma = .25;
        d21 = -38.71940424117216;
        d22 = -135.8025833007622;
        d23 = 64.51068857505875;
        d24 = -4.192663174613162;
        d25 = -2.53193205033506;
        d31 = -14.99268484949843;
        d32 = -76.30242396627033;
        d33 = 58.65928432851416;
        d34 = 16.61359034616402;
        d35 = -.6758691794084156;
    } else {
        //  Coefficients for RODAS with order 4 for linear parabolic problems
        // Gerd Steinebach (1993)
        gamma = .25;
        c2 = gamma * 3.0;
        c3 = .21;
        c4 = .63;
        bet2p = 0.;
        bet3p = c3 * c3 * (c3 / 6.0 - gamma / 2.0) / (gamma * gamma);
        bet4p = .3438;
        d1 = .25;
        d2 = -.5;
        d3 = -.023504;
        d4 = -.0362;
        a21 = 3.;
        a31 = 1.831036793486759;
        a32 = .4955183967433795;
        a41 = 2.304376582692669;
        a42 = -.05249275245743001;
        a43 = -1.176798761832782;
        a51 = -7.170454962423024;
        a52 = -4.741636671481785;
        a53 = -16.31002631330971;
        a54 = -1.062004044111401;
        c21 = -12.;
        c31 = -8.791795173947035;
        c32 = -2.207865586973518;
        c41 = 10.81793056857153;
        c42 = 6.780270611428266;
        c43 = 19.5348594464241;
        c51 = 34.19095006749676;
        c52 = 15.49671153725963;
        c53 = 54.7476087596413;
        c54 = 14.16005392148534;
        c61 = 34.62605830930532;
        c62 = 15.30084976114473;
        c63 = 56.99955578662667;
        c64 = 18.40807009793095;
        c65 = -5.714285714285717;
        d21 = 25.09876703708589;
        d22 = 11.62013104361867;
        d23 = 28.49148307714626;
        d24 = -5.664021568594133;
        d25 = 0.;
        d31 = 1.638054557396973;
        d32 = -.7373619806678748;
        d33 = 8.47791821923899;
        d34 = 15.9925314877952;
        d35 = -1.882352941176471;
    }
}


static void init(ODEIntegrator &integrator, Rodas &alg, RodasCache &cache) {
    cache.pred = alg.npred < 1;

    // Parameters for second order parameters.
    cache.nm1 = integrator.n - alg.m1;
    if (alg.m1 == 0) {
        alg.m2 = integrator.n;
    }
    if (alg.m2 == 0) {
        alg.m2 = alg.m1;
    }
    if (alg.m1 < 0 || alg.m2 < 0 || alg.m1 + alg.m2 > integrator.n) {
        throw std::runtime_error(
                "Curious inputs for m1,m2" +
                        std::to_string(alg.m1) + ", " + std::to_string(alg.m2) + ".");
    }

    if (integrator.opts.dtmax == 0.0) {
        integrator.opts.dtmax = integrator.t_final - integrator.t;
    }

    if (alg.mljac < 0) {
        alg.mljac = integrator.n;
        cache.jband = false;
    }
    if (alg.mujac < 0) {
        alg.mujac = integrator.n;
        cache.jband = false;
    }

    if (alg.mlmas < 0) {
        alg.mlmas = integrator.n;
    }
    if (alg.mumas < 0) {
        alg.mumas = integrator.n;
    }

    cache.jband = alg.mljac < cache.nm1;

    // Computation of the row-dimensions of the 2-d arrays
    // Jacobian and matrix e
    if (cache.jband) {
        cache.ldjac = alg.mljac + alg.mujac + 1;
        cache.lde = alg.mljac + cache.ldjac;
    } else {
        alg.mljac = cache.nm1;
        alg.mujac = cache.nm1;
        cache.ldjac = cache.nm1;
        cache.lde = cache.nm1;
    }
    // Mass matrix
    if (alg.implicit) {
        if (alg.mlmas != cache.nm1) {
            cache.ldmas = alg.mlmas + alg.mumas + 1;
            cache.job = cache.jband ? 4 : 3;
        } else {
            cache.ldmas = cache.nm1;
            cache.job = 5;
        }
        if (alg.mlmas > alg.mljac || alg.mumas > alg.mujac) {
            throw std::runtime_error("Rodas: Bandwidth of mass not larged than bandwidth of jacobian.");
        }
    } else {
        cache.ldmas = 0;
        cache.job = cache.jband ? 2 : 1;
    }
    cache.ldmas = std::max(1, cache.ldmas);

    // Resize the working arrays
    cache.du1.resize(integrator.n);
    cache.du.resize(integrator.n);
    cache.ak1.resize(integrator.n);
    cache.ak2.resize(integrator.n);
    cache.ak3.resize(integrator.n);
    cache.ak4.resize(integrator.n);
    cache.ak5.resize(integrator.n);
    cache.ak6.resize(integrator.n);
    cache.dfdt.resize(integrator.n);
    cache.cont.resize(4 * integrator.n);
    cache.ip.resize(integrator.n);

    cache.dfdu.resize(cache.ldjac, integrator.n);
    cache.mass_matrix.resize(cache.ldjac, integrator.n);
    cache.e.resize(cache.lde, integrator.n);

    // Additional preperations

    // Copy mass matrix
    if (alg.implicit) {
        for (int j = 0; j < integrator.n; j++) {
            for (int i = 0; i < alg.mlmas; i++) {
                cache.mass_matrix(i, j) = integrator.f->mass_matrix(i, j);
            }
        }
    }

    if (alg.m1 > 0) {
        cache.job += 10;
    }

    if (std::abs(integrator.dt) <= 10.0 * std::numeric_limits<double>::epsilon()) {
        integrator.dt = 1e-6;
    }
    integrator.dt = std::min(std::abs(integrator.dt), integrator.opts.dtmax);
    integrator.dt = std::copysign(integrator.dt, integrator.tdir);
    cache.reject = false;
    cache.last = false;
    cache.num_singular = 0;

    if (alg.autonomous) {
        cache.hd1 = 0.0;
        cache.hd2 = 0.0;
        cache.hd3 = 0.0;
        cache.hd4 = 0.0;
    }
    // Prepare bandwidths
    cache.mbdiag = alg.mumas + 1;
    if (cache.jband) {
        cache.mle = alg.mljac;
        cache.mue = alg.mujac;
        cache.mbjac = alg.mljac + alg.mujac + 1;
        cache.mbb = alg.mlmas + alg.mumas + 1;
        cache.mdiag = cache.mle + cache.mue + 1;
        cache.mdiff = cache.mle + cache.mue - alg.mumas;
    }
}


static int decomp_real(ODEIntegrator &integrator, Rodas &alg, RodasCache &cache) {
    int mm, ier = 0;
    double sum;
    //  cache.dfdu\[(.*?)\]\[(.*?)\] -> cache.e($1,$2)

    switch (cache.job) {
        case (1):
            // mass = identity, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.e(i, j) = -cache.dfdu(i, j);
                }
                cache.e(j, j) += cache.fac;
            }
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (2):
            // mass = identity, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++) {
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j);
                }
                cache.e(cache.mdiag, j) += cache.fac;
            }
            ier = decb(integrator.n, cache.e, cache.mle, cache.mue, cache.ip);
            break;

        case (3):
            // mass is a banded matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < integrator.n; i++)
                    cache.e(i, j) = -cache.dfdu(i, j);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e(i, j) += cache.fac * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (4):
            // mass is a banded matrix, Jacobian a banded matrix
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e(i + cache.mdiff, j) += cache.fac * cache.mass_matrix(i, j);
            }
            ier = decb(integrator.n, cache.e, cache.mle, cache.mue, cache.ip);
            break;

        case (5):
            // mass is a full matrix, Jacobian a full matrix
            for (int j = 0; j < integrator.n; j++)
                for (int i = 0; i < integrator.n; i++)
                    cache.e(i, j) = cache.mass_matrix(i, j) * cache.fac - cache.dfdu(i, j);
            ier = dec(integrator.n, cache.e, cache.ip);
            break;

        case (6):
            // mass is a full matrix, Jacobian a banded matrix
            // This option is not provided
            std::cout << "Not a valid job. job = " << cache.job << std::endl;
            ier = -1;
            return ier;

        case (7):
            // mass = identity, Jacobian a full matrix, Hessenberg-option
            if (cache.calhes) elmhes(integrator.n, 0, integrator.n, cache.dfdu, cache.iphes);
            cache.calhes = false;
            for (int j = 0; j < integrator.n - 1; j++) cache.e(j + 1, j) = -cache.dfdu(j + 1, j);
            for (int j = 0; j < integrator.n; j++) {
                for (int i = 0; i <= j; i++) cache.e(i, j) = -cache.dfdu(i, j);
                cache.e(j, j) += cache.fac;
            }
            ier = dech(integrator.n, cache.e, 1, cache.ip);
            break;

        case (11):
            // mass = identity, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.nm1; i++) {
                    cache.e(i, j) = -cache.dfdu(i, j + alg.m1);
                }
                cache.e(j, j) += cache.fac;
            }
            break;

        case (12):
            // mass = identity, Jacobian a banded matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                cache.e(cache.mdiag, j) += cache.fac;
            }
            break;

        case (13):
            // mass is a banded matrix, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.nm1; i++)
                    cache.e(i, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = std::max(0, j - alg.mumas); i < std::min(integrator.n, j + alg.mlmas + 1); i++)
                    cache.e(i, j) += cache.fac * cache.mass_matrix(i - j + cache.mbdiag - 1, j);
            }
            break;

        case (14):
            // mass is a banded matrix, Jacobian a banded matrix, second order
            for (int j = 0; j < cache.nm1; j++) {
                for (int i = 0; i < cache.mbjac; i++)
                    cache.e(i + cache.mle, j) = -cache.dfdu(i, j + alg.m1);
                for (int i = 0; i < cache.mbb; i++)
                    cache.e(i + cache.mdiff, j) += cache.fac * cache.mass_matrix(i, j);
            }
            break;

        case (15):
            // mass is a full matrix, Jacobian a full matrix, second order
            for (int j = 0; j < cache.nm1; j++)
                for (int i = 0; i < cache.nm1; i++)
                    cache.e(i, j) = cache.mass_matrix(i, j) * cache.fac - cache.dfdu(i, j + alg.m1);
            break;
        default:
            std::cout << "Not a value cache.job. \tcache.job = " << cache.job << std::endl;
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
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / cache.fac;
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
                        sum = (sum + cache.dfdu(i, j + k * alg.m2)) / cache.fac;
                    cache.e(i + cache.mle, j) -= sum;
                }
            }
            ier = decb(cache.nm1, cache.e, cache.mle, cache.mue, cache.ip);
            break;
        default:
            std::cout << "Not a value cache.job. \tcache.job = " << cache.job << std::endl;
            ier = -1;
            return ier;
    }

    return ier;

}


static int linearSolve(
        ODEIntegrator &integrator,
        Rodas &alg,
        RodasCache &cache,
        Vector<double> &du,
        Vector<double> &ak,
        Vector<double> &unew,
        double hd,
        bool stage1
) {
    ak = du + hd * cache.dfdt;

    if ((cache.job == 1 || cache.job == 11 || cache.job == 2 || cache.job == 12) && stage1) {
        for (int i = 0; i < integrator.n; i++) {
            ak(i) += unew(i);
        }
    }
    if ((cache.job == 13 || cache.job == 14) && stage1) {
        for (int i = 0; i < integrator.n; i++) {
            ak(i) += unew(i);
        }
        for (int i = 1; i <= cache.nm1; i++) {
            double sum = 0.0;
            for (int j = std::max(1, i - alg.mlmas); i <= std::min(cache.nm1, i + alg.mumas); i++) {
                sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j - 1) * unew(j + alg.m1 - 1);
            }
            ak(i + alg.m1 - 1) += sum;
        }
    }
    if (cache.job == 15 && stage1) {
        for (int i = 0; i < integrator.n; i++) {
            ak(i) += unew(i);
        }
        for (int i = 1; i <= cache.nm1; i++) {
            double sum = 0;
            for (int j = 1; j <= cache.nm1; j++) {
                sum += cache.mass_matrix(i - 1, j - 1) * unew(j + alg.m1 - 1);
            }
            ak(i + alg.m1 - 1) += sum;
        }
    }


    if (cache.job == 1) {
        sol(integrator.n, cache.e, ak, cache.ip);
        return 0;
    } else if (cache.job == 11 || cache.job == 13 || cache.job == 15) {
        int mm = alg.m1 / alg.m2;
        for (int j = 1; j <= alg.m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * alg.m2;
                sum = (ak(jkm - 1) + sum) / cache.fac;
                for (int i = 1; i <= cache.nm1; i++) {
                    int im1 = i + alg.m1;
                    ak(im1 - 1) = ak(im1 - 1) + cache.dfdu(i - 1, jkm - 1) * sum;
                }
            }
        }
        Vector<double> ak_seg = ak.segment(alg.m1, ak.size());
        sol(cache.nm1, cache.e, ak_seg, cache.ip);
        ak.segment(alg.m1, ak.size()) = ak_seg;
        for (int i = alg.m1; i >= 1; i--) {
            ak(i - 1) = (ak(i - 1) + ak(alg.m2 + i - 1)) / cache.fac;
        }
        return 0;
    } else if (cache.job == 2) {
        solb(integrator.n, cache.e, cache.mle, cache.mue, ak, cache.ip);
        return 0;
    } else if (cache.job == 12 || cache.job == 14) {
        int mm = alg.m1 / alg.m2;
        for (int j = 1; j <= alg.m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * alg.m2;
                sum = (ak(jkm - 1) + sum) / cache.fac;
                for (int i = std::max(1, j - alg.mujac); i <= std::min(cache.nm1, j + alg.mljac); i++) {
                    int im1 = i + alg.m1;
                    ak(im1 - 1) = ak(im1 - 1) + cache.dfdu(i + alg.mujac + 1 - j - 1, jkm - 1) * sum;
                }
            }
        }
        Vector<double> ak_seg = ak.segment(alg.m1, ak.size());
        solb(cache.nm1, cache.e, cache.mle, cache.mue, ak_seg, cache.ip);
        ak.segment(alg.m1, ak.size()) = ak_seg;
        for (int i = alg.m1; i >= 1; i--) {
            ak(i - 1) = (ak(i - 1) + ak(alg.m2 + i - 1)) / cache.fac;
        }
        return 0;
    } else if (cache.job == 3) {
        if (stage1) {
            for (int i = 1; i <= integrator.n; i++) {
                double sum = 0.0;
                for (int j = std::max(1, i - alg.mlmas); j <= std::min(integrator.n, i + alg.mumas); j++) {
                    sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j - 1) * unew(j - 1);
                }
                ak(i - 1) += sum;
            }
        }
        sol(integrator.n, cache.e, ak, cache.ip);
        return 0;
    } else if (cache.job == 4) {
        if (stage1) {
            for (int i = 1; i <= integrator.n; i++) {
                double sum = 0.0;
                for (int j = std::max(1, i - alg.mlmas); j <= std::min(integrator.n, i + alg.mumas); j++) {
                    sum += cache.mass_matrix(i - j + cache.mbdiag - 1, j - 1) * unew(j - 1);
                }
                ak(i + alg.m1 - 1) += sum;
            }
        }
        solb(integrator.n, cache.e, cache.mle, cache.mue, ak, cache.ip);
        return 0;
    } else if (cache.job == 5) {
        if (stage1) {
            for (int i = 1; i <= integrator.n; i++) {
                double sum = 0.0;
                for (int j = 1; j <= integrator.n; j++) {
                    sum += cache.mass_matrix(i - 1, j - 1) * unew(j - 1);
                }
                ak(i + alg.m1 - 1) += sum;
            }
        }
        sol(integrator.n, cache.e, ak, cache.ip);
    } else if (cache.job == 6) {
        if (stage1) {
            for (int i = 1; i <= integrator.n; i++) {
                double sum = 0.0;
                for (int j = 1; j <= integrator.n; j++) {
                    sum += cache.mass_matrix(i - 1, j - 1) * unew(j - 1);
                }
                ak(i - 1) += sum;
            }
            solb(integrator.n, cache.e, cache.mle, cache.mue, ak, cache.ip);
        }
        return 0;
    }
    return 0;
}

/**
 * Dense output for Rodas algorithm
 * @param integrator
 * @param alg
 * @param cache
 * @param i
 * @param t
 * @param dt
 * @return
 */
static double dense_output(ODEIntegrator &integrator, Rodas &alg, RodasCache &cache, int i, double t, double dt) {
    double s = (t - integrator.tprev) / dt;
    return cache.cont(i) * (1 - s) + s * (cache.cont(i + integrator.n) + (1 - s) * (cache.cont(i + integrator.n * 2)
            + s * cache.cont(i + integrator.n * 3)));
}


/**
 * Core integration for Rodas algorithm
 * @param integrator
 * @param alg
 * @param cache
 * @return
 */
static ODESolution rodas_core(ODEIntegrator &integrator, Rodas &alg, RodasCache &cache) {
    int ier = 0;

    init(integrator, alg, cache);

    ODESolution solution{};
    solution.ts.push_back(integrator.t);
    solution.us.push_back(integrator.u);

    // Basic integration step
    while (true) {
        if (integrator.num_steps > integrator.opts.max_num_steps) {
            solution.retcode = MaxIters;
            return solution;
        }
        if (0.1 * std::abs(integrator.dt) <= std::abs(integrator.t) * std::numeric_limits<double>::epsilon()) {
            solution.retcode = DtLessThanMin;
            return solution;
        }
        if (cache.last) {
            integrator.dt = cache.dtopt;
            solution.retcode = Success;
            return solution;
        }
        cache.dtopt = integrator.dt;
        if ((integrator.t + integrator.dt * 1.0001 - integrator.t_final) * integrator.tdir >= 0.0) {
            integrator.dt = integrator.t_final - integrator.t;
            cache.last = true;
        }
        // Compute jacobians
        integrator.f->dudt(cache.du1, integrator.u, integrator.t);
        integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
        if (!alg.autonomous) {
            integrator.f->dfdt(cache.dfdt, integrator.u, integrator.t);
        }

        while (true) {
            // Compute the stages
            cache.fac = 1.0 / (integrator.dt * alg.gamma);
            ier = decomp_real(integrator, alg, cache);
            if (ier != 0) {
                cache.num_singular++;
                if (cache.num_singular >= 5) {
                    //throw std::runtime_error("Rodas: Singular matrix.");
                    solution.retcode = SingularMatrix;
                    return solution;
                }
                integrator.dt *= 0.5;
                cache.reject = true;
                cache.last = false;
                continue;
            }
            integrator.num_decompositions++;

            // Prepare for the computation of the 6 stages
            double hc21 = alg.c21 / integrator.dt;
            double hc31 = alg.c31 / integrator.dt;
            double hc32 = alg.c32 / integrator.dt;
            double hc41 = alg.c41 / integrator.dt;
            double hc42 = alg.c42 / integrator.dt;
            double hc43 = alg.c43 / integrator.dt;
            double hc51 = alg.c51 / integrator.dt;
            double hc52 = alg.c52 / integrator.dt;
            double hc53 = alg.c53 / integrator.dt;
            double hc54 = alg.c54 / integrator.dt;
            double hc61 = alg.c61 / integrator.dt;
            double hc62 = alg.c62 / integrator.dt;
            double hc63 = alg.c63 / integrator.dt;
            double hc64 = alg.c64 / integrator.dt;
            double hc65 = alg.c65 / integrator.dt;

            if (!alg.autonomous) {
                cache.hd1 = integrator.dt * alg.d1;
                cache.hd2 = integrator.dt * alg.d2;
                cache.hd3 = integrator.dt * alg.d3;
                cache.hd4 = integrator.dt * alg.d4;
            }

            // The Stages

            //Solve for k1
            linearSolve(integrator, alg, cache, cache.du1, cache.ak1, cache.unew, cache.hd1, false);
            //Solve for k2
            cache.unew = integrator.u + alg.a21 * cache.ak1;
            integrator.f->dudt(cache.du, cache.unew, integrator.t + alg.c2 * integrator.dt);
            cache.unew = hc21 * cache.ak1;
            linearSolve(integrator, alg, cache, cache.du, cache.ak2, cache.unew, cache.hd2, true);
            //Solve for k3
            cache.unew = integrator.u + alg.a31 * cache.ak1 + alg.a32 * cache.ak2;
            integrator.f->dudt(cache.du, cache.unew, integrator.t + alg.c3 * integrator.dt);
            cache.unew = hc31 * cache.ak1 + hc32 * cache.ak2;
            linearSolve(integrator, alg, cache, cache.du, cache.ak3, cache.unew, cache.hd3, true);
            //Solve for k4
            cache.unew = integrator.u + alg.a41 * cache.ak1 + alg.a42 * cache.ak2 + alg.a43 * cache.ak3;
            integrator.f->dudt(cache.du, cache.unew, integrator.t + alg.c4 * integrator.dt);
            cache.unew = hc41 * cache.ak1 + hc42 * cache.ak2 + hc43 * cache.ak3;
            linearSolve(integrator, alg, cache, cache.du, cache.ak4, cache.unew, cache.hd4, true);
            //Solve for k5
            cache.unew = integrator.u + alg.a51 * cache.ak1 + alg.a52 * cache.ak2 + alg.a53 * cache.ak3 +
                    alg.a54 * cache.ak4;
            integrator.f->dudt(cache.du, cache.unew, integrator.t + integrator.dt);
            cache.ak6 = hc52 * cache.ak2 + hc54 * cache.ak4 + hc51 * cache.ak1 + hc53 * cache.ak3;
            linearSolve(integrator, alg, cache, cache.du, cache.ak5, cache.ak6, 0.0, true);

            // ------------ EMBEDDED SOLUTION ---------------
            cache.unew += cache.ak5;
            integrator.f->dudt(cache.du, cache.unew, integrator.t + integrator.dt);
            for (int i = 0; i < integrator.n; i++) {
                cache.cont[i] = hc61 * cache.ak1[i] + hc62 * cache.ak2[i] + hc65 * cache.ak5[i]
                        + hc64 * cache.ak4[i] + hc63 * cache.ak3[i];
            }
            linearSolve(integrator, alg, cache, cache.du, cache.ak6, cache.cont, 0.0, true);
            // ------------ NEW SOLUTION ---------------
            cache.unew += cache.ak6;

            integrator.num_linear_solves += 6;
            integrator.num_function_evaluations += 5;

            if (integrator.opts.dense) {
                for (int i = 0; i < integrator.n; i++) {
                    cache.cont[i] = integrator.u[i];
                    cache.cont(i + 2 * integrator.n) =
                            alg.d21 * cache.ak1(i) + alg.d22 * cache.ak2(i) + alg.d23 * cache.ak3(i) +
                                    alg.d24 * cache.ak4(i)
                                            + alg.d25 * cache.ak5(i);
                    cache.cont(i + 3 * integrator.n) =
                            alg.d31 * cache.ak1(i) + alg.d32 * cache.ak2(i) + alg.d33 * cache.ak3(i) +
                                    alg.d34 * cache.ak4(i)
                                            + alg.d35 * cache.ak5(i);
                }
            }

            /*
             *  Error Estimation
             */
            integrator.num_steps++;
            // Compute error estimation
            double err = 0.0, sk;
            for (int i = 0; i < integrator.n; i++) {
                sk = integrator.opts.abstol +
                        integrator.opts.reltol * std::max(std::abs(integrator.u[i]), std::abs(cache.unew[i]));
                err += pow(cache.ak6[i] / sk, 2);
            }
            err = std::sqrt(err / double(integrator.n));

            // Computation of new dt; we require 0.2 <= dtNew / dt <= 6
            cache.fac = std::max(alg.fac2, std::min(alg.fac1, pow(err, 0.25) / alg.safe));
            double dtNew = integrator.dt / cache.fac;

            /*
             * Is the error small enough?
             */
            if (err <= 1.0) {
                // Step is accepted
                integrator.num_accept++;
                if (cache.pred) {
                    // Use the predictive controller of Gustafsson
                    if (integrator.num_accept > 1) {
                        double facgus = (cache.dtacc / integrator.dt) * pow(err * err / cache.erracc, 0.25) / alg.safe;
                        facgus = std::max(alg.fac2, std::min(alg.fac1, facgus));
                        cache.fac = std::max(cache.fac, facgus);
                        dtNew = integrator.dt / cache.fac;
                    }
                    cache.dtacc = integrator.dt;
                    cache.erracc = std::max(1e-2, err);
                }
                integrator.u = cache.unew;
                integrator.tprev = integrator.t;
                integrator.t += integrator.dt;

                if (integrator.opts.dense) {
                    for (int i = 0; i <= integrator.n; i++) {
                        cache.cont(i + integrator.n) = integrator.u(i);
                    }
                }
                solution.ts.push_back(integrator.t);
                solution.us.push_back(integrator.u);

                if (std::abs(dtNew) > integrator.opts.dtmax) {
                    dtNew = integrator.tdir * integrator.opts.dtmax;
                }
                if (cache.reject) {
                    dtNew = integrator.tdir * std::min(std::abs(dtNew), std::abs(integrator.dt));
                }
                cache.reject = false;
                integrator.dt = dtNew;
                break;
            } else {
                // Step is rejected
                cache.reject = true;
                cache.last = false;
                integrator.dt = dtNew;
                if (integrator.num_accept >= 1) {
                    integrator.num_reject++;
                }
                continue;
            }
        }
    }
}

/**
 * Solve the ODE problem using the embedded Rosenbrock method of order (3)4
 * with step size control
 * @param prob Problem to solve
 * @param alg Algorithm to use (Rodas in this case)
 * @param opts Integration options
 * @return Solution
 */
ODESolution solve(ODEProblem &prob, Rodas &alg, ODEIntegratorOptions &opts) {
    ODEIntegrator integrator{prob, opts};
    RodasCache cache{};
    return rodas_core(integrator, alg, cache);
}

/**
 * Solve the ODE problem using the embedded Rosenbrock method of order (3)4
 * with step size control. Uses default options.
 * @param prob Problem to solve
 * @param alg Algorithm to use (Rodas in this case)
 * @return Solution
 */
ODESolution solve(ODEProblem &prob, Rodas &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}

}
}

#endif //LANRE_DIFFEQ_RODAS_HPP
