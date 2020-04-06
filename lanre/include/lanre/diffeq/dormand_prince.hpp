//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_DORMAND_PRINCE_HPP
#define LANRE_DIFFEQ_DORMAND_PRINCE_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/algorithm.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/solution.hpp"
#include "lanre/diffeq/integrator.hpp"
#include <exception>

namespace lanre {
namespace diffeq {

class DormandPrince5Stiff : public std::exception {
private:
    const char *msg;
public:

    DormandPrince5Stiff(const double t_t) : msg(
            ("DormandPrince5: The problem seems to become stiff at x " + std::to_string(t_t)).c_str()) {}

    [[nodiscard]] const char *what() const noexcept override {
        return msg;
    }
};

template<class DEIntegrator, class AlgorithmCache>
void initialize_dt(DEIntegrator &integrator, AlgorithmCache &cache) {
    double dnf = 0.0;
    double dny = 0.0;
    for (int i = 0; i < integrator.n; i++) {
        double sk = integrator.opts.abstol + integrator.opts.reltol * fabs(integrator.u[i]);
        double sqr = cache.du[i] / sk;
        dnf += sqr * sqr;
        sqr = integrator.u[i] / sk;
        dny += sqr * sqr;
    }


    if ((dnf <= 1.0E-10) || (dny <= 1.0E-10)) {
        integrator.dt = 1.0E-6;
    } else {
        integrator.dt = sqrt(dny / dnf) * 0.01;
    }

    integrator.dt = std::min(integrator.dt, integrator.opts.dtmax);
    integrator.dt = std::copysign(integrator.dt, integrator.tdir);

    // perform an explicit Euler step
    for (int i = 0; i < integrator.n; i++) {
        cache.unew[i] = integrator.u[i] + integrator.dt * cache.du[i];
    }
    integrator.f->dudt(cache.dunew, cache.unew, integrator.t + integrator.dt);

    // estimate the second derivative of the solution
    double der2 = 0.0;
    for (int i = 0; i < integrator.n; i++) {
        double sk = integrator.opts.abstol + integrator.opts.reltol * fabs(integrator.u[i]);
        double sqr = (cache.dunew[i] - cache.du[i]) / sk;
        der2 += sqr * sqr;
    }
    der2 = sqrt(der2) / integrator.dt;

    // step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01
    double der12 = std::max(fabs(der2), sqrt(dnf));
    double dtnew;
    if (der12 <= 1.0E-15) {
        dtnew = std::max(1.0E-6, fabs(integrator.dt) * 1.0E-3);
    } else {
        dtnew = pow(0.01 / der12, 1.0 / 5.0);
    }
    integrator.dt = std::min(100.0 * integrator.dt, std::min(dtnew, integrator.opts.dtmax));

    integrator.dt = std::copysign(integrator.dt, integrator.tdir);
}

/**
 * 7-Stage Dormand-Prince with 4th-order global accurracy of embedding and 5rd-order global
 * accurary of the method.
 *
 */

struct DormandPrince5Cache : public ODEAlgorithmCache {
    // Variables for stiffness detection
    int max_num_stiff;
    int num_stiff = 0;
    int num_non_stiff = 0;
    double dt_lam_b = 0;
    double next_prev_step_ratio;
    // Vectors for computing y1 = y0 + f(x0+c*h, y0+h*b1*k1+...+bs*ks)
    Vector<double> k2, k3, k4, k5, k6;
    // Vectors for continuous output
    Vector<double> rcont1, rcont2, rcont3, rcont4, rcont5;
    // States variables
    Vector<double> unew, du, dunew, uerr, ustiff;
    double dt_new;       // store proposed value of dt
    bool reject = false; // was step rejected?
    bool last = false;   // last step?

    DormandPrince5Cache(int n)
            : max_num_stiff(1000), num_stiff(0), num_non_stiff(0),
              dt_lam_b(0.0), next_prev_step_ratio(0.0),
              k2(Vector<double>{n}),
              k3(Vector<double>{n}), k4(Vector<double>{n}),
              k5(Vector<double>{n}), k6(Vector<double>{n}),
              rcont1(Vector<double>{n}), rcont2(Vector<double>{n}),
              rcont3(Vector<double>{n}), rcont4(Vector<double>{n}),
              rcont5(Vector<double>{n}),
              unew(Vector<double>{n}),
              du(Vector<double>{n}),
              dunew(Vector<double>{n}),
              uerr(Vector<double>{n}),
              ustiff(Vector<double>{n}),
              dt_new(0.0), reject(false), last(false) {}
};

struct DormandPrince5 : public ODEAlgorithm {
    // Butcher c_i's
    const double c2 = 0.2;
    const double c3 = 0.3;
    const double c4 = 0.8;
    const double c5 = 8.0 / 9.0;
    // Butcher a_ij's
    const double a21 = 0.2;
    const double a31 = 3.0 / 40.0;
    const double a32 = 9.0 / 40.0;
    const double a41 = 44.0 / 45.0;
    const double a42 = -56.0 / 15.0;
    const double a43 = 32.0 / 9.0;
    const double a51 = 19372.0 / 6561.0;
    const double a52 = -25360.0 / 2187.0;
    const double a53 = 64448.0 / 6561.0;
    const double a54 = -212.0 / 729.0;
    const double a61 = 9017.0 / 3168.0;
    const double a62 = -355.0 / 33.0;
    const double a63 = 46732.0 / 5247.0;
    const double a64 = 49.0 / 176.0;
    const double a65 = -5103.0 / 18656.0;
    // Butcher b_i's: y1 = y0 + h * (b1 * k1 + ... + bs * ks)
    const double a71 = 35.0 / 384.0;
    const double a73 = 500.0 / 1113.0;
    const double a74 = 125.0 / 192.0;
    const double a75 = -2187.0 / 6784.0;
    const double a76 = 11.0 / 84.0;
    // Error estimation constants: b_i^* - b_i
    const double e1 = 71.0 / 57600.0;
    const double e3 = -71.0 / 16695.0;
    const double e4 = 71.0 / 1920.0;
    const double e5 = -17253.0 / 339200.0;
    const double e6 = 22.0 / 525.0;
    const double e7 = -1.0 / 40.0;
    // Continuous output parameters
    const double d1 = -12715105075.0 / 11282082432.0;
    const double d3 = 87487479700.0 / 32700410799.0;
    const double d4 = -10690763975.0 / 1880347072.0;
    const double d5 = 701980252875.0 / 199316789632.0;
    const double d6 = -1453857185.0 / 822651844.0;
    const double d7 = 69997945.0 / 29380423.0;
    // Beta for stabilized step size control
    const double beta = 0.04;
    const double alpha = 0.2 - 0.0 * 0.75;
    // safety factor for stepsize prediction
    const double safe = 0.9;
    // Parameters for stepsize prediction
    const double minNextPrevStepRatio = 1.0 / 0.2;  // < dtNext / dtPrev
    const double maxNextPrevStepRatio = 1.0 / 10.0; // > dtNext / dtPrev
};

/**
 * Compute the 7 Runge-Kutta stages for 4/5 order Dormand-Prince
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 */
static void compute_stages(ODEIntegrator &integrator, DormandPrince5 &alg, DormandPrince5Cache &cache) {
    // Computation of the 6 stages
    cache.unew = integrator.u + integrator.dt * alg.a21 * cache.du;
    integrator.f->dudt(cache.k2, cache.unew, integrator.t + alg.c2 * integrator.dt);
    cache.unew = integrator.u + integrator.dt * (alg.a31 * cache.du + alg.a32 * cache.k2);
    integrator.f->dudt(cache.k3, cache.unew, integrator.t + alg.c3 * integrator.dt);
    cache.unew = integrator.u + integrator.dt * (alg.a41 * cache.du + alg.a42 * cache.k2 + alg.a43 * cache.k3);
    integrator.f->dudt(cache.k4, cache.unew, integrator.t + alg.c4 * integrator.dt);
    cache.unew = integrator.u +
            integrator.dt * (alg.a51 * cache.du + alg.a52 * cache.k2 + alg.a53 * cache.k3 + alg.a54 * cache.k4);
    integrator.f->dudt(cache.k5, cache.unew, integrator.t + alg.c5 * integrator.dt);
    cache.ustiff = integrator.u + integrator.dt *
            (alg.a61 * cache.du + alg.a62 * cache.k2 + alg.a63 * cache.k3 + alg.a64 * cache.k4 + alg.a65 * cache.k5);
    // The RK step
    double tph = integrator.t + integrator.dt;
    integrator.f->dudt(cache.k6, cache.ustiff, tph);
    cache.unew = integrator.u + integrator.dt *
            (alg.a71 * cache.du + alg.a73 * cache.k3 + alg.a74 * cache.k4 + alg.a75 * cache.k5 + alg.a76 * cache.k6);
    // Error estimate using embedded solution
    integrator.f->dudt(cache.dunew, cache.unew, tph);
    cache.uerr = integrator.dt * (alg.e1 * cache.du + alg.e3 * cache.k3 + alg.e4 * cache.k4 + alg.e5 * cache.k5 + alg
            .e6 * cache.k6 + alg.e7 * cache.dunew);
    // Update counting variables
    integrator.num_function_evaluations += 6;
}

/**
 * Prepare the continuous vectors for dense output.
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 */
static void prepare_dense(ODEIntegrator &integrator, DormandPrince5 &alg, DormandPrince5Cache &cache) {
    cache.rcont1 = integrator.u;
    cache.rcont2 = cache.unew - integrator.u;
    cache.rcont3 = integrator.dt * cache.du - cache.rcont2;
    cache.rcont4 = cache.rcont2 - integrator.dt * cache.dunew - cache.rcont3;
    cache.rcont5 = integrator.dt * (alg.d1 * cache.du +
            alg.d3 * cache.k3 +
            alg.d4 * cache.k4 +
            alg.d5 * cache.k5 +
            alg.d6 * cache.k6 +
            alg.d7 * cache.dunew);
}

/**
 * Compute the dense output for the ith component of the solution between times
 * t and t + dt
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 * @param i Component of solution to compute dense output for
 * @param t Current time value
 * @param dt Current step value
 * @return Solution value at time t
 */
static double dense_output(
        ODEIntegrator &integrator,
        DormandPrince5 &alg,
        DormandPrince5Cache &cache,
        const int i,
        const double t,
        const double dt
) {
    double s = (t - integrator.tprev) / dt;
    double s1 = 1.0 - s;
    return cache.rcont1[i] +
            s * (cache.rcont2[i] + s1 * (cache.rcont3[i] + s * (cache.rcont4[i] + s1 * cache.rcont5[i])));
}

/**
 * Compute the estimated error from the previousl performed step.
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 * @return Estimated error
 */
static double error(
        ODEIntegrator &integrator,
        DormandPrince5 &alg,
        DormandPrince5Cache &cache
) {
    double err = 0.0;
    for (int i = 0; i < integrator.n; i++) {
        double sk = integrator.opts.abstol + integrator.opts.reltol * std::max(abs(integrator.u[i]), abs(cache
                                                                                                                 .unew[i]));
        double sqr = cache.uerr[i] / sk;
        err += sqr * sqr;
    }
    return sqrt(err / double(integrator.n));
}

/**
 * Prepare the integrator and cache for the next step.
 * @param t_err Estimated error from step just taken
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 */
static void prepare_next_step(
        const double t_err,
        ODEIntegrator &integrator,
        DormandPrince5 &alg,
        DormandPrince5Cache &cache
) {
    // computation of hnew
    double fac11 = pow(t_err, alg.alpha);
    // Lund-stabilization
    double fac = fac11 / pow(cache.next_prev_step_ratio, alg.beta);
    // we require minNextPrevStepRatio <= hnew/h <= m_maxNextPrevStepRatio
    fac = std::max(alg.maxNextPrevStepRatio, std::min(alg.minNextPrevStepRatio, fac / alg.safe));
    double dtnew = integrator.dt / fac;

    if (t_err <= 1.0) {
        /* step accepted */

        cache.next_prev_step_ratio = std::max(t_err, 1.0E-4);
        integrator.num_accept++;

        /* stiffness detection */
        if (!(integrator.num_accept % cache.max_num_stiff) || (cache.num_stiff > 0)) {
            double stnum = 0.0;
            double stden = 0.0;
            for (unsigned int i = 0; i < integrator.n; i++) {
                double sqr = cache.k2[i] - cache.k6[i];
                stnum += sqr * sqr;
                sqr = cache.unew[i] - cache.ustiff[i];
                stden += sqr * sqr;
            }
            if (stden > 0.0)
                cache.dt_lam_b = integrator.dt * sqrt(stnum / stden);
            if (cache.dt_lam_b > 3.25) {
                cache.num_non_stiff = 0;
                cache.num_stiff++;
                if (cache.num_stiff == 15) {
                    throw DormandPrince5Stiff(integrator.t);
                }
            } else {
                cache.num_non_stiff++;
                if (cache.num_non_stiff == 6) {
                    cache.num_stiff = 0;
                }
            }
        }
        if (integrator.opts.dense) {
            prepare_dense(integrator, alg, cache);
        }

        cache.du = cache.dunew;
        integrator.u = cache.unew;
        integrator.tprev = integrator.t;
        integrator.t += integrator.dt;

        if (fabs(dtnew) > integrator.opts.dtmax) {
            dtnew = integrator.tdir * integrator.opts.dtmax;
        }
        if (cache.reject) {
            dtnew = integrator.tdir * std::min(fabs(dtnew), fabs(integrator.dt));
        }

        cache.reject = false;
    } else {
        /* step rejected */
        dtnew = integrator.dt / std::min(alg.minNextPrevStepRatio, fac11 / alg.safe);
        cache.reject = true;
        if (integrator.num_accept >= 1) {
            integrator.num_reject++;
        }
        cache.last = false;
    }

    integrator.dtprev = integrator.dt;
    integrator.dt = dtnew;
}

/**
 * Perform a single step using DormandPrince5
 * @param t_err Estimated error from step just taken
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 */
void step(
        ODEIntegrator &integrator,
        DormandPrince5 &alg,
        DormandPrince5Cache &cache
) {
    while (true) {
        compute_stages(integrator, alg, cache);
        double err = error(integrator, alg, cache);
        prepare_next_step(err, integrator, alg, cache);
        if (!cache.reject) {
            break;
        }
        if (abs(integrator.dt) <= abs(integrator.t) * std::numeric_limits<double>::epsilon())
            throw std::runtime_error("Step size underflow in Dormand-Prince-5");
    }
    integrator.num_steps++;
}

/**
 * Solve the ODE problem using the Dormand-Prince 5th order algorithm
 * @param prob ODE problem to solve
 * @param alg DormandPrince5 algorithm
 * @param opts Integration options
 * @return Solution
 */
ODESolution solve(
        ODEProblem &prob,
        DormandPrince5 &alg,
        ODEIntegratorOptions &opts
) {

    ODEIntegrator integrator{prob, opts};
    DormandPrince5Cache cache{integrator.n};

    integrator.sol.ts.push_back(integrator.t);
    integrator.sol.us.push_back(integrator.u);

    // Compute du
    integrator.f->dudt(cache.du, integrator.u, integrator.t);

    // Compute solution by stepping until we reach end of time interval.
    while (!cache.last) {

        if ((integrator.t + 1.01 * integrator.dt - prob.t_span.second) * integrator.tdir > 0.0) {
            integrator.dt = prob.t_span.second - integrator.t;
            cache.last = true;
        }

        step(integrator, alg, cache);
        integrator.sol.ts.push_back(integrator.t);
        integrator.sol.us.push_back(integrator.u);

        if (integrator.num_steps > integrator.opts.max_num_steps) {
            throw std::runtime_error("too many steps..");
        }
    }
    return integrator.sol;
}

/**
 * Solve the ODE problem using the Dormand-Prince 5th order algorithm
 * @param prob ODE problem to solve
 * @param alg DormandPrince5 algorithm
 * @return Solution
 */
ODESolution solve(ODEProblem &prob, DormandPrince5 &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}


/**
 * 8th-order Dormand-Prince with 5th-order + 3rd-order estimators yielding
 * 7th-order dense output.
 *
 */

struct DormandPrince8Cache : public ODEAlgorithmCache {
    // Variables for stiffness detection
    int max_num_stiff;
    int num_stiff = 0;
    int num_non_stiff = 0;
    double dt_lam_b = 0;
    double next_prev_step_ratio;
    // Vectors for computing y1 = y0 + f(x0+c*h, y0+h*b1*k1+...+bs*ks)
    Vector<double> k2, k3, k4, k5, k6, k7, k8, k9, k10;
    // Vectors for continuous output
    Vector<double> rcont1, rcont2, rcont3, rcont4, rcont5, rcont6, rcont7, rcont8;
    // States variables
    Vector<double> unew, du, dunew, uerr, ustiff;
    double dt_new;       // store proposed value of dt
    bool reject = false; // was step rejected?
    bool last = false;   // last step?

    DormandPrince8Cache(int n)
            : max_num_stiff(1000), num_stiff(0), num_non_stiff(0),
              dt_lam_b(0.0), next_prev_step_ratio(0.0),
              k2(Vector<double>{n}),
              k3(Vector<double>{n}), k4(Vector<double>{n}),
              k5(Vector<double>{n}), k6(Vector<double>{n}),
              k7(Vector<double>{n}), k8(Vector<double>{n}),
              k9(Vector<double>{n}), k10(Vector<double>{n}),
              rcont1(Vector<double>{n}), rcont2(Vector<double>{n}),
              rcont3(Vector<double>{n}), rcont4(Vector<double>{n}),
              rcont5(Vector<double>{n}), rcont6(Vector<double>{n}),
              rcont7(Vector<double>{n}), rcont8(Vector<double>{n}),
              unew(Vector<double>{n}),
              du(Vector<double>{n}),
              dunew(Vector<double>{n}),
              uerr(Vector<double>{n}),
              ustiff(Vector<double>{n}),
              dt_new(0.0), reject(false), last(false) {}
};

struct DormandPrince8 : public ODEAlgorithm {
    // Butcher c_i's
    const double c2 = 0.526001519587677318785587544488E-01;
    const double c3 = 0.789002279381515978178381316732E-01;
    const double c4 = 0.118350341907227396726757197510E+00;
    const double c5 = 0.281649658092772603273242802490E+00;
    const double c6 = 0.333333333333333333333333333333E+00;
    const double c7 = 0.25E+00;
    const double c8 = 0.307692307692307692307692307692E+00;
    const double c9 = 0.651282051282051282051282051282E+00;
    const double c10 = 0.6E+00;
    const double c11 = 0.857142857142857142857142857142E+00;
    const double c14 = 0.1E+00;
    const double c15 = 0.2E+00;
    const double c16 = 0.777777777777777777777777777778E+00;
    // Butcher b_i's
    const double b1 = 5.42937341165687622380535766363E-2;
    const double b6 = 4.45031289275240888144113950566E0;
    const double b7 = 1.89151789931450038304281599044E0;
    const double b8 = -5.8012039600105847814672114227E0;
    const double b9 = 3.1116436695781989440891606237E-1;
    const double b10 = -1.52160949662516078556178806805E-1;
    const double b11 = 2.01365400804030348374776537501E-1;
    const double b12 = 4.47106157277725905176885569043E-2;
    const double bhh1 = 0.244094488188976377952755905512E+00;
    const double bhh2 = 0.733846688281611857341361741547E+00;
    const double bhh3 = 0.220588235294117647058823529412E-01;
    // Butcher a_ij's
    const double a21 = 5.26001519587677318785587544488E-2;
    const double a31 = 1.97250569845378994544595329183E-2;
    const double a32 = 5.91751709536136983633785987549E-2;
    const double a41 = 2.95875854768068491816892993775E-2;
    const double a43 = 8.87627564304205475450678981324E-2;
    const double a51 = 2.41365134159266685502369798665E-1;
    const double a53 = -8.84549479328286085344864962717E-1;
    const double a54 = 9.24834003261792003115737966543E-1;
    const double a61 = 3.7037037037037037037037037037E-2;
    const double a64 = 1.70828608729473871279604482173E-1;
    const double a65 = 1.25467687566822425016691814123E-1;
    const double a71 = 3.7109375E-2;
    const double a74 = 1.70252211019544039314978060272E-1;
    const double a75 = 6.02165389804559606850219397283E-2;
    const double a76 = -1.7578125E-2;
    const double a81 = 3.70920001185047927108779319836E-2;
    const double a84 = 1.70383925712239993810214054705E-1;
    const double a85 = 1.07262030446373284651809199168E-1;
    const double a86 = -1.53194377486244017527936158236E-2;
    const double a87 = 8.27378916381402288758473766002E-3;
    const double a91 = 6.24110958716075717114429577812E-1;
    const double a94 = -3.36089262944694129406857109825E0;
    const double a95 = -8.68219346841726006818189891453E-1;
    const double a96 = 2.75920996994467083049415600797E1;
    const double a97 = 2.01540675504778934086186788979E1;
    const double a98 = -4.34898841810699588477366255144E1;
    const double a101 = 4.77662536438264365890433908527E-1;
    const double a104 = -2.48811461997166764192642586468E0;
    const double a105 = -5.90290826836842996371446475743E-1;
    const double a106 = 2.12300514481811942347288949897E1;
    const double a107 = 1.52792336328824235832596922938E1;
    const double a108 = -3.32882109689848629194453265587E1;
    const double a109 = -2.03312017085086261358222928593E-2;
    const double a111 = -9.3714243008598732571704021658E-1;
    const double a114 = 5.18637242884406370830023853209E0;
    const double a115 = 1.09143734899672957818500254654E0;
    const double a116 = -8.14978701074692612513997267357E0;
    const double a117 = -1.85200656599969598641566180701E1;
    const double a118 = 2.27394870993505042818970056734E1;
    const double a119 = 2.49360555267965238987089396762E0;
    const double a1110 = -3.0467644718982195003823669022E0;
    const double a121 = 2.27331014751653820792359768449E0;
    const double a124 = -1.05344954667372501984066689879E1;
    const double a125 = -2.00087205822486249909675718444E0;
    const double a126 = -1.79589318631187989172765950534E1;
    const double a127 = 2.79488845294199600508499808837E1;
    const double a128 = -2.85899827713502369474065508674E0;
    const double a129 = -8.87285693353062954433549289258E0;
    const double a1210 = 1.23605671757943030647266201528E1;
    const double a1211 = 6.43392746015763530355970484046E-1;
    const double a141 = 5.61675022830479523392909219681E-2;
    const double a147 = 2.53500210216624811088794765333E-1;
    const double a148 = -2.46239037470802489917441475441E-1;
    const double a149 = -1.24191423263816360469010140626E-1;
    const double a1410 = 1.5329179827876569731206322685E-1;
    const double a1411 = 8.20105229563468988491666602057E-3;
    const double a1412 = 7.56789766054569976138603589584E-3;
    const double a1413 = -8.298E-3;
    const double a151 = 3.18346481635021405060768473261E-2;
    const double a156 = 2.83009096723667755288322961402E-2;
    const double a157 = 5.35419883074385676223797384372E-2;
    const double a158 = -5.49237485713909884646569340306E-2;
    const double a1511 = -1.08347328697249322858509316994E-4;
    const double a1512 = 3.82571090835658412954920192323E-4;
    const double a1513 = -3.40465008687404560802977114492E-4;
    const double a1514 = 1.41312443674632500278074618366E-1;
    const double a161 = -4.28896301583791923408573538692E-1;
    const double a166 = -4.69762141536116384314449447206E0;
    const double a167 = 7.68342119606259904184240953878E0;
    const double a168 = 4.06898981839711007970213554331E0;
    const double a169 = 3.56727187455281109270669543021E-1;
    const double a1613 = -1.39902416515901462129418009734E-3;
    const double a1614 = 2.9475147891527723389556272149E0;
    const double a1615 = -9.15095847217987001081870187138E0;
    // Error estimation constants
    const double er1 = 0.1312004499419488073250102996E-01;
    const double er6 = -0.1225156446376204440720569753E+01;
    const double er7 = -0.4957589496572501915214079952E+00;
    const double er8 = 0.1664377182454986536961530415E+01;
    const double er9 = -0.3503288487499736816886487290E+00;
    const double er10 = 0.3341791187130174790297318841E+00;
    const double er11 = 0.8192320648511571246570742613E-01;
    const double er12 = -0.2235530786388629525884427845E-01;
    // Continuous output parameters
    const double d41 = -0.84289382761090128651353491142E+01;
    const double d46 = 0.56671495351937776962531783590E+00;
    const double d47 = -0.30689499459498916912797304727E+01;
    const double d48 = 0.23846676565120698287728149680E+01;
    const double d49 = 0.21170345824450282767155149946E+01;
    const double d410 = -0.87139158377797299206789907490E+00;
    const double d411 = 0.22404374302607882758541771650E+01;
    const double d412 = 0.63157877876946881815570249290E+00;
    const double d413 = -0.88990336451333310820698117400E-01;
    const double d414 = 0.18148505520854727256656404962E+02;
    const double d415 = -0.91946323924783554000451984436E+01;
    const double d416 = -0.44360363875948939664310572000E+01;
    const double d51 = 0.10427508642579134603413151009E+02;
    const double d56 = 0.24228349177525818288430175319E+03;
    const double d57 = 0.16520045171727028198505394887E+03;
    const double d58 = -0.37454675472269020279518312152E+03;
    const double d59 = -0.22113666853125306036270938578E+02;
    const double d510 = 0.77334326684722638389603898808E+01;
    const double d511 = -0.30674084731089398182061213626E+02;
    const double d512 = -0.93321305264302278729567221706E+01;
    const double d513 = 0.15697238121770843886131091075E+02;
    const double d514 = -0.31139403219565177677282850411E+02;
    const double d515 = -0.93529243588444783865713862664E+01;
    const double d516 = 0.35816841486394083752465898540E+02;
    const double d61 = 0.19985053242002433820987653617E+02;
    const double d66 = -0.38703730874935176555105901742E+03;
    const double d67 = -0.18917813819516756882830838328E+03;
    const double d68 = 0.52780815920542364900561016686E+03;
    const double d69 = -0.11573902539959630126141871134E+02;
    const double d610 = 0.68812326946963000169666922661E+01;
    const double d611 = -0.10006050966910838403183860980E+01;
    const double d612 = 0.77771377980534432092869265740E+00;
    const double d613 = -0.27782057523535084065932004339E+01;
    const double d614 = -0.60196695231264120758267380846E+02;
    const double d615 = 0.84320405506677161018159903784E+02;
    const double d616 = 0.11992291136182789328035130030E+02;
    const double d71 = -0.25693933462703749003312586129E+02;
    const double d76 = -0.15418974869023643374053993627E+03;
    const double d77 = -0.23152937917604549567536039109E+03;
    const double d78 = 0.35763911791061412378285349910E+03;
    const double d79 = 0.93405324183624310003907691704E+02;
    const double d710 = -0.37458323136451633156875139351E+02;
    const double d711 = 0.10409964950896230045147246184E+03;
    const double d712 = 0.29840293426660503123344363579E+02;
    const double d713 = -0.43533456590011143754432175058E+02;
    const double d714 = 0.96324553959188282948394950600E+02;
    const double d715 = -0.39177261675615439165231486172E+02;
    const double d716 = -0.14972683625798562581422125276E+03;
    // Beta for stabilized step size control
    const double beta = 0.0;
    const double alpha = 1.0 / 8.0 - 0.0 * 0.2;
    // safety factor for stepsize prediction
    const double safe = 0.9;
    // Parameters for stepsize prediction
    const double minNextPrevStepRatio = 1.0 / 0.333; // < dtNext / dtPrev
    const double maxNextPrevStepRatio = 1.0 / 6.0; // > dtNext / dtPrev
};


/**
 * Compute the 10 Runge-Kutta stages for 8/5 order Dormand-Prince
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince8 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince8
 */
static void compute_stages(ODEIntegrator &integrator, DormandPrince8 &alg, DormandPrince8Cache &cache) {
    cache.unew = integrator.u + integrator.dt * alg.a21 * cache.du;
    integrator.f->dudt(cache.k2, cache.unew, integrator.t + alg.c2 * integrator.dt);

    cache.unew = integrator.u + integrator.dt * (alg.a31 * cache.du + alg.a32 * cache.k2);
    integrator.f->dudt(cache.k3, cache.unew, integrator.t + alg.c3 * integrator.dt);

    cache.unew = integrator.u + integrator.dt * (alg.a41 * cache.du + alg.a43 * cache.k3);
    integrator.f->dudt(cache.k4, cache.unew, integrator.t + alg.c4 * integrator.dt);

    cache.unew = integrator.u + integrator.dt * (alg.a51 * cache.du + alg.a53 * cache.k3 + alg.a54 * cache.k4);
    integrator.f->dudt(cache.k5, cache.unew, integrator.t + alg.c5 * integrator.dt);

    cache.unew = integrator.u + integrator.dt * (alg.a61 * cache.du + alg.a64 * cache.k4 + alg.a65 * cache.k5);
    integrator.f->dudt(cache.k6, cache.unew, integrator.t + alg.c6 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a71 * cache.du + alg.a74 * cache.k4 + alg.a75 * cache.k5 + alg.a76 * cache.k6);
    integrator.f->dudt(cache.k7, cache.unew, integrator.t + alg.c7 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a81 * cache.du + alg.a84 * cache.k4 + alg.a85 * cache.k5 + alg.a86 * cache.k6 +
                    alg.a87 * cache.k7);
    integrator.f->dudt(cache.k8, cache.unew, integrator.t + alg.c8 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a91 * cache.du + alg.a94 * cache.k4 + alg.a95 * cache.k5 + alg.a96 * cache.k6 +
                    alg.a97 * cache.k7 + alg.a98 * cache.k8);
    integrator.f->dudt(cache.k9, cache.unew, integrator.t + alg.c9 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a101 * cache.du + alg.a104 * cache.k4 + alg.a105 * cache.k5 + alg.a106 * cache.k6 +
                    alg.a107 * cache.k7 + alg.a108 * cache.k8 + alg.a109 * cache.k9);
    integrator.f->dudt(cache.k10, cache.unew, integrator.t + alg.c10 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a111 * cache.du + alg.a114 * cache.k4 + alg.a115 * cache.k5 + alg.a116 * cache.k6 +
                    alg.a117 * cache.k7 + alg.a118 * cache.k8 + alg.a119 * cache.k9 + alg.a1110 * cache.k10);
    integrator.f->dudt(cache.k2, cache.unew, integrator.t + alg.c11 * integrator.dt);

    cache.unew = integrator.u +
            integrator.dt * (alg.a121 * cache.du + alg.a124 * cache.k4 + alg.a125 * cache.k5 + alg.a126 * cache.k6 +
                    alg.a127 * cache.k7 + alg.a128 * cache.k8 + alg.a129 * cache.k9 +
                    alg.a1210 * cache.k10 + alg.a1211 * cache.k2);
    integrator.f->dudt(cache.k3, cache.unew, integrator.t + integrator.dt);
    integrator.num_function_evaluations += 11;


    cache.k4 = alg.b1 * cache.du + alg.b6 * cache.k6 + alg.b7 * cache.k7 + alg.b8 * cache.k8 + alg.b9 * cache.k9 +
            alg.b10 * cache.k10 + alg.b11 * cache.k2 + alg.b12 * cache.k3;
    cache.k5 = integrator.u + integrator.dt * cache.k4;
}

/**
 * Compute the estimated error from the previousl performed step.
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince8 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince8
 * @return Estimated error
 */
static double error(ODEIntegrator &integrator, DormandPrince8 &alg, DormandPrince8Cache &cache) {
    double err = 0.0;
    double err2 = 0.0;
    for (int i = 0; i < integrator.n; i++) {
        double sk =
                integrator.opts.abstol + integrator.opts.reltol * std::max(fabs(integrator.u[i]), fabs(cache.k5[i]));
        double erri = cache.k4[i] - alg.bhh1 * cache.du[i] - alg.bhh2 * cache.k9[i] - alg.bhh3 * cache.k3[i];
        double sqr = erri / sk;
        err2 += sqr * sqr;
        erri = alg.er1 * cache.du[i] + alg.er6 * cache.k6[i] + alg.er7 * cache.k7[i] + alg.er8 * cache.k8[i] +
                alg.er9 * cache.k9[i] + alg.er10 * cache.k10[i] + alg.er11 * cache.k2[i] + alg.er12 * cache.k3[i];
        sqr = erri / sk;
        err += sqr * sqr;
    }
    double deno = err + 0.01 * err2;
    if (deno <= 0.0) {
        deno = 1.0;
    }
    return fabs(integrator.dt) * err * sqrt(1.0 / (deno * double(integrator.n)));
}

/**
 * Prepare the continuous vectors for dense output.
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince8 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince8
 */
static void prepare_dense(ODEIntegrator &integrator, DormandPrince8 &alg, DormandPrince8Cache &cache) {
    double ydiff, bspl;
    for (int i = 0; i < integrator.n; i++) {
        cache.rcont1[i] = integrator.u[i];
        ydiff = cache.k5[i] - integrator.u[i];
        cache.rcont2[i] = ydiff;
        bspl = integrator.dt * cache.du[i] - ydiff;
        cache.rcont3[i] = bspl;
        cache.rcont4[i] = ydiff - integrator.dt * cache.k4[i] - bspl;
        cache.rcont5[i] = alg.d41 * cache.du[i] + alg.d46 * cache.k6[i] + alg.d47 * cache.k7[i] +
                alg.d48 * cache.k8[i] + alg.d49 * cache.k9[i] + alg.d410 * cache.k10[i] +
                alg.d411 * cache.k2[i] + alg.d412 * cache.k3[i];
        cache.rcont6[i] = alg.d51 * cache.du[i] + alg.d56 * cache.k6[i] + alg.d57 * cache.k7[i] +
                alg.d58 * cache.k8[i] + alg.d59 * cache.k9[i] + alg.d510 * cache.k10[i] +
                alg.d511 * cache.k2[i] + alg.d512 * cache.k3[i];
        cache.rcont7[i] = alg.d61 * cache.du[i] + alg.d66 * cache.k6[i] + alg.d67 * cache.k7[i] +
                alg.d68 * cache.k8[i] + alg.d69 * cache.k9[i] + alg.d610 * cache.k10[i] +
                alg.d611 * cache.k2[i] + alg.d612 * cache.k3[i];
        cache.rcont8[i] = alg.d71 * cache.du[i] + alg.d76 * cache.k6[i] + alg.d77 * cache.k7[i] +
                alg.d78 * cache.k8[i] + alg.d79 * cache.k9[i] + alg.d710 * cache.k10[i] +
                alg.d711 * cache.k2[i] + alg.d712 * cache.k3[i];
    }

    /* the next three function evaluations */
    for (int i = 0; i < integrator.n; i++)
        cache.unew[i] = integrator.u[i] + integrator.dt * (alg.a141 * cache.du[i] + alg.a147 * cache.k7[i] +
                alg.a148 * cache.k8[i] + alg.a149 * cache.k9[i] + alg.a1410 * cache.k10[i] +
                alg.a1411 * cache.k2[i] + alg.a1412 * cache.k3[i] + alg.a1413 * cache.k4[i]);
    integrator.f->dudt(cache.k10, cache.unew, integrator.t + alg.c14 * integrator.dt);
    for (int i = 0; i < integrator.n; i++)
        cache.unew[i] = integrator.u[i] + integrator.dt * (alg.a151 * cache.du[i] + alg.a156 * cache.k6[i] +
                alg.a157 * cache.k7[i] + alg.a158 * cache.k8[i] + alg.a1511 * cache.k2[i] +
                alg.a1512 * cache.k3[i] + alg.a1513 * cache.k4[i] + alg.a1514 * cache.k10[i]);
    integrator.f->dudt(cache.k2, cache.unew, integrator.t + alg.c15 * integrator.dt);
    for (int i = 0; i < integrator.n; i++)
        cache.unew[i] = integrator.u[i] + integrator.dt * (alg.a161 * cache.du[i] + alg.a166 * cache.k6[i] +
                alg.a167 * cache.k7[i] + alg.a168 * cache.k8[i] + alg.a169 * cache.k9[i] +
                alg.a1613 * cache.k4[i] + alg.a1614 * cache.k10[i] + alg.a1615 * cache.k2[i]);
    integrator.f->dudt(cache.k3, cache.unew, integrator.t + alg.c16 * integrator.dt);
    integrator.num_function_evaluations += 3;

    // Final preperation
    for (int i = 0; i < integrator.n; i++) {
        cache.rcont5[i] = integrator.dt * (cache.rcont5[i] + alg.d413 * cache.k4[i] + alg.d414 * cache.k10[i] +
                alg.d415 * cache.k2[i] + alg.d416 * cache.k3[i]);
        cache.rcont6[i] = integrator.dt * (cache.rcont6[i] + alg.d513 * cache.k4[i] + alg.d514 * cache.k10[i] +
                alg.d515 * cache.k2[i] + alg.d516 * cache.k3[i]);
        cache.rcont7[i] = integrator.dt * (cache.rcont7[i] + alg.d613 * cache.k4[i] + alg.d614 * cache.k10[i] +
                alg.d615 * cache.k2[i] + alg.d616 * cache.k3[i]);
        cache.rcont8[i] = integrator.dt * (cache.rcont8[i] + alg.d713 * cache.k4[i] + alg.d714 * cache.k10[i] +
                alg.d715 * cache.k2[i] + alg.d716 * cache.k3[i]);
    }
}

/**
 * Compute the dense output for the ith component of the solution between times
 * t and t + dt
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince5 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince5
 * @param i Component of solution to compute dense output for
 * @param t Current time value
 * @param dt Current step value
 * @return Solution value at time t
 */
static double dense_output(
        ODEIntegrator &integrator,
        DormandPrince8 &alg,
        DormandPrince8Cache &cache,
        int i,
        double t_t,
        double t_dt
) {
    double s = (t_t - integrator.tprev) / t_dt;
    double s1 = 1.0 - s;

    return cache.rcont1[i] +
            s * (cache.rcont2[i] + s1 * (cache.rcont3[i] + s * (cache.rcont4[i] + s1 * (cache.rcont5[i] +
                    s * (cache.rcont6[i] + s1 * (cache.rcont7[i] + s * cache.rcont8[i]))))));
}


/**
 * Prepare the integrator and cache for the next step.
 * @param t_err Estimated error from step just taken
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince8 algorithm
 * @param cache Cache for vectors/values needed for DormandPrince8
 */
static void prepare_next_step(
        const double t_err,
        ODEIntegrator &integrator,
        DormandPrince8 &alg,
        DormandPrince8Cache &cache
) {
    /* computation of hnew */
    double fac11 = pow(t_err, alg.alpha);
    /* Lund-stabilization */
    double fac = fac11 / pow(cache.next_prev_step_ratio, alg.beta);
    /* we require fac1 <= hnew/h <= fac2 */
    fac = std::max(alg.maxNextPrevStepRatio, std::min(alg.minNextPrevStepRatio, fac / alg.safe));
    double dtnew = integrator.dt / fac;

    if (t_err <= 1.0) {
        // step accepted
        cache.next_prev_step_ratio = std::max(t_err, 1.0E-4);
        integrator.num_accept++;
        integrator.f->dudt(cache.k4, cache.k5, integrator.t + integrator.dt);
        integrator.num_function_evaluations++;

        /* stiffness detection */
        if (!(integrator.num_accept % cache.max_num_stiff) || (cache.num_stiff > 0)) {
            double stnum = 0.0;
            double stden = 0.0;
            for (int i = 0; i < integrator.n; i++) {
                double sqr = cache.k4[i] - cache.k3[i];
                stnum += sqr * sqr;
                sqr = cache.k5[i] - cache.unew[i];
                stden += sqr * sqr;
            }
            if (stden > 0.0) {
                cache.dt_lam_b = integrator.dt * sqrt(stnum / stden);
            }
            if (cache.dt_lam_b > 6.1) {
                cache.num_non_stiff = 0;
                cache.num_stiff++;
                if (cache.num_stiff == 15) {
                    throw std::runtime_error("ODE seems to be stiff.");
                }
            } else {
                cache.num_non_stiff++;
                if (cache.num_non_stiff == 6) {
                    cache.num_stiff = 0;
                }
            }
        }

        if (integrator.opts.dense) {
            prepare_dense(integrator, alg, cache);
        }

        cache.du = cache.k4;
        integrator.u = cache.k5;
        integrator.tprev = integrator.t;
        integrator.t += integrator.dt;

        if (fabs(dtnew) > integrator.opts.dtmax) {
            dtnew = integrator.tdir * integrator.opts.dtmax;
        }
        if (cache.reject) {
            dtnew = integrator.tdir * std::min(fabs(dtnew), fabs(integrator.dt));
        }
        cache.reject = false;
    } else {
        /* step rejected */
        dtnew = integrator.dt / std::min(cache.next_prev_step_ratio, fac11 / alg.safe);
        cache.reject = true;
        if (integrator.num_accept >= 1) {
            integrator.num_reject++;
        }
        cache.last = false;
    }
    integrator.dtprev = integrator.dt;
    integrator.dt = dtnew;
}


/**
 * Perform a single step using DormandPrince8
 * @param t_err Estimated error from step just taken
 * @param integrator Driver for ODE integration
 * @param alg DormandPrince8 algorithm
 */
static void step(ODEIntegrator &integrator, DormandPrince8 &alg, DormandPrince8Cache &cache) {
    while (true) {
        compute_stages(integrator, alg, cache);
        double err = error(integrator, alg, cache);
        prepare_next_step(err, integrator, alg, cache);
        if (!cache.reject) {
            break;
        }
        if (abs(integrator.dt) <= abs(integrator.t) * std::numeric_limits<double>::epsilon())
            throw std::runtime_error("Step size underflow in Dormand-Prince-8");
    }
    integrator.num_steps++;
}


/**
 * Solve the ODE problem using the Dormand-Prince 5th order algorithm
 * @param prob ODE problem to solve
 * @param alg DormandPrince8 algorithm
 * @param opts Integration options
 * @return Solution
 */
ODESolution solve(
        ODEProblem &prob,
        DormandPrince8 &alg,
        ODEIntegratorOptions &opts
) {

    ODEIntegrator integrator{prob, opts};
    DormandPrince8Cache cache{integrator.n};


    integrator.sol.ts.push_back(integrator.t);
    integrator.sol.us.push_back(integrator.u);

    // Compute du
    integrator.f->dudt(cache.du, integrator.u, integrator.t);

    // Compute solution by stepping until we reach end of time interval.
    while (!cache.last) {

        if ((integrator.t + 1.01 * integrator.dt - prob.t_span.second) * integrator.tdir > 0.0) {
            integrator.dt = prob.t_span.second - integrator.t;
            cache.last = true;
        }

        step(integrator, alg, cache);
        integrator.sol.ts.push_back(integrator.t);
        integrator.sol.us.push_back(integrator.u);

        if (integrator.num_steps > integrator.opts.max_num_steps) {
            throw std::runtime_error("too many steps..");
        }
    }
    return integrator.sol;
}

/**
 * Solve the ODE problem using the Dormand-Prince 5th order algorithm
 * @param prob ODE problem to solve
 * @param alg DormandPrince8 algorithm
 * @return Solution
 */
ODESolution solve(ODEProblem &prob, DormandPrince8 &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}

}
}

#endif //LANRE_DIFFEQ_DORMAND_PRINCE_HPP
