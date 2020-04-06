//
// Created by Logan Morrison on 3/17/20.
//

#ifndef LANRE_DIFFEQ_SEULEX_HPP
#define LANRE_DIFFEQ_SEULEX_HPP

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

struct Seulex {
    const int kmax = 12;
    const int imax = kmax + 1;
    const double fac1 = 0.6;
    const double fac2 = 0.93;
    const double fac3 = 0.1;
    const double fac4 = 4.0;
    const double fac5 = 0.5;
    const double kfac1 = 0.7;
    const double kfac2 = 0.9;
};


struct SeulexCache {
    double theta;
    int ipt, kright, ktarg;

    bool first_step = true, last_step = false, jac_redo;
    bool reject = false, prev_reject = false, caljac = false;
    double errold, dtnext, dtdid;

    Vector<double> uend, scale, fsave, cost;
    Vector<double> utemp, du, dutemp, del, dens, factrl;
    Vector<double> dtopt, work, usav, useq, umid;
    Vector<double> dfdt;
    Vector<int> nseq;
    Matrix<double> a;
    Matrix<double> dfdu, table, coeff;
};


static void init(ODEIntegrator &integrator, Seulex &alg, SeulexCache &cache) {
    cache.uend.resize(integrator.n);
    cache.scale.resize(integrator.n);
    cache.fsave.resize((alg.imax * alg.imax - 1) / 2 + 2, integrator.n);
    cache.cost.resize(alg.imax);
    cache.utemp.resize(integrator.n);
    cache.du.resize(integrator.n);
    cache.dutemp.resize(integrator.n);
    cache.del.resize(integrator.n);
    cache.dens.resize((alg.imax + 2) * integrator.n);
    cache.factrl.resize(alg.imax);
    cache.dtopt.resize(alg.imax);
    cache.work.resize(alg.imax);
    cache.usav.resize(integrator.n);
    cache.useq.resize(integrator.n);
    cache.umid.resize(integrator.n);
    cache.dfdt.resize(integrator.n);
    cache.nseq.resize(alg.imax);
    cache.a.resize(integrator.n, integrator.n);
    cache.dfdu.resize(integrator.n, integrator.n);
    cache.table.resize(alg.kmax, integrator.n);
    cache.coeff.resize(alg.imax, alg.imax);


    static const double costfunc = 1.0, costjac = 5.0, costlu = 1.0, costsolve = 1.0;
    cache.jac_redo = std::min(1.0e-4, integrator.opts.reltol);
    cache.theta = 2.0 * cache.jac_redo;
    cache.nseq[0] = 2;
    cache.nseq[1] = 3;
    for (int i = 2; i < alg.imax; i++) {
        cache.nseq[i] = 2 * cache.nseq[i - 2];
    }
    cache.cost[0] = costjac + costlu + cache.nseq[0] * (costfunc + costsolve);
    for (int k = 0; k < alg.kmax; k++) {
        cache.cost[k + 1] = cache.cost[k] + (cache.nseq[k + 1] - 1) * (costfunc + costsolve) + costlu;
    }
    cache.dtnext = -1.0e99;
    double logfact = -log10(integrator.opts.reltol + integrator.opts.abstol) * 0.6 + 0.5;
    cache.ktarg = std::max(1, std::min(alg.kmax - 1, int(logfact)));
    for (int k = 0; k < alg.imax; k++) {
        for (int l = 0; l < k; l++) {
            double ratio = double(cache.nseq[k]) / cache.nseq[l];
            cache.coeff(k, l) = 1.0 / (ratio - 1.0);
        }
    }
    cache.factrl[0] = 1.0;
    for (int k = 0; k < alg.imax - 1; k++) {
        cache.factrl[k + 1] = (k + 1) * cache.factrl[k];
    }
}


static bool compute_du(
        ODEIntegrator &integrator,
        Seulex &alg,
        SeulexCache &cache,
        Vector<double> &u,
        const double dttot,
        const int k,
        Vector<double> &uend,
        int &ipt,
        Vector<double> &scale
) {
    Vector<double> del(integrator.n);
    Vector<double> utemp(integrator.n);
    Vector<double> dutemp(integrator.n);

    int nstep = cache.nseq[k];
    double dt = dttot / nstep;

    for (int i = 0; i < integrator.n; i++) {
        for (int j = 0; j < integrator.n; j++) cache.a(i, j) = -cache.dfdu(i, j);
        cache.a(i, i) += 1.0 / dt;
    }
    //LUdcmp alu(cache.a);
    auto alu = cache.a.lu();
    double xnew = integrator.t + dt;
    integrator.f->dudt(del, u, xnew);
    // derivs\((.*?),(.*?),(.*?)\)
    // integrator.f->dudt($3,$2,$1)
    for (int i = 0; i < integrator.n; i++)
        utemp[i] = u[i];
    del = alu.solve(del);
    if (integrator.opts.dense && nstep == k + 1) {
        ipt++;
        for (int i = 0; i < integrator.n; i++)
            cache.fsave(ipt, i) = del[i];
    }
    for (int nn = 1; nn < nstep; nn++) {
        for (int i = 0; i < integrator.n; i++)
            utemp[i] += del[i];
        xnew += dt;
        integrator.f->dudt(uend, utemp, xnew);
        if (nn == 1 && k <= 1) {
            double del1 = 0.0;
            for (int i = 0; i < integrator.n; i++)
                del1 += boost::math::pow<2>(del[i] / scale[i]);
            del1 = sqrt(del1);
            integrator.f->dudt(dutemp, utemp, integrator.t + dt);
            for (int i = 0; i < integrator.n; i++)
                del[i] = dutemp[i] - del[i] / dt;
            del = alu.solve(del);
            double del2 = 0.0;
            for (int i = 0; i < integrator.n; i++)
                del2 += boost::math::pow<2>(del[i] / scale[i]);
            del2 = sqrt(del2);
            cache.theta = del2 / std::max(1.0, del1);
            if (cache.theta > 1.0)
                return false;
        }
        del = alu.solve(uend);
        if (integrator.opts.dense && nn >= nstep - k - 1) {
            ipt++;
            for (int i = 0; i < integrator.n; i++)
                cache.fsave(ipt, i) = del[i];
        }
    }
    for (int i = 0; i < integrator.n; i++)
        uend[i] = utemp[i] + del[i];
    return true;
}

static void polynomial_extrapolate(
        ODEIntegrator &integrator,
        Seulex &alg,
        SeulexCache &cache,
        const int k,
        Matrix<double> &table,
        Vector<double> &last
) {
    int l = last.size();
    for (int j = k - 1; j > 0; j--)
        for (int i = 0; i < l; i++)
            table(j - 1, i) = table(j, i) + cache.coeff(k, j) * (table(j, i) - table(j - 1, i));
    for (int i = 0; i < l; i++)
        last[i] = table(0, i) + cache.coeff(k, 0) * (table(0, i) - last[i]);
}


static double prepare_dense(
        ODEIntegrator &integrator,
        Seulex &alg,
        SeulexCache &cache,
        Vector<double> &usav,
        Vector<double> &scale,
        const int k,
        double &error
) {
    cache.kright = k;
    for (int i = 0; i < integrator.n; i++) {
        cache.dens[i] = usav[i];
        cache.dens[integrator.n + i] = integrator.u[i];
    }
    for (int klr = 0; klr < cache.kright; klr++) {
        if (klr >= 1) {
            for (int kk = klr; kk <= k; kk++) {
                int lbeg = ((kk + 3) * kk) / 2;
                int lend = lbeg - kk + 1;
                for (int l = lbeg; l >= lend; l--)
                    for (int i = 0; i < integrator.n; i++)
                        cache.fsave(l, i) -= -cache.fsave(l - 1, i);
            }
        }
        for (int kk = klr; kk <= k; kk++) {
            double facnj = cache.nseq[kk];
            facnj = pow(facnj, klr + 1) / cache.factrl[klr + 1];
            int ipt = ((kk + 3) * kk) / 2;
            int krn = (kk + 2) * integrator.n;
            for (int i = 0; i < integrator.n; i++) {
                cache.dens[krn + i] = cache.fsave(ipt, i) * facnj;
            }
        }
        for (int j = klr + 1; j <= k; j++) {
            double dblenj = cache.nseq[j];
            for (int l = j; l >= klr + 1; l--) {
                double factor = dblenj / cache.nseq[l - 1] - 1.0;
                for (int i = 0; i < integrator.n; i++) {
                    int krn = (l + 2) * integrator.n + i;
                    cache.dens[krn - integrator.n] =
                            cache.dens[krn] + (cache.dens[krn] - cache.dens[krn - integrator.n]) / factor;
                }
            }
        }
    }
    for (int in = 0; in < integrator.n; in++) {
        for (int j = 1; j <= cache.kright + 1; j++) {
            int ii = integrator.n * j + in;
            cache.dens[ii] = cache.dens[ii] - cache.dens[ii - integrator.n];
        }
    }
}

static double dense_out(
        ODEIntegrator &integrator,
        Seulex &alg,
        SeulexCache &cache,
        const int i,
        const double t,
        const double dt
) {
    double theta = (t - integrator.tprev) / dt;
    int k = cache.kright;
    double yinterp = cache.dens[(k + 1) * integrator.n + i];
    for (int j = 1; j <= k; j++)
        yinterp = cache.dens[(k + 1 - j) * integrator.n + i] + yinterp * (theta - 1.0);
    return cache.dens[i] + yinterp * theta;
}

static void step(
        ODEIntegrator &integrator,
        Seulex &alg,
        SeulexCache &cache,
        double dttry
) {
    int i, k;
    double fac, dt, dtnew, err;
    bool firstk;

    Vector<double> dtopt(alg.imax);
    Vector<double> work(alg.imax);
    Vector<double> usav(integrator.n);
    Vector<double> useq(integrator.n);
    Vector<double> umid(integrator.n);
    Vector<double> scale(integrator.n);

    work[0] = 1.e30;
    dt = dttry;
    for (i = 0; i < integrator.n; i++) {
        usav[i] = integrator.u[i];
    }
    if (dt != cache.dtnext && !cache.first_step) {
        cache.last_step = true;
    }
    if (cache.reject) {
        cache.prev_reject = true;
        cache.last_step = false;
        cache.theta = 2.0 * cache.jac_redo;
    }
    for (i = 0; i < integrator.n; i++)
        scale[i] = integrator.opts.abstol + integrator.opts.reltol * abs(integrator.u[i]);
    cache.reject = false;
    firstk = true;
    dtnew = abs(dt);
    compute_jac:
    if (cache.theta > cache.jac_redo && !cache.caljac) {
        integrator.f->dfdu(cache.dfdu, integrator.u, integrator.t);
        integrator.f->dfdt(cache.dfdt, integrator.u, integrator.t);
        cache.caljac = true;
    }
    while (firstk || cache.reject) {
        dt = integrator.tdir * dtnew;
        firstk = false;
        cache.reject = false;
        if (abs(dt) <= abs(integrator.t) * std::numeric_limits<double>::epsilon())
            throw std::runtime_error("Seulex: Step size too small.");
        int ipt = -1;
        for (k = 0; k <= cache.ktarg + 1; k++) {
            bool success = compute_du(integrator, alg, cache, cache.usav, dt, k, cache.useq, ipt, scale);
            if (!success) {
                cache.reject = true;
                dtnew = abs(dt) * alg.fac5;
                break;
            }
            if (k == 0)
                integrator.u = cache.useq;
            else
                for (i = 0; i < integrator.n; i++)
                    cache.table(k - 1, i) = cache.useq[i];
            if (k != 0) {
                polynomial_extrapolate(integrator, alg, cache, k, cache.table, integrator.u);
                err = 0.0;
                for (i = 0; i < integrator.n; i++) {
                    scale[i] = integrator.opts.abstol + integrator.opts.reltol * abs(cache.usav[i]);
                    err += boost::math::pow<2>((integrator.u[i] - cache.table(0, i)) / scale[i]);
                }
                err = sqrt(err / integrator.n);
                if (err > 1.0 / std::numeric_limits<double>::epsilon() || (k > 1 && err >= cache.errold)) {
                    cache.reject = true;
                    dtnew = abs(dt) * alg.fac5;
                    break;
                }
                cache.errold = std::max(4.0 * err, 1.0);
                double expo = 1.0 / (k + 1);
                double facmin = pow(alg.fac3, expo);
                if (err == 0.0)
                    fac = 1.0 / facmin;
                else {
                    fac = alg.fac2 / pow(err / alg.fac1, expo);
                    fac = std::max(facmin / alg.fac4, std::min(1.0 / facmin, fac));
                }
                dtopt[k] = abs(dt * fac);
                work[k] = cache.cost[k] / dtopt[k];
                if ((cache.first_step || cache.last_step) && err <= 1.0)
                    break;
                if (k == cache.ktarg - 1 && !cache.prev_reject && !cache.first_step && !cache.last_step) {
                    if (err <= 1.0)
                        break;
                    else if (err > cache.nseq[cache.ktarg] * cache.nseq[cache.ktarg + 1] * 4.0) {
                        cache.reject = true;
                        cache.ktarg = k;
                        if (cache.ktarg > 1 && work[k - 1] < alg.kfac1 * work[k])
                            cache.ktarg--;
                        dtnew = dtopt[cache.ktarg];
                        break;
                    }
                }
                if (k == cache.ktarg) {
                    if (err <= 1.0)
                        break;
                    else if (err > cache.nseq[k + 1] * 2.0) {
                        cache.reject = true;
                        if (cache.ktarg > 1 && work[k - 1] < alg.kfac1 * work[k])
                            cache.ktarg--;
                        dtnew = dtopt[cache.ktarg];
                        break;
                    }
                }
                if (k == cache.ktarg + 1) {
                    if (err > 1.0) {
                        cache.reject = true;
                        if (cache.ktarg > 1 && work[cache.ktarg - 1] < alg.kfac1 * work[cache.ktarg])
                            cache.ktarg--;
                        dtnew = dtopt[cache.ktarg];
                    }
                    break;
                }
            }
        }
        if (cache.reject) {
            cache.prev_reject = true;
            if (!cache.caljac) {
                cache.theta = 2.0 * cache.jac_redo;
                goto compute_jac;
            }
        }
    }
    cache.caljac = false;
    if (integrator.opts.dense)
        prepare_dense(integrator, alg, cache, cache.usav, scale, k, err);
    integrator.tprev = integrator.t;
    integrator.t += dt;
    cache.dtdid = dt;
    cache.first_step = false;
    int kopt;
    if (k == 1)
        kopt = 2;
    else if (k <= cache.ktarg) {
        kopt = k;
        if (work[k - 1] < alg.kfac1 * work[k])
            kopt = k - 1;
        else if (work[k] < alg.kfac2 * work[k - 1])
            kopt = std::min(k + 1, alg.kmax - 1);
    } else {
        kopt = k - 1;
        if (k > 2 && work[k - 2] < alg.kfac1 * work[k - 1])
            kopt = k - 2;
        if (work[k] < alg.kfac2 * work[kopt])
            kopt = std::min(k, alg.kmax - 1);
    }
    if (cache.prev_reject) {
        cache.ktarg = std::min(kopt, k);
        dtnew = std::min(abs(dt), dtopt[cache.ktarg]);
        cache.prev_reject = false;
    } else {
        if (kopt <= k)
            dtnew = dtopt[kopt];
        else {
            if (k < cache.ktarg && work[k] < alg.kfac2 * work[k - 1])
                dtnew = dtopt[k] * cache.cost[kopt + 1] / cache.cost[k];
            else
                dtnew = dtopt[k] * cache.cost[kopt] / cache.cost[k];
        }
        cache.ktarg = kopt;
    }
    cache.dtnext = integrator.tdir * dtnew;
}


static ODESolution seulex_core(ODEIntegrator &integrator, Seulex &alg, SeulexCache &cache) {
    ODESolution solution{};
    init(integrator, alg, cache);

    integrator.f->dudt(cache.du, integrator.u, integrator.t);

    solution.ts.push_back(integrator.t);
    solution.us.push_back(integrator.u);

    double tstart = integrator.t;

    for (integrator.num_steps = 0; integrator.num_steps < integrator.opts.max_num_steps; integrator.num_steps++) {
        if ((integrator.t + integrator.dt * 1.0001 - integrator.t_final) * (integrator.t_final - tstart) > 0.0)
            integrator.dt = integrator.t_final - integrator.t;
        step(integrator, alg, cache, integrator.dt);
        if (cache.dtdid == integrator.dt)
            ++integrator.num_accept;
        else
            ++integrator.num_reject;

        solution.ts.push_back(integrator.t);
        solution.us.push_back(integrator.u);

        if ((integrator.t - integrator.t_final) * (integrator.t_final - tstart) >= 0.0) {
            return solution;
        }
        if (abs(cache.dtnext) <= integrator.opts.dtmin)
            throw std::runtime_error("Seulex: Step size too small.");
        integrator.dtprev = integrator.dt;
        integrator.dt = cache.dtnext;
    }
    throw std::runtime_error("Seulex: Too many steps.");
}


ODESolution solve(ODEProblem &prob, Seulex &alg, ODEIntegratorOptions &opts) {
    ODEIntegrator integrator{prob, opts};
    SeulexCache cache{};
    return seulex_core(integrator, alg, cache);
}

ODESolution solve(ODEProblem &prob, Seulex &alg) {
    ODEIntegratorOptions opts{};
    return solve(prob, alg, opts);
}

}
}

#endif //LANRE_DIFFEQ_SEULEX_HPP
