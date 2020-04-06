//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_INTEGRATOR_HPP
#define LANRE_DIFFEQ_INTEGRATOR_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/solution.hpp"

namespace lanre {
namespace diffeq {

/*
 * Struct containing the common solver options
 */
struct ODEIntegratorOptions {
    double reltol = 0.0;
    double abstol = 0.0;
    bool dense = false;
    double dtstart = 0.0;
    double dtmax = 0.0, dtmin = std::numeric_limits<double>::epsilon();
    unsigned int max_num_steps = 100000;
};

/**
 * Driving structure for solving ODEs
 */
struct ODEIntegrator {
    // Required
    std::shared_ptr<ODEFunction> f;
    int n;
    Vector<double> u, uprev;

    // Integrator options
    ODEIntegratorOptions opts;

    double t, tprev, tdir, dt, dtprev, t_final;
    ODESolution sol;

    // Counting variables
    unsigned int num_steps;
    unsigned int num_function_evaluations;
    unsigned int num_jacobian_evaluations;
    unsigned int num_decompositions;
    unsigned int num_linear_solves;
    unsigned int num_accept;
    unsigned int num_reject;

    ODEIntegrator(ODEProblem &t_prob)
            : f(t_prob.ode_function), n(t_prob.u_init.size()),
              u(t_prob.u_init), uprev(Vector<double>{n}),
              opts(ODEIntegratorOptions{}),
              t(t_prob.t_span.first), tprev(t_prob.t_span.first),
              tdir(std::copysign(1.0, t_prob.t_span.second - t_prob.t_span.first)),
              dt(opts.dtstart), sol(ODESolution{}), dtprev(dt), t_final(t_prob.t_span.second),
              num_steps(0), num_function_evaluations(0), num_jacobian_evaluations(0),
              num_linear_solves(0), num_accept(0), num_reject(0) {
        if (dt <= 0.0) {
            dt = 1e-6;
        }
        if (opts.reltol <= 0.0) {
            opts.reltol = 1e-5;
        }
        if (opts.abstol <= 0.0) {
            opts.abstol = 1e-5;
        }
        if (opts.dtmax <= 0.0) {
            opts.dtmax = t_prob.t_span.second - t_prob.t_span.first;
        }
    }

    ODEIntegrator(ODEProblem &t_prob, ODEIntegratorOptions &t_opts)
            : f(t_prob.ode_function), n(t_prob.u_init.size()),
              u(t_prob.u_init), uprev(Vector<double>{n}),
              t(t_prob.t_span.first), tprev(t_prob.t_span.first),
              tdir(std::copysign(1.0, t_prob.t_span.second - t_prob.t_span.first)),
              dt(t_opts.dtstart), sol(ODESolution{}), dtprev(dt),
              t_final(t_prob.t_span.second), opts(t_opts),
              num_steps(0), num_function_evaluations(0), num_jacobian_evaluations(0),
              num_linear_solves(0), num_accept(0), num_reject(0) {
        if (dt <= 0.0) {
            dt = 1e-6;
        }
        if (opts.reltol <= 0.0) {
            opts.reltol = 1e-5;
        }
        if (opts.abstol <= 0.0) {
            opts.abstol = 1e-5;
        }
        if (opts.dtmax <= 0.0) {
            opts.dtmax = t_prob.t_span.second - t_prob.t_span.first;
        }
    }
};

}
}

#endif //LANRE_DIFFEQ_INTEGRATOR_HPP
