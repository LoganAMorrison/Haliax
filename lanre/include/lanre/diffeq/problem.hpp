//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_PROBLEM_HPP
#define LANRE_DIFFEQ_PROBLEM_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/function.hpp"
#include <memory>

namespace lanre {
namespace diffeq {

struct ODEProblem {
    std::shared_ptr<ODEFunction> ode_function;
    Vector<double> u_init;
    std::pair<double, double> t_span;
};

}
}

#endif //LANRE_DIFFEQ_PROBLEM_HPP
