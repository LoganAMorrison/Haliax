//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_SOLUTION_HPP
#define LANRE_DIFFEQ_SOLUTION_HPP

#include "lanre/diffeq/base.hpp"
#include <vector>

namespace lanre {
namespace diffeq {

enum Retcode {
    Default,
    Success,
    MaxIters,
    Unstable,
    DtLessThanMin,
    LinearAlgError,
    SingularMatrix,
    Failure
};

struct ODESolution {
    std::vector<double> ts{};
    std::vector<Vector < double>> us{};
    Retcode retcode = Default;

    ODESolution() {
        ts.reserve(500);
        us.reserve(500);
    }

};

}
}

#endif //LANRE_DIFFEQ_SOLUTION_HPP
