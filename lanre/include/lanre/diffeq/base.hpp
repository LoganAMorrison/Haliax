//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_BASE_HPP
#define LANRE_DIFFEQ_BASE_HPP

#include <Eigen/Dense>

namespace lanre {
namespace diffeq {

template<class T>
using Matrix =Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<class T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

}
}

#endif //LANRE_DIFFEQ_BASE_HPP
