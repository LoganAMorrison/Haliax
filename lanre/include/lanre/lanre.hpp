/*
 * Created by Logan Morrison on 4/11/20.
 *
 * Contains common definitions used throughout this library.
 */

#ifndef LANRE_HPP
#define LANRE_HPP

#include <Eigen/Dense>

namespace lanre {

template<class Type, int n = Eigen::Dynamic>
using Vector = Eigen::Matrix<Type, n, 1>;

template<class Type, int n = Eigen::Dynamic>
using RowVector = Eigen::Matrix<Type, 1, n>;

template<class Type, int n = Eigen::Dynamic, int m = Eigen::Dynamic>
using Matrix = Eigen::Matrix<Type, n, m>;

template<class Type, int n = Eigen::Dynamic, int m = Eigen::Dynamic>
using Array = Eigen::Array<Type, n, m>;

}


#endif //LANRE_HPP
