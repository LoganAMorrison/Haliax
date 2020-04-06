//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_BASE_HPP
#define LANRE_INTEGRATE_BASE_HPP

#include <limits>

namespace lanre {
namespace integrate {

static const int kLIMIT = 500;
static const double kEPMACH = std::numeric_limits<double>::epsilon();
static const double kUFLOW = std::numeric_limits<double>::min();
static const double kOFLOW = std::numeric_limits<double>::max();
static const int kCOSINE = 1;
static const int kSINE = 2;

}
}

#endif //LANRE_INTEGRATE_BASE_HPP
