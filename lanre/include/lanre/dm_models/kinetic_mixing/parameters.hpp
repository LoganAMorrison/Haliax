//
// Created by Logan Morrison on 4/11/20.
//


#ifndef LANRE_KINETIC_MIXING_PARAMETERS_HPP
#define LANRE_KINETIC_MIXING_PARAMETERS_HPP

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

struct Parameters {
    double mx;
    double mv;
    double gvxx;
    double eps;
    double widthv;
};

}
}
}

#endif //LANRE_KINETIC_MIXING_PARAMETERS_HPP

