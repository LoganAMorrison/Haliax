//
// Created by Logan Morrison on 4/12/20.
//

#ifndef LANRE_COSMOLOGY_THERMODYNAMIC_FUNCTIONS_HPP
#define LANRE_COSMOLOGY_THERMODYNAMIC_FUNCTIONS_HPP

#include "lanre/autodiff/dual.hpp"
#include "lanre/special_functions/besselk.hpp"
#include <cmath>
#include <utility>

namespace lanre {
namespace cosmology {

/**
 * Number density scaled by the particle's temperature and internal d.o.f.
 * @param x mass divided by temperature: m / T
 * @param spin2 Two times the spin of the particle.
 * @return neq(T) / (g T^3)
 */
template<class Type>
static Type neq_scaled(Type x, unsigned int spin2, double abserr = 1e-12, int maxiter = 15) {
    using lanre::special_functions::besselk2;
    if (x == 0.0)
        return (spin2 % 2 == 0) ? static_cast<Type>(0.121793828233573) : static_cast<Type>(0.0913453711751798);
    else {
        Type eta = static_cast<Type>((spin2 % 2 == 0) ? 1.0 : -1.0);

        Type abserr_t = static_cast<Type>(abserr);
        Type sum = static_cast<Type>(0.0);
        Type sum_old, err;
        int n = 1;
        do {
            sum_old = sum;
            sum += pow(eta, n + 1) * besselk2(n * x) / n;
            err = abs(sum - sum_old);
            n++;
        } while (err > abserr_t && n <= maxiter);

        return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
    }
}

/**
 * Energy density scaled by the particle's temperature and internal d.o.f.
 * @param x mass divided by temperature: m / T
 * @param spin2 Two times the spin of the particle.
 * @return energy_density(T) / (g T^4)
 */
template<class Type>
static Type energy_density_scaled(Type x, unsigned int spin2, double abserr = 1e-12, int maxiter = 5) {
    using lanre::special_functions::besselk1;
    using lanre::special_functions::besselk2;
    if (x == 0.0)
        return static_cast<Type>(M_PI * M_PI / 30.0 * ((spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0));
    else {
        Type eta = static_cast<Type>((spin2 % 2 == 0) ? 1.0 : -1.0);

        Type abserr_t = static_cast<Type>(abserr);
        Type sum = static_cast<Type>(0.0);
        Type sum_old, err;
        int n = 1;
        do {
            sum_old = sum;
            sum += pow(eta, n + 1) / (n * n) * (x * n * besselk1(n * x) + 3.0 * besselk2(n * x));
            err = abs(sum - sum_old);
            n++;
        } while (err > abserr_t && n <= maxiter);

        return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
    }
}

/**
 * Pressure density scaled by the particle's temperature and internal d.o.f.
 * @param x mass divided by temperature: m / T
 * @param spin2 Two times the spin of the particle.
 * @return pressure_density(T) / (g T^4)
 */
template<class Type>
static Type pressure_density_scaled(Type x, unsigned int spin2, double abserr = 1e-12, int maxiter = 5) {
    using lanre::special_functions::besselk2;
    if (x == 0.0)
        return static_cast<Type>(M_PI * M_PI / 90.0 * ((spin2 % 2 == 0) ? 1 : 7.0 / 8.0));
    else {
        Type eta = static_cast<Type>((spin2 % 2 == 0) ? 1.0 : -1.0);

        Type abserr_t = static_cast<Type>(abserr);
        Type sum = static_cast<Type>(0.0);
        Type sum_old, err;
        int n = 1;
        do {
            sum_old = sum;
            sum += pow(eta, n + 1) / (n * n) * besselk2(n * x);
            err = abs(sum - sum_old);
            n++;
        } while (err > abserr_t && n <= maxiter);

        return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
    }
}

/**
 * Entropy density scaled by the particle's temperature and internal d.o.f.
 * @param x mass divided by temperature: m / T
 * @param spin2 Two times the spin of the particle.
 * @return entropy_density(T) / (g T^3)
 */
template<class Type>
static Type entropy_density_scaled(Type x, unsigned int spin2, double abserr = 1e-12, int maxiter = 5) {
    using lanre::special_functions::besselk3;
    if (x == static_cast<Type>(0.0))
        return static_cast<Type>(2.0 * M_PI * M_PI / 45.0 * ((spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0));
    else {
        Type eta = static_cast<Type>((spin2 % 2 == 0) ? 1.0 : -1.0);

        Type abserr_t = static_cast<Type>(abserr);
        Type sum = static_cast<Type>(0.0);
        Type sum_old, err;
        int n = 1;
        do {
            sum_old = sum;
            sum += pow(eta, n + 1) / n * besselk3(n * x);
            err = abs(sum - sum_old);
            n++;
        } while (err > abserr_t && n <= maxiter);

        return x * x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
    }
}


/**
 * Compute the equilibrium number density of particle with temperature T.
 * @param T Temperature of particle.
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return Number density
 */
template<class TypeT, class TypeM>
auto neq(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return g * neq_scaled(mass / T, spin2, abserr, maxiter) * T * T * T;
}


/**
 * Compute the equilibrium energy density of particle with temperature T.
 * @param T Temperature of particle
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return Energy density
 */
template<class TypeT, class TypeM>
auto energy_density(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return g * energy_density_scaled(mass / T, spin2, abserr, maxiter) * T * T * T * T;
}

/**
 * Compute the equilibrium pressure density of particle with temperature T.
 * @param T Temperature of particle
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return Pressure density
 */
template<class TypeT, class TypeM>
auto pressure_density(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return g * pressure_density_scaled(mass / T, spin2, abserr, maxiter) * T * T * T * T;
}

/**
 * Compute the equilibrium entropy density of particle with temperature T.
 * @param T Temperature of particle
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return Entropy density
 */
template<class TypeT, class TypeM>
auto entropy_density(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return g * entropy_density_scaled(mass / T, spin2, abserr, maxiter) * T * T * T;
}

/**
 * Compute the equilibrium number of d.o.f. stored in energy of a particle
 * with temperature T
 * @param T Temperature of particle
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return D.o.f. stored in energy
 */
template<class TypeT, class TypeM>
auto geff(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return 30.0 / (M_PI * M_PI) * g * energy_density_scaled(mass / T, spin2, abserr, maxiter);
}

/**
 * Compute the equilibrium number of d.o.f. stored in entropy of a particle
 * with temperature T
 * @param T Temperature of particle
 * @param mass Mass of the particle.
 * @param g Internal degrees of freedom of the particle.
 * @param spin2 Two times the spin of the particle.
 * @return D.o.f. stored in entropy
 */
template<class TypeT, class TypeM>
auto heff(
        TypeT T,
        TypeM mass,
        double g,
        unsigned int spin2,
        double abserr = 1e-12,
        int maxiter = 5
) -> decltype(mass * T) {
    return g * 45.0 / (2.0 * M_PI * M_PI) * entropy_density_scaled(mass / T, spin2, abserr, maxiter);
}


}
}

#endif //LANRE_COSMOLOGY_THERMODYNAMIC_FUNCTIONS_HPP
