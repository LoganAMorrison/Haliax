/**
 * Created by Logan Morrison on 2/10/20.
 *
 * Class for computing various thermodynamic quantities of a particle in
 * thermal equilibrium.
 *
 * Attributes:
 *      double mass: mass of the particle in GeV
 *      double g: internal dof of the particle
 *      double spin2: 2 times the particles spin (i.e., for spin=1/2, spin2 = 1)
 *
 * Methods:
 *      double neq(double T): equilibrium number density
 *      double energy_density(double T): equilibrium energy density
 *      double pressure_density(double T): equilibrium pressure density
 *      double entropy_density(double T): equilibrium entropy density
 *      double dof_energy(double T): equilibrium dof in energy
 *      double dof_entropy(double T): equilibrium dof in entropy
 *
 */

#ifndef LANRE_COSMOLOGY_THERMODYNAMIC_PARTICLE_HPP
#define LANRE_COSMOLOGY_THERMODYNAMIC_PARTICLE_HPP

#include "lanre/autodiff/dual.hpp"
#include "lanre/special_functions/besselk.hpp"
#include <cmath>

namespace lanre::cosmology {

class ThermodynamicParticle {
private:
    double m_mass;
    double m_g;
    unsigned int m_spin2;

    /**
     * Number density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return neq(T) / (g T^3)
     */
    template<class Type>
    Type neqScaled(Type x) const {
        using lanre::special_functions::besselk2;
        if (x == 0.0)
            return (m_spin2 % 2 == 0) ? static_cast<Type>(121793828233573) : static_cast<Type>(0.0913453711751798);
        else {
            using namespace boost::math;
            Type eta = static_cast<Type>((m_spin2 % 2 == 0) ? 1.0 : -1.0);

            Type sum = static_cast<Type>(0.0);
            for (int n = 1; n <= 5; n++) {
                sum += pow(eta, n + 1) * besselk2(n * x) / n;
            }
            return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
        }
    }

    /**
     * Energy density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return energy_density(T) / (g T^4)
     */
    template<class Type>
    Type energyDensityScaled(Type x) const {
        using lanre::special_functions::besselk1;
        using lanre::special_functions::besselk2;
        if (x == 0.0)
            return static_cast<Type>(M_PI * M_PI / 30.0 * ((m_spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0));
        else {
            Type eta = static_cast<Type>((m_spin2 % 2 == 0) ? 1.0 : -1.0);

            Type sum = 0.0;
            for (int n = 1; n <= 5; n++) {
                sum += pow(eta, n + 1) / (n * n) * (x * n * besselk1(n * x) + 3.0 * besselk2(n * x));
            }
            return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
        }
    }

    /**
     * Pressure density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return pressure_density(T) / (g T^4)
     */
    template<class Type>
    Type pressureDensityScaled(Type x) const {
        using lanre::special_functions::besselk2;
        if (x == 0.0)
            return static_cast<Type>(M_PI * M_PI / 90.0 * ((m_spin2 % 2 == 0) ? 1 : 7.0 / 8.0));
        else {
            Type eta = static_cast<Type>((m_spin2 % 2 == 0) ? 1.0 : -1.0);

            Type sum = static_cast<Type>(0.0);
            for (int n = 1; n <= 5; n++) {
                sum += pow(eta, n + 1) / (n * n) * besselk2(n * x);
            }
            return x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
        }
    }

    /**
     * Entropy density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return entropy_density(T) / (g T^3)
     */
    template<class Type>
    Type entropyDensityScaled(Type x) const {
        using lanre::special_functions::besselk3;
        if (x == static_cast<Type>(0.0))
            return static_cast<Type>(2.0 * M_PI * M_PI / 45.0 * ((m_spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0));
        else {
            Type eta = static_cast<Type>((m_spin2 % 2 == 0) ? 1.0 : -1.0);

            Type sum = 0.0;
            for (int n = 1; n <= 5; n++) {
                sum += pow(eta, n + 1) / n * besselk3(n * x);
            }
            return x * x * x * sum / static_cast<Type>(2 * M_PI * M_PI);
        }
    }

public:
    ThermodynamicParticle(double t_mass, double t_g, unsigned int t_spin2)
            : m_mass(t_mass), m_g(t_g), m_spin2(t_spin2) {}

    ~ThermodynamicParticle() = default;

    double get_mass() const { return m_mass; }

    double get_g() const { return m_g; }

    unsigned int get_spin2() const { return m_spin2; }

    void set_mass(double t_mass) { m_mass = t_mass; }

    void set_g(double t_g) { m_g = t_g; }

    void set_spin2(unsigned int t_spin2) { m_spin2 = t_spin2; }

    /**
     * Compute the equilibrium number density of particle with temperature T.
     * @param T Temperature of particle
     * @return Number density
     */
    double neq(double T) const {
        return m_g * neqScaled(m_mass / T) * T * T * T;
    }

    /**
     * Compute the derivative of the equilibrium number density of a particle
     * with temperature T.
     * @param T Temperature of particle.
     * @return Derivative of the number density.
     */
    double neq_deriv(double T) const {
        using lanre::autodiff::Dual;
        Dual<double> TT{T, 1.0};
        Dual<double> res = m_g * neqScaled(m_mass / TT) * TT * TT * TT;
        return res.eps;
    }

    /**
     * Compute the equilibrium energy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Energy density
     */
    double energy_density(double T) const {
        return m_g * energyDensityScaled(m_mass / T) * T * T * T * T;
    }

    /**
     * Compute the equilibrium pressure density of particle with temperature T.
     * @param T Temperature of particle
     * @return Pressure density
     */
    double pressure_density(double T) const {
        return m_g * pressureDensityScaled(m_mass / T) * T * T * T * T;
    }

    /**
     * Compute the equilibrium entropy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Entropy density
     */
    double entropy_density(double T) const {
        return m_g * entropyDensityScaled(m_mass / T) * T * T * T;
    }

    /**
     * Compute the equilibrium number of d.o.f. stored in energy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in energy
     */
    double geff(double T) const {
        return 30.0 / (M_PI * M_PI) * m_g * energyDensityScaled(m_mass / T);
    }

    /**
     * Compute the equilibrium number of d.o.f. stored in entropy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in entropy
     */
    double heff(double T) const {
        return m_g * 45.0 / (2.0 * M_PI * M_PI) * entropyDensityScaled(m_mass / T);
    }
};

} // namespace lanre::cosmology

#endif //LANRE_COSMOLOGY_THERMODYNAMIC_PARTICLE_HPP
