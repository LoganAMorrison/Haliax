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

#include "lanre/cosmology/thermodynamic_functions.hpp"
#include <cmath>

namespace lanre {
namespace cosmology {

class ThermodynamicParticle {
private:
    double m_mass;
    double m_g;
    unsigned int m_spin2;

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
        return cosmology::neq(T, m_mass, m_g, m_spin2);
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
        auto res = cosmology::neq(TT, m_mass, m_g, m_spin2);
        return res.eps;
    }

    /**
     * Compute the equilibrium energy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Energy density
     */
    double energy_density(double T) const {
        return cosmology::energy_density(T, m_mass, m_g, m_spin2);
    }

    /**
     * Compute the equilibrium pressure density of particle with temperature T.
     * @param T Temperature of particle
     * @return Pressure density
     */
    double pressure_density(double T) const {
        return cosmology::pressure_density(T, m_mass, m_g, m_spin2);
    }

    /**
     * Compute the equilibrium entropy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Entropy density
     */
    double entropy_density(double T) const {
        return cosmology::entropy_density(T, m_mass, m_g, m_spin2);
    }

    /**
     * Compute the equilibrium number of d.o.f. stored in energy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in energy
     */
    double geff(double T) const {
        return cosmology::geff(T, m_mass, m_g, m_spin2);
    }

    /**
     * Compute the equilibrium number of d.o.f. stored in entropy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in entropy
     */
    double heff(double T) const {
        return cosmology::heff(T, m_mass, m_g, m_spin2);
    }
};


} // namespace cosmology
} // namespace lanre

#endif //LANRE_COSMOLOGY_THERMODYNAMIC_PARTICLE_HPP
