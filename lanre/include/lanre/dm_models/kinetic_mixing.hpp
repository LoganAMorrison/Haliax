//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DM_MODELS_KINETIC_MIXING_HPP
#define LANRE_DM_MODELS_KINETIC_MIXING_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/cross_sections.hpp"
#include "lanre/dm_models/kinetic_mixing/widths.hpp"
#include "lanre/dm_models/kinetic_mixing/thermal_cross_section.hpp"
#include "lanre/dm_models/kinetic_mixing/boltzmann.hpp"
#include "lanre/dm_models/kinetic_mixing/relic_density_gondolo_gelmini.hpp"
#include "lanre/dm_models/kinetic_mixing/relic_density_mpu.hpp"
#include "lanre/dm_models/kinetic_mixing/relic_density_bender.hpp"
#include <string>
#

namespace lanre {
namespace dm_models {

using namespace diffeq;
using namespace cosmology;

class KineticMixing {
private:
    kinetic_mixing::Parameters m_params;

public:
    KineticMixing(double mx, double mv, double gvxx, double eps)
            : m_params{mx, mv, gvxx, eps} {
        m_params.widthv = kinetic_mixing::vector_mediator_width(m_params, "all");
    }

    ~KineticMixing() = default;

    double get_mx() const { return m_params.mx; }

    double get_mv() const { return m_params.mv; }

    double get_gvxx() const { return m_params.gvxx; }

    double get_eps() const { return m_params.eps; }

    double get_width_v() const { return m_params.widthv; }

    void set_mx(double mx) {
        m_params.mx = mx;
        m_params.widthv = kinetic_mixing::vector_mediator_width(m_params, "all");
    }

    void set_mv(double mv) {
        m_params.mv = mv;
        m_params.widthv = kinetic_mixing::vector_mediator_width(m_params, "all");
    }

    void set_gvxx(double gvxx) {
        m_params.gvxx = gvxx;
        m_params.widthv = kinetic_mixing::vector_mediator_width(m_params, "all");
    }

    void set_eps(double eps) {
        m_params.eps = eps;
        m_params.widthv = kinetic_mixing::vector_mediator_width(m_params, "all");
    }

    double vector_mediator_width(const std::string &state) const {
        return kinetic_mixing::vector_mediator_width(m_params, state);
    }

    double annihilation_cross_section(double Q, const std::string &state, const std::string &channel) const {
        return kinetic_mixing::annihilation_cross_section(m_params, Q, state, channel);
    }

    double thermal_cross_section(double x, const std::string &state, const std::string &channel) const {
        return kinetic_mixing::thermal_cross_section(m_params, x, state, channel);
    }

    ODESolution solve_boltzmann(double xstart, double xend, double reltol, double abstol, const std::string &t_alg) {
        return kinetic_mixing::solve_boltzmann(m_params, xstart, xend, reltol, abstol, t_alg);
    }

    double relic_density(double xstart, double xend, double reltol, double abstol, const std::string &t_alg) {
        if (t_alg == "gg") {
            return kinetic_mixing::relic_density_gondolo_gelmini(m_params);
        } else if (t_alg == "mpu") {
            return kinetic_mixing::relic_density_mpu(m_params);
        } else if (t_alg == "bender") {
            return kinetic_mixing::relic_density_bender(m_params);
        } else {
            return kinetic_mixing::relic_density(m_params, xstart, xend, reltol, abstol, t_alg);
        }
    }
};

}
}

#endif //LANRE_DM_MODELS_KINETIC_MIXING_HPP
