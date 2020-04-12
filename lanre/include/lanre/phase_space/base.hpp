//
// Created by Logan Morrison on 3/24/20.
//

#ifndef LANRE_PHASE_SPACE_BASE_HPP
#define LANRE_PHASE_SPACE_BASE_HPP

#include <array>
#include <random>
#include <iostream>

namespace lanre {
namespace phase_space {

struct FourVector {
    double t;
    double x;
    double y;
    double z;
};

FourVector operator+(const FourVector &fv1, const FourVector &fv2) {
    return FourVector{
            fv1.t + fv2.t,
            fv1.x + fv2.x,
            fv1.y + fv2.y,
            fv1.z + fv2.z
    };
}

FourVector operator-(const FourVector &fv1, const FourVector &fv2) {
    return FourVector{
            fv1.t - fv2.t,
            fv1.x - fv2.x,
            fv1.y - fv2.y,
            fv1.z - fv2.z
    };
}

std::ostream &operator<<(std::ostream &os, const FourVector &fv) {
    os << "FourVector(" << fv.t << ", " << fv.x << ", " << fv.y << ", " << fv.z << ")";
    return os;
}

double mass(FourVector &fv) {
    double m = fv.t * fv.t - fv.x * fv.x - fv.y * fv.y - fv.z * fv.z;
    return m >= 0 ? std::sqrt(m) : -std::sqrt(-m);
}


double scalar_product(const FourVector &fv1, const FourVector &fv2) {
    return fv1.t * fv2.t - fv1.x * fv2.x - fv1.y * fv2.y - fv1.z * fv2.z;
}


struct PhaseSpaceEvent {
    std::vector<FourVector> momenta;
    double weight;
};

class PhaseSpaceGenerator {
protected:
    const size_t phase_space_dim;
    /* mutex lock for locking access to m_events */
    std::mutex m_mtx;

    /* common weight factor to all events */
    double m_base_weight{};

    /* private storage container for the events produced by generate_events */
    std::vector<PhaseSpaceEvent> m_events;

    /**
     * Uniform random number generator that is thread-safe.
     * @return random number between (0,1)
     */
    static double phase_space_uniform_rand() {
        static thread_local std::random_device rd{};
        static thread_local std::mt19937 generator{rd()};
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(generator);
    }

public:
    std::vector<double> isp_masses{};
    std::vector<double> fsp_masses{};
    double cme{};
    std::function<double(const std::vector<FourVector> &)> mat_squared;

    // Full constructor
    PhaseSpaceGenerator(std::vector<double> isp_masses, std::vector<double> fsp_masses, double cme,
                        std::function<double(const std::vector<FourVector> &mat_squared)> mat_squared)
            : phase_space_dim(3 * fsp_masses.size() - 4), isp_masses(std::move(isp_masses)),
              fsp_masses(std::move(fsp_masses)), cme(cme), mat_squared(std::move(mat_squared)) {}

    // Constructor assuming constant matrix element
    PhaseSpaceGenerator(std::vector<double> isp_masses, std::vector<double> fsp_masses, double cme)
            : phase_space_dim(3 * fsp_masses.size() - 4),
              isp_masses(std::move(isp_masses)), fsp_masses(std::move(fsp_masses)),
              cme(cme), mat_squared([](const std::vector<FourVector> &) { return 1.0; }) {}

    virtual ~PhaseSpaceGenerator() = default;

    /**
     * Generate several events.
     * @return Vector of PhaseSpaceEvent's containing four-momenta and weight.
     */
    virtual std::vector<PhaseSpaceEvent> generate_events(size_t num_events) = 0;

    std::pair<double, double> compute_width_cross_section(size_t num_events);
};

/**
 * Compute the decay width or scattering cross-section.
 * @param num_events number of events to generate.
 * @return average and standard-deviation.
 */
std::pair<double, double> PhaseSpaceGenerator::compute_width_cross_section(size_t num_events) {
    generate_events(num_events);
    auto num_events_d = (double) m_events.size();

    // Compute average: <w_i> and average of squares: <w_i^2>
    double avg = 0.0, avg2 = 0.0;
    for (auto &event: m_events) {
        double weight = event.weight;
        avg += weight;
        avg2 += weight * weight;
    }
    avg /= num_events_d;
    avg2 /= num_events_d;

    /* Compute the pre-factor of width or cross-section based on the number
     * of initial state particles.
     */
    double pre_factor;
    if (isp_masses.size() == 2) {
        double m1 = isp_masses[0];
        double m2 = isp_masses[1];
        double eng1 = (cme * cme + m1 * m1 - m2 * m2) / (2.0 * cme);
        double eng2 = (cme * cme - m1 * m1 + m2 * m2) / (2.0 * cme);
        double p = sqrt((m1 - m2 - cme) * (m1 + m2 - cme) * (m1 - m2 + cme) *
                                (m1 + m2 + cme)) / (2.0 * cme);

        double v1 = p / eng1, v2 = p / eng2;
        double v_rel = v1 + v2;

        pre_factor = 1.0 / (2.0 * eng1 * 2.0 * eng2 * v_rel);

    } else {
        pre_factor = 1.0 / (2.0 * cme);
    }

    /* Compute standard deviation:
     *  var = <x^2> - <x>^2
     *  sig = sqrt(var / N)
     */
    double var = avg2 - avg * avg;
    double sig = sqrt(var / num_events_d);
    if (isnan(sig)) sig = avg * 1e-12;
    return std::make_pair(pre_factor * avg, pre_factor * sig);
}


}
}

#endif //LANRE_PHASE_SPACE_BASE_HPP