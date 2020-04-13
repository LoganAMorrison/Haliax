//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_PARTICLES_HPP
#define LANRE_PARTICLES_HPP

#include "lanre/constants.hpp"

namespace lanre {

class Particle {
private:
    // Fundamentals
    double m_mass;
    int m_spin2;
    double m_width = 0.0;
    // Quantum numbers
    int m_col;
    double m_charge = 0.0;
    int m_id = -1;
public:
    Particle(
            double mass,
            int spin2
    ) : m_mass(mass), m_spin2(spin2) {}

    double mass() { return m_mass; }

    int spin2() { return m_spin2; }


};

// Leptons

const Particle electron = Particle{kELECTRON_MASS, 1};

const Particle muon = Particle{kMUON_MASS, 1};

const Particle tau = Particle{kTAU_MASS, 1};

// Quarks

const Particle UpQuark = Particle{kUP_QUARK_MASS, 1};

const Particle CharmQuark = Particle{kCHARM_QUARK_MASS, 1};

const Particle TopQuark = Particle{kTOP_QUARK_MASS, 1};

const Particle DownQuark = Particle{kDOWN_QUARK_MASS, 1};

const Particle StrangeQuark = Particle{kSTRANGE_QUARK_MASS, 1};

const Particle BottomQuark = Particle{kBOTTOM_QUARK_MASS, 1};


}

#endif //LANRE_PARTICLES_HPP
