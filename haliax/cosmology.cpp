

#include <lanre/cosmology/standard_model.hpp>
#include <lanre/cosmology/thermodynamic_particle.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lanre::cosmology;

static std::string tp_init_ds = R"Doc(
Construct a thermal particle with a given mass, number of internal 
degrees of freedom and spin.

Parameters
----------
mass: float
    Mass of the particle. Units can be anything. But temperature used in 
    function calls must match.
g: float
    Internal degrees of freedom of the particle; i.e. spin, charge, color,
    ect.
spin2: int
    Twice the spin of the particle; i.e. 0 for scalar, 1 for fermion, 2 for
    vector and so on.
)Doc";

static std::string tp_neq_ds = R"Doc(
Compute the equillibrium number density given a temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
neq: float
    Equillibrium number density.
)Doc";

static std::string tp_ed_ds = R"Doc(
Compute the equillibrium energy density given a temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
rho: float
    Equillibrium energy density.
)Doc";

static std::string tp_pd_ds = R"Doc(
Compute the equillibrium pressure density given a temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
p: float
    Equillibrium pressure density.
)Doc";

static std::string tp_entropy_den_ds = R"Doc(
Compute the equillibrium entropy density given a temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
s: float
    Equillibrium entropy density.
)Doc";

static std::string tp_geff_ds = R"Doc(
Compute the effective number of degrees of freedom stored in energy given 
the paticles temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
g: float
    Effective number of d.o.f. in energy.
)Doc";

static std::string tp_heff_ds = R"Doc(
Compute the effective number of degrees of freedom stored in entropy 
given the paticles temperature.

Parameters
----------
T: float
    Temperature of the particle. Units must match the mass of the particle.

Returns
-------
h: float
    Effective number of d.o.f. in entropy.
)Doc";

static std::string sm_geff_ds = R"Doc(
Compute the effective number of degrees of freedom stored in energy in the
SM bath for a given temperature.

Parameters
----------
T: float
    Temperature of the SM bath in GeV.

Returns
-------
g: float
    Effective number of d.o.f. in energy.
)Doc";

static std::string sm_heff_ds = R"Doc(
Compute the effective number of degrees of freedom stored in entropy in the
SM bath for a given temperature.

Parameters
----------
T: float
    Temperature of the SM bath in GeV.

Returns
-------
h: float
    Effective number of d.o.f. in entropy.
)Doc";

static std::string sm_sqrt_gs_ds = R"Doc(
Compute the square root of g-star, which is given by:
    sqrt(gstar) = heff / sqrt(geff) * (1 + T/(3heff) dheff/dT)

Parameters
----------
T: float
    Temperature of the SM bath in GeV.

Returns
-------
sqrt_gstar: float
    Square root of g-star.
)Doc";

static std::string sm_sqrt_energy_dens_ds = R"Doc(
Compute the energy density of the standard model for a given temperature.

Parameters
----------
T: float
    Temperature of the SM bath in GeV.

Returns
-------
rho: float
    Energy density.
)Doc";

static std::string sm_sqrt_entropy_dens_ds = R"Doc(
Compute the entropy density of the standard model for a given temperature.

Parameters
----------
T: float
    Temperature of the SM bath in GeV.

Returns
-------
s: float
    Entropy density.
)Doc";

PYBIND11_MODULE(cosmology, m)
{
    py::class_<ThermodynamicParticle>(m, "ThermodynamicParticle")
        .def(py::init<double, double, unsigned int>())
        .def_property("mass", &ThermodynamicParticle::get_mass, &ThermodynamicParticle::set_mass)
        .def_property("g", &ThermodynamicParticle::get_g, &ThermodynamicParticle::set_g)
        .def_property("spin2", &ThermodynamicParticle::get_spin2, &ThermodynamicParticle::set_spin2)
        .def("neq", &ThermodynamicParticle::neq, tp_neq_ds.c_str(), py::arg("T"))
        .def("energy_density", &ThermodynamicParticle::energy_density, tp_ed_ds.c_str(), py::arg("T"))
        .def("pressure_density", &ThermodynamicParticle::pressure_density, tp_pd_ds.c_str(), py::arg("T"))
        .def("entropy_density", &ThermodynamicParticle::entropy_density, tp_entropy_den_ds.c_str(), py::arg("T"))
        .def("geff", &ThermodynamicParticle::geff, tp_geff_ds.c_str(), py::arg("T"))
        .def("heff", &ThermodynamicParticle::heff, tp_heff_ds.c_str(), py::arg("T"));

    m.def("sm_geff", &sm_geff, sm_geff_ds.c_str(), py::arg("T"));
    m.def("sm_heff", &sm_heff, sm_heff_ds.c_str(), py::arg("T"));
    m.def("sm_sqrt_gstar", &sm_sqrt_gstar, sm_sqrt_gs_ds.c_str(), py::arg("T"));
    m.def("sm_energy_density", &sm_energy_density, sm_sqrt_energy_dens_ds.c_str(), py::arg("T"));
    m.def("sm_entropy_density", &sm_entropy_density, sm_sqrt_entropy_dens_ds.c_str(), py::arg("T"));
}
