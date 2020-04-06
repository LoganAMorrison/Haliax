#include <lanre/dm_models/higgs_portal.hpp>
#include <lanre/dm_models/kinetic_mixing.hpp>
#include <lanre/dm_models/constant_thermal_cross_section.hpp>
#include <lanre/dm_models/darksun.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace lanre::dm_models;

static std::string hp_init_ds = R"Doc(
Construct a Higgs Portal model.

Parameters
----------
mx: float
    DM mass.
ms: float
    Scalar mediator mass.
gsxx: float
    Coupling of DM to scalar mediator.
smix: float
    Sine of the mixing angle between scalar
    mediator and Higgs.
)Doc";

static std::string hp_scalar_width_ds = R"Doc(
Compute the partial width to a particular final state.

Parameters
----------
state: str
    Final state to compute partial width for. I.e., `state` = 'e e' 
    will compute Gamma(S-> e^+ e^-).
)Doc";

static std::string hp_ann_cs_ds = R"Doc(
Compute the annihilation cross section.

Parameters
----------
Q: float
    Center of mass energy.
state: {str, optional}
    Final state to compute cross section for. I.e., `state`='e e' will
    compute the cross section for x + xbar -> e^+ + e^-.
channel: {str, optional}
    String specifying the restriction on channels. I.e. `channel`='s' 
    will only include s-channel diagrams. Other options are 'tu' 
    (for t,u-channels) and 'all' (for all channels, including 
    interference.)

Returns
-------
sigma: float
    Annihilation cross section for x + xbar -> `state` through the 
    `channel` given a center of mass energy `Q`.
)Doc";

static std::string hp_therm_cs_ds = R"Doc(
Compute the thermally averaged annihilation cross section.

Parameters
----------
x: float
    DM mass divided by its temperature.
state: {str, optional}
    Final state to compute cross section for. I.e., `state`='e e' will
    compute the cross section for x + xbar -> e^+ + e^-.
channel: {str, optional}
    String specifying the restriction on channels. I.e. `channel`='s' 
    will only include s-channel diagrams. Other options are 'tu' 
    (for t,u-channels) and 'all' (for all channels, including 
    interference.)

Returns
-------
<sigma v>: float
    Thermally averaged annihilation cross section for x + xbar -> `state`
    through the `channel` given `x = m / T`.
)Doc";

static std::string hp_rd_ds = R"Doc(
Compute the relic density of the DM.

Parameters
----------
x_start: {float, optional}
    Initial value of the x = mass / T. Default is 1.0.
x_end: {float, optional}
    Final value of the x = mass / T. Default is 100.0.
reltol: {float, optional}
    Relative tolerance used by integrator. Default is 1e-5.
abstol: {float, optional}
    Absolute tolerance used by integrator. Default is 1e-5.

Returns
-------
rd: float
    Relic density of the DM in units Omega_cdm * h^2.
)Doc";


PYBIND11_MODULE(dm_models, m) {
    py::class_<HiggsPortal>(m, "HiggsPortal")
            .def(py::init<double, double, double, double>())
            .def_property("mx", &HiggsPortal::get_mx, &HiggsPortal::set_mx)
            .def_property("ms", &HiggsPortal::get_ms, &HiggsPortal::set_ms)
            .def_property("gsxx", &HiggsPortal::get_gsxx, &HiggsPortal::set_gsxx)
            .def_property("smix", &HiggsPortal::get_smix, &HiggsPortal::set_smix)
            .def_property("width_s", &HiggsPortal::get_width_s, nullptr)
            .def("scalar_partial_width",
                 &HiggsPortal::scalar_partial_width,
                 hp_scalar_width_ds.c_str(),
                 py::arg("state") = "all")
            .def("annihilation_cross_section",
                 &HiggsPortal::annihilation_cross_section,
                 hp_ann_cs_ds.c_str(),
                 py::arg("Q"),
                 py::arg("state") = "all",
                 py::arg("channel") = "all")
            .def("thermal_cross_section",
                 &HiggsPortal::thermal_cross_section,
                 hp_therm_cs_ds.c_str(),
                 py::arg("x"),
                 py::arg("state") = "all",
                 py::arg("channel") = "all")
            .def("relic_density", &HiggsPortal::relic_density,
                 hp_rd_ds.c_str(),
                 py::arg("xstart") = 1.0,
                 py::arg("xend") = 500.0,
                 py::arg("reltol") = 1e-3,
                 py::arg("abstol") = 1e-9,
                 py::arg("alg") = "radau");

    py::class_<KineticMixing>(m, "KineticMixing")
            .def(py::init<double, double, double, double>())
            .def_property("mx", &KineticMixing::get_mx, &KineticMixing::set_mx)
            .def_property("mv", &KineticMixing::get_mv, &KineticMixing::set_mv)
            .def_property("gvxx", &KineticMixing::get_gvxx, &KineticMixing::set_gvxx)
            .def_property("eps", &KineticMixing::get_eps, &KineticMixing::set_eps)
            .def_property("width_v", &KineticMixing::get_width_v, nullptr)
            .def("vector_partial_width",
                 &KineticMixing::vector_mediator_width,
                 "compute width",
                 py::arg("state") = "all")
            .def("annihilation_cross_section",
                 &KineticMixing::annihilation_cross_section,
                 "annihilation cs",
                 py::arg("Q"),
                 py::arg("state") = "all",
                 py::arg("channel") = "all")
            .def("thermal_cross_section",
                 &KineticMixing::thermal_cross_section,
                 "compute tcs",
                 py::arg("x"),
                 py::arg("state") = "all",
                 py::arg("channel") = "all")
            .def("relic_density", &KineticMixing::relic_density,
                 "compute rd",
                 py::arg("xstart") = 1.0,
                 py::arg("xend") = 500.0,
                 py::arg("reltol") = 1e-3,
                 py::arg("abstol") = 1e-9,
                 py::arg("alg") = "rodas")
            .def("solve_boltzmann", &KineticMixing::solve_boltzmann,
                 "Compute the solution to the Boltzmann equation",
                 py::arg("xstart") = 1.0,
                 py::arg("xend") = 500.0,
                 py::arg("reltol") = 1e-3,
                 py::arg("abstol") = 1e-9,
                 py::arg("alg") = "rodas");

    py::class_<ConstantThermalCrossSection>(m, "ConstantThermalCrossSection")
            .def(py::init<double, double, int>())
            .def_property("mx", &ConstantThermalCrossSection::get_mx, &ConstantThermalCrossSection::set_mx)
            .def_property("sigmav", &ConstantThermalCrossSection::get_sigmav, &ConstantThermalCrossSection::set_sigmav)
            .def_property("n", &ConstantThermalCrossSection::get_n, &ConstantThermalCrossSection::set_n)
            .def("relic_density", &ConstantThermalCrossSection::relic_density,
                 "compute rd",
                 py::arg("xstart"),
                 py::arg("xend"),
                 py::arg("reltol"),
                 py::arg("abstol"));

    py::class_<DarkSUN>(m, "DarkSUN")
            .def(py::init<double, unsigned int, double, double, double, double, double, double, bool>(),
                 py::arg("lam"),
                 py::arg("N"),
                 py::arg("L1"),
                 py::arg("c"),
                 py::arg("a"),
                 py::arg("mu_eta"),
                 py::arg("mu_delta"),
                 py::arg("xi_inf"),
                 py::arg("has_dp")
            )
            .def_property("lam", &DarkSUN::get_lam, &DarkSUN::set_lam)
            .def_property("N", &DarkSUN::get_N, &DarkSUN::set_N)
            .def_property("L1", &DarkSUN::get_L1, &DarkSUN::set_L1)
            .def_property("c", &DarkSUN::get_c, &DarkSUN::set_c)
            .def_property("a", &DarkSUN::get_a, &DarkSUN::set_a)
            .def_property("mu_eta", &DarkSUN::get_mu_eta, &DarkSUN::set_mu_eta)
            .def_property("mu_delta", &DarkSUN::get_mu_delta, &DarkSUN::set_mu_delta)
            .def_property("xi_inf", &DarkSUN::get_xi_inf, &DarkSUN::set_xi_inf)
            .def_property("has_dp", &DarkSUN::get_has_dp, &DarkSUN::set_has_dp)
            .def_property("xi_fo", &DarkSUN::get_xi_fo, nullptr)
            .def_property("Tsm_fo", &DarkSUN::get_Tsm_fo, nullptr)
            .def_property("xi_bbn", &DarkSUN::get_xi_bbn, nullptr)
            .def_property("xi_cmb", &DarkSUN::get_xi_cmb, nullptr)
            .def_property("rd_eta", &DarkSUN::get_rd_eta, nullptr)
            .def_property("rd_delta", &DarkSUN::get_rd_delta, nullptr)
            .def_property("dneff_cmb", &DarkSUN::get_dneff_cmb, nullptr)
            .def_property("dneff_bbn", &DarkSUN::get_dneff_bbn, nullptr)
            .def_property("eta_si_per_mass", &DarkSUN::get_eta_si_per_mass, nullptr)
            .def_property("delta_si_per_mass", &DarkSUN::get_delta_si_per_mass, nullptr)
            .def_property("solution", &DarkSUN::get_solution, nullptr)
            .def(py::pickle(
                    [](const DarkSUN &p) {//__getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(
                                p.get_lam(),
                                p.get_N(),
                                p.get_L1(),
                                p.get_c(),
                                p.get_a(),
                                p.get_mu_eta(),
                                p.get_mu_delta(),
                                p.get_xi_inf(),
                                p.get_has_dp(),
                                p.get_xi_fo(),
                                p.get_Tsm_fo(),
                                p.get_xi_bbn(),
                                p.get_xi_cmb(),
                                p.get_rd_eta(),
                                p.get_rd_delta(),
                                p.get_dneff_bbn(),
                                p.get_dneff_cmb(),
                                p.get_eta_si_per_mass(),
                                p.get_delta_si_per_mass(),
                                p.get_solution()
                        );
                    },
                    [](const py::tuple &t) {//__setstate__
                        if (t.size() != 20) {
                            throw std::runtime_error("Invalid state!");
                        }
                        DarkSUN p(
                                t[0].cast<double>(), // lam
                                t[1].cast<unsigned int>(), // N
                                t[2].cast<double>(), // L1
                                t[3].cast<double>(), // c
                                t[4].cast<double>(), // a
                                t[5].cast<double>(), // mu_eta
                                t[6].cast<double>(), // mu_delta
                                t[7].cast<double>(), // xi_inf
                                t[8].cast<bool>() // has_dp
                        );
                        p.set_xi_fo(t[9].cast<double>());
                        p.set_Tsm_fo(t[10].cast<double>());
                        p.set_xi_bbn(t[11].cast<double>());
                        p.set_xi_cmb(t[12].cast<double>());

                        p.set_rd_eta(t[13].cast<double>());
                        p.set_rd_delta(t[14].cast<double>());
                        p.set_dneff_bbn(t[15].cast<double>());
                        p.set_dneff_cmb(t[16].cast<double>());
                        p.set_eta_si_per_mass(t[17].cast<double>());
                        p.set_delta_si_per_mass(t[18].cast<double>());
                        p.set_solution(t[19].cast<ODESolution>());
                        return p;
                    }
            ))
            .def("dark_heff",
                 &DarkSUN::dark_heff,
                 R"Doc(
Compute the effective number of degrees of freedom stored
in entropy in the dark SU(N) model.

Parameters
----------
Td: float
    Dark sector temperature.

Returns
-------
heff: float
    Effective number of d.o.f. in entropy of the dark
    sector.
)Doc",
                 py::arg("Td")
            )
            .def("dark_geff",
                 &DarkSUN::dark_geff,
                 R"Doc(
Compute the effective number of degrees of freedom stored
in energy in the dark SU(N) model.

Parameters
----------
Td: float
    Dark sector temperature.

Returns
-------
geff: float
    Effective number of d.o.f. in energy of the dark
    sector.
)Doc",
                 py::arg("Td")
            )
            .def("sqrt_gstar",
                 &DarkSUN::sqrt_gstar,
                 R"Doc(
Compute the effective value of sqrt(gstar) including
the dark degrees of freedom stored in enetropy.

Parameters
----------
Td: float
    Dark sector temperature.
xi: float
    Ratio of the dark to SM temperatures.

Returns
-------
sqrt_gstar: float
    Effective value of sqrt(gstar).
)Doc",
                 py::arg("Td"),
                 py::arg("xi")
            )
            .def("compute_xi",
                 &DarkSUN::compute_xi,
                 R"Doc(
Compute the ratio of dark to standard model
temperature given a standard model temperature.

Parameters
----------
Tsm: float
    Standard model temperature.

Returns
-------
xi: float
    Ratio of dark to SM temperature.
)Doc",
                 py::arg("Tsm")
            )
            .def("thermal_cross_section_2eta_4eta",
                 &DarkSUN::thermal_cross_section_2eta_4eta,
                 R"Doc(
Compute the thermally averaged cross section for
eta + eta -> eta + eta + eta + eta.

Parameters
----------
x: float
    Mass of the dark eta' divided by its
    temperature.

Returns
-------
<sigmav>: float
    Thermally average cross section for 2eta->4eta.
)Doc",
                 py::arg("x")
            )
            .def("thermal_cross_section_2eta_2delta",
                 &DarkSUN::thermal_cross_section_2eta_2delta,
                 R"Doc(
Compute the thermally averaged cross section for
eta + eta -> delta + delta.

Parameters
----------
x: float
    Mass of the dark eta' divided by its
    temperature.

Returns
-------
<sigmav>: float
    Thermally averaged cross section for 2eta->2delta.
)Doc",
                 py::arg("x")
            )
            .def("cross_section_2eta_2eta",
                 &DarkSUN::cross_section_2eta_2eta,
                 R"Doc(
Compute the self-interation cross section for
eta' + eta' -> eta' + eta'

Returns
-------
sigma: float
    self-interation cross section for 2eta->2eta.
)Doc"
            )
            .def("cross_section_2delta_2delta",
                 &DarkSUN::cross_section_2delta_2delta,
                 R"Doc(
Compute the self-interation cross section for
delta + delta -> delta + delta

Returns
-------
sigma: float
    self-interation cross section for 2delta->2delta.
)Doc"
            )
            .def("delta_n_eff_cmb",
                 &DarkSUN::delta_n_eff_cmb,
                 R"Doc(
Compute the value of dNeff at CMB.

Notes
-----
This should only be ran AFTER 'solve_boltzmann'
or 'relic_densities'.

Returns
-------
dNeff: float
    Delta Neff at CMB.
)Doc"
            )
            .def("delta_n_eff_bbn",
                 &DarkSUN::delta_n_eff_bbn,
                 R"Doc(
Compute the value of dNeff at BBN.

Notes
-----
This should only be ran AFTER 'solve_boltzmann'
or 'relic_densities'.

Returns
-------
dNeff: float
    Delta Neff at BBN.
)Doc"
            )
            .def("solve_boltzmann",
                 &DarkSUN::solve_boltzmann,
                 R"Doc(
Solve the Boltzmann equation

Parameters
----------
reltol: {float, optional}
    Relative tolerance used in ODE integrator.
    Default value is 1e-6.
abstol: {float, optional}
    Absolute tolerance used in ODE integrator.
    Default value is 1e-6.
alg: {str, optional}
    ODE integrator to use. Options are "radau" or
    "rodas". Default is "radau".

Returns
-------
sol: ODESolution
    Solution object. Note that the 'us' represent
    Y = n / s for the eta' and delta respectively.
)Doc",
                 py::arg("reltol") = 1e-6,
                 py::arg("abstol") = 1e-6,
                 py::arg("alg") = "radau"
            );
}