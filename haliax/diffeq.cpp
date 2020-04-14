

#include <lanre/lanre.hpp>
#include <lanre/diffeq/solution.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_MODULE(diffeq, m) {
    py::enum_<lanre::diffeq::Retcode>(m, "Retcode")
            .value("Default", lanre::diffeq::Retcode::Default)
            .value("Success", lanre::diffeq::Retcode::Success)
            .value("MaxIters", lanre::diffeq::Retcode::MaxIters)
            .value("Unstable", lanre::diffeq::Retcode::Unstable)
            .value("DtLessThanMin", lanre::diffeq::Retcode::DtLessThanMin)
            .value("LinearAlgError", lanre::diffeq::Retcode::LinearAlgError)
            .value("SingularMatrix", lanre::diffeq::Retcode::SingularMatrix)
            .value("Failure", lanre::diffeq::Retcode::Failure)
            .export_values();

    py::class_<lanre::diffeq::ODESolution>(m, "ODESolution")
            .def(py::init<>())
            .def_readonly("t", &lanre::diffeq::ODESolution::ts)
            .def_readonly("u", &lanre::diffeq::ODESolution::us)
            .def_readonly("retcode", &lanre::diffeq::ODESolution::retcode)
            .def(py::pickle(
                    [](const lanre::diffeq::ODESolution &p) {//__getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.ts, p.us, p.retcode);
                    },
                    [](const py::tuple &t) {//__setstate__
                        if (t.size() != 3) {
                            throw std::runtime_error("Invalid state!");
                        }
                        lanre::diffeq::ODESolution p{};
                        p.ts = t[0].cast<std::vector<double>>();
                        p.us = t[1].cast<std::vector<lanre::Vector<double>>>();
                        p.retcode = t[2].cast<lanre::diffeq::Retcode>();
                        return p;
                    }
            ));
}
