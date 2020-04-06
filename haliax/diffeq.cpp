

#include <lanre/diffeq/solution.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace lanre::diffeq;

PYBIND11_MODULE(diffeq, m) {
    py::enum_<Retcode>(m, "Retcode")
            .value("Default", Retcode::Default)
            .value("Success", Retcode::Success)
            .value("MaxIters", Retcode::MaxIters)
            .value("Unstable", Retcode::Unstable)
            .value("DtLessThanMin", Retcode::DtLessThanMin)
            .value("LinearAlgError", Retcode::LinearAlgError)
            .value("SingularMatrix", Retcode::SingularMatrix)
            .value("Failure", Retcode::Failure)
            .export_values();

    py::class_<ODESolution>(m, "ODESolution")
            .def(py::init<>())
            .def_readonly("t", &ODESolution::ts)
            .def_readonly("u", &ODESolution::us)
            .def_readonly("retcode", &ODESolution::retcode)
            .def(py::pickle(
                    [](const ODESolution &p) {//__getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.ts, p.us, p.retcode);
                    },
                    [](const py::tuple &t) {//__setstate__
                        if (t.size() != 3) {
                            throw std::runtime_error("Invalid state!");
                        }
                        ODESolution p{};
                        p.ts = t[0].cast<std::vector<double>>();
                        p.us = t[1].cast<std::vector<Vector<double>>>();
                        p.retcode = t[2].cast<Retcode>();
                        return p;
                    }
            ));
}
