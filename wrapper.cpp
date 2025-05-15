
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/module.hpp>
#include <boost/python/call.hpp>

#include "biqbin.h"

using namespace boost::python;

namespace py = boost::python;
namespace np = boost::python::numpy;

int run_py(py::list py_args)
{
    int argc = len(py_args);
    char **argv = new char *[argc];

    for (int i = 0; i < argc; ++i)
    {
        std::string arg = py::extract<std::string>(py_args[i]);
        argv[i] = strdup(arg.c_str()); // Be sure to free this memory later
    }

    int result = solver_main(argc, argv);

    // Clean up
    for (int i = 0; i < argc; ++i)
    {
        free(argv[i]);
    }
    delete[] argv;
    return result;
}

// Global Python override function (if any)
py::object python_GW_override;

// Setter from Python
void set_GW_heuristic_override(py::object func)
{
    python_GW_override = func;
}

double heuristic(int org_problem_size, double *subproblem_L, int subproblem_n,
                            int *xfixed, int *sol_x, int *x, int num)
{
    if (python_GW_override && !python_GW_override.is_none())
    {
        np::initialize();

        // Wrap C++ pointers as NumPy arrays **without copying**
        np::ndarray np_subL = np::from_data(subproblem_L, np::dtype::get_builtin<double>(),
                                            py::make_tuple(subproblem_n),
                                            py::make_tuple(sizeof(double)), py::object());

        np::ndarray np_xfixed = np::from_data(xfixed, np::dtype::get_builtin<int>(),
                                              py::make_tuple(subproblem_n),
                                              py::make_tuple(sizeof(int)), py::object());

        np::ndarray np_solx = np::from_data(sol_x, np::dtype::get_builtin<int>(),
                                            py::make_tuple(subproblem_n),
                                            py::make_tuple(sizeof(int)), py::object());

        np::ndarray np_x = np::from_data(x, np::dtype::get_builtin<int>(),
                                         py::make_tuple(org_problem_size),
                                         py::make_tuple(sizeof(int)), py::object());

        // Call Python version
        double result = py::extract<double>(python_GW_override(
            org_problem_size, np_subL, subproblem_n, np_xfixed, np_solx, np_x, num));

        // No need to copy back — memory is shared
        return result;
    }

    // Default C++ fallback
    return GW_heuristic(org_problem_size, subproblem_L, subproblem_n,
                        xfixed, sol_x, x, num);
}

double GW_wrapped(int org_problem_size, np::ndarray &subproblem_L, int subproblem_n, np::ndarray &xfixed, np::ndarray &sol_x, np::ndarray &x, int num) {
    // Convert numpy arrays to raw pointers with proper casting
    double* subproblem_L_ptr = reinterpret_cast<double*>(subproblem_L.get_data());
    int* xfixed_ptr = reinterpret_cast<int*>(xfixed.get_data());
    int* sol_x_ptr = reinterpret_cast<int*>(sol_x.get_data());
    int* x_ptr = reinterpret_cast<int*>(x.get_data());

    return GW_heuristic(org_problem_size, subproblem_L_ptr, subproblem_n,
        xfixed_ptr, sol_x_ptr, x_ptr, num);
}

// Python module exposure
BOOST_PYTHON_MODULE(solver)
{
    np::initialize();
    py::def("set_heuristic", set_GW_heuristic_override);
    def("GW_heuristic", &GW_wrapped);
    def("run", &run_py);
}
