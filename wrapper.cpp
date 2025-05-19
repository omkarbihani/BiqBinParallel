
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/module.hpp>
#include <boost/python/call.hpp>

#include "biqbin_cpp_api.h"
#include "wrapper.h"

using namespace boost::python;

namespace py = boost::python;
namespace np = boost::python::numpy;

int run_py(py::list py_args)
{
    int argc = len(py_args);
    // +1 so we can do argv[argc] = nullptr
    char **argv = new char *[argc + 1];

    for (int i = 0; i < argc; ++i) {
        std::string arg = py::extract<std::string>(py_args[i]);
        argv[i] = strdup(arg.c_str());
        // std::cout << argv[i] << std::endl;
    }
    argv[argc] = nullptr;    // <-- sentinel

    int result = wrapped_main(argc, argv);

    for (int i = 0; i < argc; ++i) {
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


double wrapped_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num)
{
        if (python_GW_override && !python_GW_override.is_none())
        {
            np::initialize();
            
        // Wrap C++ pointers as NumPy arrays **without copying**
        np::ndarray np_subL = np::from_data(P->L, np::dtype::get_builtin<double>(),
                                            py::make_tuple(P->n),
                                            py::make_tuple(sizeof(double)), py::object());

        np::ndarray np_xfixed = np::from_data(node->xfixed, np::dtype::get_builtin<int>(),
                                              py::make_tuple(P->n),
                                              py::make_tuple(sizeof(int)), py::object());
                                              
        np::ndarray np_solx = np::from_data(node->sol.X, np::dtype::get_builtin<int>(),
                                            py::make_tuple(P->n),
                                            py::make_tuple(sizeof(int)), py::object());

        np::ndarray np_x = np::from_data(x, np::dtype::get_builtin<int>(),
                                         py::make_tuple(P0->n),
                                         py::make_tuple(sizeof(int)), py::object());

        // Call Python version
        double result = py::extract<double>(python_GW_override(P0->n, np_subL, P->n, np_xfixed, np_solx, np_x, num));
        // No need to copy back — memory is shared
        return result;
    }

    // Default C++ fallback
    return GW_heuristic(P0, P, node, x, num);
}

py::object py_read_data_override;

void set_read_data_override(py::object func) {
    py_read_data_override = func;
}

int wrapped_read_data(const char *instance) {
    if (py_read_data_override && !py_read_data_override.is_none())
    {
        int success = py::extract<int>(py_read_data_override(instance));
        std::cout << "success ?= " << success << std::endl;
        return success;
    }
    return readData(instance);
}

// Python module exposure
BOOST_PYTHON_MODULE(solver)
{
    np::initialize();
    py::def("set_heuristic", set_GW_heuristic_override);
    py::def("set_read_data", set_read_data_override);
    def("read_bqp_data", &read_data_BQP);
    def("run", &run_py);
}