
#include <iostream>

/* Boost library, some are unnecessary */
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/module.hpp>
#include <boost/python/call.hpp>

#include "biqbin_cpp_api.h"
#include "blas_laplack.h"
#include "wrapper.h"

/* biqbin's global variables from global_var.h */
extern Problem *SP;
extern Problem *PP;
extern BabSolution *BabSol;
extern int BabPbSize;
extern double *X;
extern double *Z; // stores Cholesky decomposition: X = ZZ^T



/* solver functions */
using namespace boost::python;

namespace py = boost::python;
namespace np = boost::python::numpy;

/* final solution */
std::vector<int> selected_nodes;
double spent_time;

// Global Python override function (if any)
py::object python_GW_override;
py::object py_read_data_override;

/// @brief set heuristic function from Python
/// @param func 
void set_GW_heuristic_override(py::object func)
{
    python_GW_override = func;
}
/// @brief set instance reading function from Python
/// @param func 
void set_read_data_override(py::object func)
{
    py_read_data_override = func;
}

/// @brief Creates a numpy array of the solution, returned after biqbin is done solving
/// @return np.ndarray(dtype = np.int32) of the final solution (node names in a np list)
np::ndarray get_selected_nodes_np_array()
{
    // Create NumPy array to return to user
    np::dtype dtype = np::dtype::get_builtin<int>();
    py::tuple shape = py::make_tuple(selected_nodes.size());

    // Allocate NumPy array
    np::ndarray result = np::zeros(shape, dtype);
    for (int i = 0; i < static_cast<int>(selected_nodes.size()); ++i)
    {
        result[i] = selected_nodes[i];
    }
    return result;
}

/// @brief Run the solver, retrieve the solution
/// @param py_args argv in a Python string list ["biqbin", "graph_instance_path", "parameters_path"]
/// @return dictionary of "max_val" value of maximum cut and "solution" vertices
py::dict run_py(py::list py_args)
{
    int argc = len(py_args);
    // + 1 so we can do argv[argc] = nullptr
    char **argv = new char *[argc + 1];

    for (int i = 0; i < argc; ++i)
    {
        std::string arg = py::extract<std::string>(py_args[i]);
        argv[i] = strdup(arg.c_str());
    }
    argv[argc] = nullptr; // <-- sentinel

    wrapped_main(argc, argv);

    // Free argv from memory
    for (int i = 0; i < argc; ++i)
    {
        free(argv[i]);
    }
    delete[] argv;

    // Build result dictionary
    py::dict result_dict;
    result_dict["max_val"] = Bab_LBGet();
    result_dict["solution"] = get_selected_nodes_np_array();
    result_dict["time"] = spent_time;
    return result_dict;
}




/// @brief GW_heuristic contains some necessary execution that needs to be run regardless of used heuristic, 
/// should probably be copied into run_heuristic function rather
void pre_heuristic_computes()
{
    /* GW necesities for upper bound computation */
    int n = PP->n;
    int inc = 1;
    char UPLO = 'L';
    int info = 0;
    int nn = n * n;

    // Z = X
    dcopy_(&nn, X, &inc, Z, &inc);

    // compute Cholesky factorization
    dpotrf_(&UPLO, &n, Z, &n, &info);

    // set lower triangle of Z to zero
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            Z[j + i * n] = 0.0;
        }
    }
}

double wrapped_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num)
{
    if (python_GW_override && !python_GW_override.is_none())
    {
        pre_heuristic_computes();
        std::vector<int> sol(P0->n - 1);

        double best_value = - 1e+9; // best lower bound found

        /* convert to numpy arrays */
        // Subproblems L matrix
        np::ndarray np_subL = np::from_data(
            P->L, np::dtype::get_builtin<double>(),
            py::make_tuple(P->n, P->n),
            py::make_tuple(sizeof(double) * P->n, sizeof(double)),
            py::object());

        // Call Python heuristics, heuristic solution should be stored in np_temp_x!
        np::ndarray np_temp_x = boost::python::extract<np::ndarray>(python_GW_override(np_subL));

        // Bellow is taken from GW as well, it works the same as GW
        // store local cut temp_x into global cut sol
        int index = 0;
        for (int i = 0; i < P0->n - 1; ++i)
        {
            if (node->xfixed[i]) {
                sol[i] = node->sol.X[i];
            }
            else
            {
                int value = boost::python::extract<int>(np_temp_x[index]);
                sol[i] = (value + 1) / 2;
                ++index;
            }
        }
        // Update the x solution from temp solution sol, if we found better
        update_best(x, sol.data(), &best_value, P0);

        return best_value;
    }

    // Default C++ fallback
    return GW_heuristic(P0, P, node, x, num);
}

int wrapped_read_data(const char *instance)
{
    if (py_read_data_override && !py_read_data_override.is_none())
    {
        // py_read_data_override(instance);
        np::initialize();

        // Call Python function and expect a NumPy array
        np::ndarray arr = py::extract<np::ndarray>(py_read_data_override(instance));

        int n = arr.shape(0);
        int m = arr.shape(1);

        alloc(SP, Problem);
        alloc(PP, Problem);

        SP->n = n;
        PP->n = n;
        BabPbSize = n - 1;
        // allocate memory for objective matrices for SP and PP
        alloc_matrix(SP->L, SP->n, double);
        alloc_matrix(PP->L, SP->n, double);
        double *numpy_data = reinterpret_cast<double *>(arr.get_data());
        std::memcpy(SP->L, numpy_data, n * m * sizeof(double));
        std::memcpy(PP->L, numpy_data, n * m * sizeof(double));

        return 0;
    }
    return readData(instance);
}

void clean_python_references(void)
{
    py_read_data_override = py::object(); // Clear the callback
    python_GW_override = py::object();
}

// Python module exposure
BOOST_PYTHON_MODULE(solver)
{
    Py_Initialize();
    np::initialize();

    py::def("set_heuristic", &set_GW_heuristic_override);
    py::def("set_read_data", &set_read_data_override);
    def("read_bqp_data", &read_data_BQP);
    def("run", &run_py);
}

/// @brief Copy the solution before memory is freed, so it can be retrieved in Python
void copy_solution() {
    for (int i = 0; i < BabPbSize; ++i)
    {
        if (BabSol->X[i] == 1)
        {
            selected_nodes.push_back(i + 1); // 1-based indexing
        }
    }
}

/// @brief record time at the end
/// @param time_taken 
void record_time(double time_taken) {
    spent_time = time_taken;
}