
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


namespace py = boost::python;


/* biqbin's global variables from global_var.h */
extern Problem *SP;
extern Problem *PP;
extern BabSolution *BabSol;
extern int BabPbSize;
//extern double *X;
//extern double *Z; // stores Cholesky decomposition: X = ZZ^T



/* solver functions */
using namespace boost::python;

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace p = boost::python;

/* final solution */
std::vector<int> selected_nodes;
double spent_time;
int rank;

// Global Python override function (if any)
py::object python_heuristic_override;
py::object py_read_data_override;

/// @brief set heuristic function from Python
/// @param func 
void set_heuristic_override(py::object func)
{
    python_heuristic_override = func;
}


/// @brief set instance reading function from Python
/// @param func 
void set_read_data_override(py::object func)
{
    py_read_data_override = func;
}

void set_rank(int r) {
    rank = r;
}

int get_rank() {
    return rank;
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
py::dict run_py(const char * prog_name, const char *problem_instance_name, const char *params_file_name)
{
    const char* argv[3] = {prog_name, problem_instance_name, params_file_name};

    wrapped_main(3, argv);


    // Build result dictionary
    py::dict result_dict;
    result_dict["max_val"] = Bab_LBGet();
    result_dict["solution"] = get_selected_nodes_np_array();
    result_dict["time"] = spent_time;
    return result_dict;
}





double run_heuristic_python(
    const np::ndarray & P0_L_array, 
    const np::ndarray & P_L_array,
    const np::ndarray & xfixed_array,
    const np::ndarray & node_sol_X_array,
    const np::ndarray & x_array)
{
    double * P0_L = reinterpret_cast<double*>(P0_L_array.get_data());
    double * P_L = reinterpret_cast<double*>(P_L_array.get_data());
    int * xfixed = reinterpret_cast<int*>(xfixed_array.get_data());
    int * node_sol_X = reinterpret_cast<int*>(node_sol_X_array.get_data());
    int * x = reinterpret_cast<int*>(x_array.get_data());
 
    return runHeuristic_unpacked(P0_L, P0_L_array.shape(0) , P_L, P_L_array.shape(0), xfixed, node_sol_X, x);
}

/// @brief Called where GW_heuristic was called in the original biqbin, if python_GW_override was not overridden still functions like the original.
/// @param P0 is the original Problem *SP
/// @param P  current subproblem Problem *PP
/// @param node current branch and bound node
/// @param x stores the best solution nodes found the by the heuristic function
/// @param num GW receives this, it is set to SP->n or the number of vertices in the graph
/// @return best lower bound of the current subproblem found by the heuristic used
double wrapped_heuristic(Problem *P0, Problem *P, BabNode *node, int *x)
{

//    return GW_heuristic(P0->L, P0->n, P->L, P->n, node->xfixed, node->sol.X, P0->n, sol);

    np::ndarray P0_L_array = np::from_data(P0->L, np::dtype::get_builtin<double>(),
                                     p::make_tuple(P0->n, P0->n),
                                     p::make_tuple(sizeof(double) * P0->n, sizeof(double)),
                                     p::object());

    np::ndarray P_L_array = np::from_data(P->L, np::dtype::get_builtin<double>(),
                                     p::make_tuple(P->n, P->n),
                                     p::make_tuple(sizeof(double) * P->n, sizeof(double)),
                                     p::object());

    np::ndarray xfixed_array = np::from_data(node->xfixed, np::dtype::get_builtin<int>(),
                                     p::make_tuple(P0->n-1),
                                     p::make_tuple(sizeof(int)),
                                     p::object());
    
    np::ndarray sol_X_array = np::from_data(node->sol.X, np::dtype::get_builtin<int>(),
                                     p::make_tuple(P0->n-1),
                                     p::make_tuple(sizeof(int)),
                                     p::object());

    np::ndarray x_array = np::from_data(x, np::dtype::get_builtin<int>(),
                                     p::make_tuple(BabPbSize),
                                     p::make_tuple(sizeof(int)),
                                     p::object());

//    std::cout << "P" << std::endl;
    return p::extract<double>(python_heuristic_override(P0_L_array, P_L_array, xfixed_array, sol_X_array, x_array));

}


//np::ndarray np_adj;

np::ndarray read_data_python(const char *instance)
{

    double *adj;
    int adj_N;
    adj = readData(instance, &adj_N); // possible mem leak ???
    return np::from_data(adj, np::dtype::get_builtin<double>(),
                                     p::make_tuple(adj_N, adj_N),
                                     p::make_tuple(sizeof(double) * adj_N, sizeof(double)),
                                     p::object());
}

/// @brief Called instead of readData of the original biqbin which could only take edge weight lists in a specific format
/// @param instance path to the instance file
/// @return 0 if parsing was successful 1 if not
double* wrapped_read_data(const char *instance, int *adj_N)
{
    np::ndarray np_adj = py::extract<np::ndarray>(py_read_data_override());
    *adj_N = np_adj.shape(0);

    return reinterpret_cast<double*>(np_adj.get_data());
}

//     free(Adj);



//     if (py_read_data_override && !py_read_data_override.is_none())
//     {
//         // Call Python function and expect a NumPy array
//         np::ndarray np_adj = py::extract<np::ndarray>(py_read_data_override(instance));

//         int n = np_adj.shape(0);
//         BabPbSize = n - 1;

//         double* adj;
//         alloc_matrix(adj, n, double);
//         std::memcpy(adj, np_adj.get_data(), sizeof(double) * n * n);
//         return adj;
//     }
//     return readData(instance);
// }

void clean_python_references(void)
{
    py_read_data_override = py::object(); // Clear the callback
    python_heuristic_override = py::object();
}

// Python module exposure
BOOST_PYTHON_MODULE(solver)
{
    np::initialize();

    py::def("set_heuristic", &set_heuristic_override);
    py::def("set_read_data", &set_read_data_override);
    def("read_bqp_data", &read_data_BQP);
    def("run", &run_py);
    def("default_heuristic", &run_heuristic_python);
    def("default_read_data", &read_data_python);
    def("get_rank", &get_rank);
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