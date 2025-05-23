import sys
import numpy as np
import joblib

from qubo_helper_functions import maxcut_to_L, qubo_to_maxcut, maxcut_solution_to_qubo_solution
from solver import run, set_heuristic, set_read_data, read_bqp_data


qubo_info = joblib.load(filename=sys.argv[1])


def read_data(filepath: str) -> np.ndarray:
    """Python function that overrides the original biqbin readData, needs to return a "Laplacean"

    Args:
        filepath (str): function recieves the file path to the instance.

    Returns:
        np.ndarray: L matrix (idk what to call this, it's not Laplacean exactly)
    """
    qubo_info = joblib.load(filename=filepath)

    # print(qubo_info["QUBO"])
    qubo = qubo_info["QUBO"]
    maxcut = qubo_to_maxcut(qubo)
    La = maxcut_to_L(maxcut)
    return La


set_read_data(read_data)  # need to set the function in C++ wrapper

# run needs a Python list ["biqbin", "path/to/instance", "path/to/params"]
# returns a dictionary of the solution nodes ["solution"], maximum value ["max_val"], time taken ["time"]
data = run(sys.argv)

# Only the master process returns the solution
if len(data["solution"]) > 0:
    print(data)

    qubo_solution = maxcut_solution_to_qubo_solution(
        data["solution"],
        len(qubo_info["QUBO"]) + 1
    )

    print("Qubo solution = " +
          str(float(qubo_info["QUBO"]@qubo_solution@qubo_solution)))
