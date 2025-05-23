import sys
import numpy as np
from neal import SimulatedAnnealingSampler

from solver import run, set_heuristic


sampler = SimulatedAnnealingSampler()


def get_max_cut_solution_SA(L: np.ndarray):
    """Example of a Python heuristic

    Args:
        L (np.ndarray): subproblem L matrix (Problem *PP->L in biqbin)

    Returns:
        np.ndarray: solution nodes provided by the heuristc, should be in 0, 1 form (1 node is chosen, 0 it is not chosen)
    """
    sampleset = sampler.sample_qubo(-L, num_reads=100)

    heuristic_solution = np.array(
        list(sampleset.first.sample.values()), dtype=np.int32)
    return heuristic_solution


# Need to set the new heuristic function
set_heuristic(get_max_cut_solution_SA)

# run needs a Python list ["biqbin", "path/to/instance", "path/to/params"]
# returns a dictionary of the solution nodes ["solution"], maximum value ["max_val"], time taken ["time"]
data = run(sys.argv)

# Only the master process receives the solution
if len(data["solution"]) > 0:
    print(data)
