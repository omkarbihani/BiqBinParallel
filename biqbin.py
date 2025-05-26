import sys
import numpy as np

from biqbin_base import Solver, default_heuristic

class BiqbinSolver(Solver):

    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        result = default_heuristic(L0, L, xfixed, sol_X, x)
        return result


if __name__ == '__main__':

    _, problem_instance_file_name, params = sys.argv

    solver = BiqbinSolver(problem_instance_file_name, params)
    result = solver.run()
    
    #  We need to expose node ID, so we know if we are on master !!!
    rank =  solver.get_rank()
    if rank == 0:
        print(f"{rank = }")
        print(result)


