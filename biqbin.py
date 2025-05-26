import sys
import numpy as np

from biqbin_base import Solver, default_heuristic


class BiqbinSolver(Solver):

    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):

        #print(xfixed.sum())
        #print(L.shape)
        #print(L)
        #print(L0)
        #print(x)

        result = default_heuristic(L0, L, xfixed, sol_X, x)

        #print(x)

        # exit(1)

        return result


if __name__ == '__main__':
    solver = BiqbinSolver()
    result = solver.run(sys.argv)
    
    #  We need to expose node ID, so we know if we are on master !!!
    print(result)
