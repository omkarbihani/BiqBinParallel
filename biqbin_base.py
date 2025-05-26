from abc import abstractmethod

import numpy as np


from solver import (run, set_heuristic, 
                    default_heuristic) # We need this to be exposed for import from biqbin_base.py
        


class Solver:

    @abstractmethod
    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        ...

    def run(self, argv):
        set_heuristic(self.heuristic)
        return run(argv)