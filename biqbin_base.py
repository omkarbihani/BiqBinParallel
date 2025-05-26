from abc import abstractmethod

import numpy as np


from solver import (run, set_heuristic, 
                    default_heuristic, get_rank, set_read_data, default_read_data) # We need this to be exposed for import from biqbin_base.py
        


class Solver:
    def __init__(self, problem_instance_file_name, params):
        self.problem_instance_file_name = problem_instance_file_name
        self.params = params
        set_read_data(self.read_data)
        set_heuristic(self.heuristic)

    def read_data(self, fname):
        result = default_read_data(fname)
        return result

    @abstractmethod
    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        ...

    def run(self):
        return run('biqbin', self.problem_instance_file_name, self.params)
    
    def get_rank(self) -> int:
        return get_rank()
