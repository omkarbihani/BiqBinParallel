from abc import abstractmethod

import numpy as np


from solver import (run, set_heuristic, 
                    default_heuristic, # We need this to be exposed for import from biqbin_base.py
                    get_rank, set_read_data, default_read_data) 
        


class Solver:
    def __init__(self, script_name, problem_instance_file_name, params):
        self.script_name = script_name
        self.problem_instance_file_name = problem_instance_file_name
        self.params = params
        set_read_data(self.read_data)
        set_heuristic(self.heuristic)

    # def convert_sol_from_maxcat_to_qubo(self, ....):
    #     ....

    def read_data(self):
        #qubo = get_qubo_from_somewhere(some_params)
        #adj = convert_to_maxcut(qubo)
        #return adj
        
        result = default_read_data(self.problem_instance_file_name)
        return result

    @abstractmethod
    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        ...

    def run(self):
        result =  run(self.script_name, self.problem_instance_file_name, self.params)
        return result
        #return self.convert_sol_from_maxcat_to_qubo(result)
    
    def get_rank(self) -> int:
        return get_rank()
