__version__ = '1.0.0'

from abc import ABC, abstractmethod
import numpy as np
import json

from solver import (run, set_heuristic, 
                    default_heuristic, 
                    get_rank, set_read_data, 
                    default_read_data) 
        


class MaxCutSolver:
    solver_name = f'PyBiqBin-MaxCut {__version__}'

    def __init__(self, problem_instance_name, params):
        self.problem_instance_name = problem_instance_name
        self.params = params
        set_read_data(self.read_data)
        set_heuristic(self.heuristic)

    def read_data(self):
        result = default_read_data(self.problem_instance_name)
        return result

    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        return default_heuristic(L0, L, xfixed, sol_X, x)

    def run(self):
        return run(self.solver_name, self.problem_instance_name, self.params)
    
    def get_rank(self) -> int:
        return get_rank()


class DataGetter(ABC):
    @abstractmethod
    def problem_instance_name(self):
        ...

    @abstractmethod
    def problem_instance(self):
        ...


class DataGetterJson(DataGetter):
    def __init__(self, filename):
         self.filename = filename
         with open(filename, "r") as f:
            self.qubo = np.array(json.load(f))

    def problem_instance_name(self):
        return self.filename

    def problem_instance(self):
        return self.qubo


class QUBOSolver(MaxCutSolver):
    solver_name = f'PyBiqBin-QUBO {__version__}'

    def __init__(self, data_getter, params):
        self.data_getter = data_getter
        super().__init__(data_getter.problem_instance_name(), params)

    def _qubo2maxcut(self, qubo: np.ndarray) -> np.ndarray:
        """Convert QUBO to adjacency matrix that biqbin can read

        Args:
            qubo (np.ndarray): qubo as 2d numpy array
        Returns:
            2d np.ndarray: adjacency matrix of max cut problem
        """
        q_sym = 1/2*(qubo.T + qubo)
        Qe_plus_c = -np.array([(np.sum(q_sym, 1))])
        np.fill_diagonal(q_sym, 0)

        return np.block([
            [np.zeros((1, 1)), Qe_plus_c],
            [Qe_plus_c.T,     q_sym]
        ])

    def _maxcut_solution2qubo_solution(self, maxcut_solution: np.ndarray) -> np.ndarray:
        """Convert maxcut solution nodes to qubo solution node

        Args:
            maxcut_solution (np.ndarray): maxcut solution found by biqbin
            maxcut_num_vertices (int): number of vertices in maxcut

        Returns:
            np.ndarray: qubo solution nodes
        """

        n, _ = self.data_getter.problem_instance().shape

        _x_mc = np.array(maxcut_solution, dtype=int)-1
        x_mc_sol = np.ones(n + 1)
        x_mc_sol[_x_mc] = -1
        x_mc_sol *= -x_mc_sol[0]
        x_mc_sol
        y = 1/2*(x_mc_sol+1)[1:]
        qubo_solution = np.nonzero(y)[0] + 1
        return qubo_solution, y

    def read_data(self) -> np.ndarray:
        """Read qubo json file, return an adjacency matrix for maxcut

        Returns:
            np.ndarray: adjacency matrix
        """
        return self._qubo2maxcut(self.data_getter.problem_instance())

    def run(self) -> dict:
        """Runs the original biqbin then adds the qubo solution nodes to the result dict

        Returns:
            dictionary: "solution": maxcut solution nodes, "max_val": maxcut maximum value, "time": solving time, "qubo_solution": qubo solution nodes
        """
        result = super().run()
        if (self.get_rank() == 0):
            qubo_solution, qubo_x = self._maxcut_solution2qubo_solution(result["solution"])
            qubo_min_value = self.data_getter.problem_instance().dot(qubo_x).dot(qubo_x)
            
            return {'maxcut': result, 'qubo': {'solution': qubo_solution, 'x': qubo_x, 'min_value': float(qubo_min_value)}}
        else:
            return None


