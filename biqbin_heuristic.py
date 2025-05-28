import numpy as np
import sys
from neal import SimulatedAnnealingSampler
from biqbin_base import MaxCutSolver
import warnings

"""
*** TESTING ONLY, THIS FILE IS STILL A WORK IN PROGRESS ***
"""
warnings.warn(
    "This functionality is still under development and needs excessive testing!", UserWarning
)


class MaxCutSimulatedAnnealing(MaxCutSolver):
    def __init__(self, problem_instance_name: str, params: str, num_reads: int = 100):
        super().__init__(problem_instance_name, params)
        self.sampler = SimulatedAnnealingSampler()
        self.num_reads = num_reads

    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        """Heuristc with D-Waves simulated annealing sampler

        Args:
            L0 (np.ndarray): main Problem *SP->L matrix
            L (np.ndarray): subproblem *PP->L matrix
            xfixed (np.array): BabNode xfixed variables array
            sol_X (np.array): Solution.X array in BabNode
            x (np.array): stores current best solution

        Returns:
            np.ndarray: solution nodes provided by the heuristc, should be in 0, 1 form (1 node is chosen, 0 it is not chosen)
        """
        self.heuristic_counter += 1
        # Sampler needs -L
        sampler_result = self.sampler.sample_qubo(-L, num_reads=self.num_reads)
        # get the results from sampler
        result = np.array(
            list(sampler_result.first.sample.values()),
            dtype=np.int32
        )
        # copy the results into x, might want to make it return it instead, and we can evaluate it in cpp?
        j = 0
        for i in range(len(x)):
            if xfixed[i] == 0:
                x[i] = result[j]
                j += 1
            else:
                x[i] = sol_X[i]
        sol_value = self.evaluate_solution(L0, x)
        return sol_value

    def evaluate_solution(self, L0: np.ndarray, sol: np.ndarray) -> float:
        """Calculate the lowerbound value of heuristic solution

        Args:
            L0 (np.ndarray): main Problem *SP->L matrix
            sol (np.ndarray): current solution

        Returns:
            float: value of the solution
        """
        sol_val = 0
        for i in range(len(sol)):
            for j in range(len(sol)):
                sol_val += L0[i][j] * sol[i] * sol[j]
        return sol_val


if __name__ == '__main__':

    _, problem_instance_file_name, params = sys.argv

    solver = MaxCutSimulatedAnnealing(problem_instance_file_name, params)
    result = solver.run()

    #  We need to expose node ID, so we know if we are on master !!!
    rank = solver.get_rank()
    print(f"{rank=} heuristics ran {solver.heuristic_counter} times")
    if rank == 0:
        print(result)
