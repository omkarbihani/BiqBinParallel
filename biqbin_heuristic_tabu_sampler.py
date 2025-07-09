import numpy as np
import sys
from neal import SimulatedAnnealingSampler
from biqbin_base import QUBOSolver, DataGetterJson, default_heuristic
import logging
from copy import deepcopy


import subprocess
import json
import os
from dwave.samplers import TabuSampler


logger = logging.getLogger(__name__)


def custom_sa(Q, max_iter=None, initial_temp=None, final_temp=0.1, alpha=0.999):

    n = Q.shape[0]
    if initial_temp is None:
        initial_temp = n

    x = np.random.randint(0, 2, size=n)

    Q_diag = np.diag(Q)

    energy = x @ Q @ x

    best_solution = x.copy()
    best_energy = energy

    T = initial_temp
    iteration = 0

    while T >= final_temp:
        i = np.random.randint(n)

        # ΔE = (1 - 2*x[i]) * (Q[i,i] + 2 * sum_j != i Q[i,j] * x[j])
        delta_energy = (1 - 2 * x[i]) * (Q_diag[i] + 2 * (Q[i, :] @ x - Q[i, i] * x[i]))

        if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / T):
            # Flip bit i
            x[i] = 1 - x[i]
            energy += delta_energy

            if energy < best_energy:
                best_energy = energy
                best_solution = x.copy()

        T *= alpha
        iteration += 1

        if max_iter is not None and iteration >= max_iter:
            break

    return best_solution.tolist()


class QuboCustomSA(QUBOSolver):

    def heuristic(self, L0, L, xfixed, sol_X, x):
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
        custom_sa_x = custom_sa(-L[:-1, :-1])

        _x = np.array(
            custom_sa_x,
            dtype=np.int32
        )

        #  This can be likely simplified.
        # _x = 2*_x-1                # normalize between -1, 1
        # _x = np.hstack([_x, [-1]]) # append -1
        # _x = 1/2*(_x+1)            # normalize back to 0, 1

        _x = np.hstack([_x, [0]]) # simplification for above

        j = 0
        for i in range(len(x)):
            if xfixed[i] == 0:
                x[i] = _x[j]
                j += 1
            else:
                x[i] = sol_X[i]

        sol_value = self.evaluate_solution(L0, x)

        if logger.isEnabledFor(logging.DEBUG):
            her_value = default_heuristic(L0, L, xfixed, sol_X, deepcopy(x))
            logger.debug(
                f'Custom heuristic: {sol_value}, default heuristic: {her_value}')

        return sol_value

class QUBOTabuSearch(QUBOSolver):

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
        self.sampler = TabuSampler()

        _x = np.array(
            list(self.sampler.sample_qubo(-L[:-1, :-1],
                 num_reads=self.num_reads).first.sample.values()),
            dtype=np.int32
        )


        #  This can be likely simplified.
        # _x = 2*_x-1                # normalize between -1, 1
        # _x = np.hstack([_x, [-1]]) # append -1
        # _x = 1/2*(_x+1)            # normalize back to 0, 1

        _x = np.hstack([_x, [0]]) # simplification for above

        j = 0
        for i in range(len(x)):
            if xfixed[i] == 0:
                x[i] = _x[j]
                j += 1
            else:
                x[i] = sol_X[i]

        sol_value = self.evaluate_solution(L0, x)

        if logger.isEnabledFor(logging.DEBUG):
            her_value = default_heuristic(L0, L, xfixed, sol_X, deepcopy(x))
            logger.debug(
                f'Custom heuristic: {sol_value}, default heuristic: {her_value}')

        return sol_value

if __name__ == '__main__':

    # https://stackoverflow.com/questions/7016056/python-logging-not-outputting-anything
    logging.basicConfig()
    logging.root.setLevel(logging.DEBUG)

    problem_instance_file_name = "qubo_example_2.json"
    params = "params"

    data_getter = DataGetterJson(problem_instance_file_name)

    # solver = QuboQBSolve(data_getter, params=params, num_reads=5)
    # solver = QuboCustomSA(data_getter, params=params, num_reads=5)
    solver = QUBOTabuSearch(data_getter, params=params, num_reads=5)
    result = solver.run()

    rank = solver.get_rank()
    print(f"{rank=} heuristics ran {solver.heuristic_counter} times")
    if rank == 0:
        # Master rank prints the results
        print(result)
        solver.save_result(result)