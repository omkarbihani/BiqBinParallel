import sys
import json
import numpy as np

from biqbin_base import Solver, default_heuristic, set_read_data


class BiqbinQuboSolver(Solver):
    def __init__(self, script_name, problem_instance_file_name, params):
        self.qubo = None  # needed to convert the maxcut solution to regular solution
        super().__init__(script_name, problem_instance_file_name, params)

    def run(self) -> dict:
        """Runs the original biqbin then adds the qubo solution nodes to the result dict

        Returns:
            dictionary: "solution": maxcut solution nodes, "max_val": maxcut maximum value, "time": solving time, "qubo_solution": qubo solution nodes
        """
        result = super().run()
        if (self.get_rank() == 0):
            result["qubo_solution"] = self.maxcut_solution_to_qubo_solution(
                result["solution"])

        return result

    def heuristic(self, L0: np.ndarray, L: np.ndarray, xfixed: np.array, sol_X: np.array, x: np.array):
        result = default_heuristic(L0, L, xfixed, sol_X, x)
        return result

    def read_data(self) -> np.ndarray:
        """Read qubo json file, return an adjacency matrix for maxcut

        Returns:
            np.ndarray: adjacency matrix
        """
        with open(self.problem_instance_file_name, "r") as f:
            qubo = np.array(json.load(f))
        self.qubo = qubo
        adj = self.qubo_to_maxcut_adj(qubo)
        return adj

    def qubo_to_maxcut_adj(self, q: np.ndarray) -> np.ndarray:
        """Convert general QUBO to adjacency matrix that biqbin can read

        Args:
            q (np.ndarray): qubo numpy array
        Returns:
            np.ndarray: adjacency matrix
        """
        q_sym = 1/2*(q.T + q)
        Qe_plus_c = -np.array([(np.sum(q_sym, 1))])
        np.fill_diagonal(q_sym, 0)

        return np.block([
            [np.zeros((1, 1)), Qe_plus_c],
            [Qe_plus_c.T,     q_sym]
        ])

    def maxcut_solution_to_qubo_solution(self, maxcut_solution: np.ndarray) -> np.ndarray:
        """Convert maxcut solution nodes to qubo solution node

        Args:
            maxcut_solution (np.ndarray): maxcut solution found by biqbin
            maxcut_num_vertices (int): number of vertices in maxcut

        Returns:
            np.ndarray: qubo solution nodes
        """
        _x_mc = np.array(maxcut_solution, dtype=int)-1
        x_mc_sol = np.ones(self.qubo.shape[0] + 1)
        x_mc_sol[_x_mc] = -1
        x_mc_sol *= -x_mc_sol[0]
        x_mc_sol
        y = 1/2*(x_mc_sol+1)[1:]
        qubo_solution = np.nonzero(y)[0] + 1
        return qubo_solution


if __name__ == "__main__":
    script_name, problem_instance_file_name, params = sys.argv
    qubo_solver = BiqbinQuboSolver(
        script_name, problem_instance_file_name, params)

    result = qubo_solver.run()

    # Only the master process returns the solution
    if qubo_solver.get_rank() == 0:
        print("*** PYTHON PRINTS ***")
        for k in result.keys():
            print(f"{k} = {result[k]}")
