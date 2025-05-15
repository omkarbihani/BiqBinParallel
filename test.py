import numpy as np

import solver
from solver import run, GW_heuristic, set_heuristic

path = 'rudy/g05_'
argv = ['biqbin', path + '60.0', 'params']


def my_heuristic(org_problem_size, subproblem_L, subproblem_n, xfixed, sol_x, x, num):
    # Modify x and sol_x — this affects the original C++ memory!
    solution = GW_heuristic(org_problem_size, subproblem_L, subproblem_n, xfixed, sol_x, x, num)
    print(f"Python heur {org_problem_size}; {solution = }; {x = }")
    return solution

set_heuristic(my_heuristic)


run(argv)
