import numpy as np

import solver
from solver import run, set_heuristic, set_read_data, read_bqp_data

path = 'rudy/'
argv = ['biqbin', path + 'g05_60.4', 'params']


# def my_heuristic(org_problem_size, subproblem_L, subproblem_n, xfixed, sol_x, x, num):
#     # Modify x and sol_x — this affects the original C++ memory!
#     # solution = GW_heuristic(org_problem_size, subproblem_L, subproblem_n, xfixed, sol_x, x, num)
#     solution = 5000
#     print(f"Python heur {org_problem_size}; {sum(sol_x)}; {sum(x) = }")
#     return solution

# set_read_data(read_bqp_data)
# set_heuristic(my_heuristic)


run(argv)
