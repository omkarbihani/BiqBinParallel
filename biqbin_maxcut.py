import sys

from biqbin_base import MaxCutSolver


if __name__ == '__main__':

    _, problem_instance_file_name, params = sys.argv

    solver = MaxCutSolver(problem_instance_file_name, params)
    result = solver.run()

    #  We need to expose node ID, so we know if we are on master !!!
    rank = solver.get_rank()
    print(f"{rank=} heuristics ran {solver.heuristic_counter} times")
    if rank == 0:
        print(result)
