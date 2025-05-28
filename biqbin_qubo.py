import sys

from biqbin_base import QUBOSolver, DataGetterJson


if __name__ == '__main__':

    _, problem_instance_file_name, params = sys.argv

    data_getter = DataGetterJson(problem_instance_file_name)

    solver = QUBOSolver(data_getter=data_getter, params=params)
    result = solver.run()

    rank = solver.get_rank()
    if rank == 0:
        print(f"{rank=}")
        print(result)
        print(f"Qubo Minimal Value = {round(result["qubo"]["min_value"], 4)}")
