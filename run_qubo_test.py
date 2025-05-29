import sys
from biqbin_base import QUBOSolver, DataGetterJson


if __name__ == '__main__':
    _, problem_instance_file_name, params = sys.argv

    data_getter = DataGetterJson(problem_instance_file_name)

    solver = QUBOSolver(data_getter=data_getter, params=params)
    result = solver.run()

    rank = solver.get_rank()
    if rank == 0:
        qubo_data = solver.data_getter.qubo_data
        expected_optimum = qubo_data["optimum"]
        expected_x = qubo_data["x"]
        expected_len_x = sum(expected_x)

        qubo_result = result["qubo"]
        computed_optimum = qubo_result["computed_val"]
        computed_x = [int(i) for i in qubo_result["x"]]
        computed_len_x = sum(computed_x)

        if computed_optimum == expected_optimum:
            print(
                f"OK! - {qubo_data["name"]} Computed value == expected value = {computed_optimum}; nodes diff = {computed_len_x - expected_len_x}; Time = {result["maxcut"]["time"]}"
            )
            exit(0)
        else:
            print(
                f"FAILED! - {qubo_data["name"]} Computed value: {computed_optimum} != Expected value: {expected_optimum} diff = {computed_optimum - expected_optimum};  Time = {result["maxcut"]["time"]}"
            )
            print(
                f"Computed x num nodes = {computed_len_x} - Expected x num nodes = {expected_len_x}; diff = {computed_len_x - expected_len_x}")
            print(f"Computed x = {computed_x}")
            print(f"Expected x = {expected_x}")
            exit(1)
