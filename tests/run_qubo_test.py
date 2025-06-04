import sys
import json
from biqbin_base import QUBOSolver, DataGetterJson

"""
    Compares expected output in qubo json file to the computed one
"""
if __name__ == '__main__':
    _, problem_instance_file_name, params = sys.argv

    data_getter = DataGetterJson(problem_instance_file_name)

    solver = QUBOSolver(data_getter=data_getter, params=params)
    result = solver.run()

    rank = solver.get_rank()
    if rank == 0:
        with open(solver.problem_instance_name + "-expected_output.json",  "r") as f:
            expected_result = json.load(f)
        
        test_failed = False
        diffs = {
            "time": expected_result["time"] - result["time"]
        }
        for problem in ["maxcut", "qubo"]:
            com_val = result[problem]["computed_val"]
            exp_val = expected_result[problem]["computed_val"]

            val_diff = com_val - exp_val
            diffs[f"{problem}_val"] = val_diff

            com_x = list(float(i) for i in result[problem]["x"])
            com_len_x = sum(com_x)
            exp_x = expected_result[problem]["x"]
            exp_len_x = sum(exp_x)
            if abs(val_diff) > 0.0001:
                print(f"Computed x = {com_x}")
                print(f"Expected x = {exp_x}")
                test_failed = True
        if not test_failed:
            print(
                f"OK! - {solver.problem_instance_name} Qubo val diff = {diffs["qubo_val"]}; nodes diff = {com_len_x - exp_len_x}; Time diff = {diffs["time"]}"
            )
            exit(0)
        else:
            print(
            f"FAILED! - {solver.problem_instance_name} Maxcut diff = {diffs["maxcut_val"]}; Qubo diff = {diffs["qubo_val"]};  Time diff = {diffs["time"]}"
            )
            exit(1)