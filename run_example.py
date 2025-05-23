import sys

from solver import run

# run needs a Python list ["biqbin", "path/to/instance", "path/to/params"]
# returns a dictionary of the solution nodes ["solution"], maximum value ["max_val"], time taken ["time"]
data = run(sys.argv)

# Only the master process receives the solution
if len(data["solution"]) > 0:
    print(data)
