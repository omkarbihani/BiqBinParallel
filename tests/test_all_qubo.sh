#!/bin/bash
# run_all_qubo.sh
dir=$1
nproc=${2:-8}  # default to 8 processes if not passed

if [ -z "$dir" ]; then
  echo "Usage: $0 <directory> [num_processes]"
  exit 1
fi

for file in $(find "$dir" -type f -name "*.json"); do
  PYTHONPATH=.. tests/qubo_test.sh "mpiexec -n $nproc python3 -m tests.run_qubo_test" "$file" params
done