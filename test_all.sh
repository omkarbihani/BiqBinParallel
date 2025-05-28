#!/bin/bash
# run_all_qubo.sh

SCRIPT_NAME=$1        # e.g. "python biqbin_maxcut.py"
SIZE=$2               # e.g. 60
NPROC=${3:-8}         # Default to 8 if not given

if [[ -z "$SCRIPT_NAME" || -z "$SIZE" ]]; then
    echo "Usage: $0 <script_name> <problem_size> [num_procs]"
    exit 1
fi

for i in $(seq 0 9); do
    INPUT="tests/rudy/g05_${SIZE}.${i}"
    EXPECTED="${INPUT}-expected_output"

    ./test.sh "mpiexec -n $NPROC $SCRIPT_NAME" "$INPUT" "$EXPECTED" params
done