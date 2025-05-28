#!/bin/bash

output=$($1 $2 $3 2>&1)  # capture both stdout and stderr
status=$?

echo "$output" | sed '1,/Maximum number of workers used:/d'
exit $status