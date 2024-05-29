#!/bin/bash -e

# Define the array of values
values=(40 60 80 90 110 120)

# Define the function that will run the command
command_function() {
    local value=$1
    python3 ../match_index.py df_crit.csv subset_${value}.txt subset_${value}.csv /dev/null 
    head -n 1 subset_${value}.txt | python3 -c 'import sys; import re; print(f"k={sys.argv[1]}", re.search(r"discrepancy=([\d\.]+)", input())[0])' ${value}
}

# Run the commands sequentially
for value in "${values[@]}"; do
    command_function $value
done
