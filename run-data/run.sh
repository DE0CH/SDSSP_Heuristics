#!/bin/bash

# Define the array of values
values=(40 60 80 90 110 120)

# Define the function that will run the command
command_function() {
    local value=$1
    ./a.out df_crit.txt 3 2188 ${value} subset_${value}.txt 2>&1 | ts > log_${value}.txt
}

# Export the function to be used in parallel
export -f command_function

# Run the commands in parallel with 3 cores
parallel -j 2 command_function ::: "${values[@]}"
