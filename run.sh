#!/bin/bash -e

# Function to show usage
usage() {
    echo "Usage: $0 executable input_file dimension num_points values [-j num_cores]"
    echo "  executable   The executable to run (e.g., ./a.out)"
    echo "  input_file   Input file name (e.g., df_crit.txt)"
    echo "  dimension    Dimension value (e.g., 3)"
    echo "  num_points   Number of points in the input file (e.g., 2188)"
    echo "  values       A comma-separated list of values (number of points in the subset) (e.g., 40,60,80,90)"
    echo "  -j num_cores (Optional) Number of cores to use in parallel (default is 1)"
    exit 1
}

# Ensure at least 5 positional arguments are provided
if [ "$#" -lt 5 ]; then
    usage
fi

# Default number of cores
num_cores=1

# Parse positional arguments
executable="$1"
input_file="$2"
dimension="$3"
num_points="$4"
IFS=',' read -r -a values <<< "$5"
shift 5

# Parse optional arguments
while getopts "j:" opt; do
    case $opt in
        j)
            num_cores="$OPTARG"
            ;;
        *)
            usage
            ;;
    esac
done

# Define the function that will run the command
command_function() {
    local value=$1
    ${executable} ${input_file} ${dimension} ${num_points} ${value} subset_${value}.txt 2>&1 | ts > log_${value}.txt
}


# Export the function if parallel execution is used
if [ "$num_cores" -gt 1 ]; then
    export -f command_function
    # Run the commands in parallel with the specified number of cores
    parallel -j "$num_cores" command_function ::: "${values[@]}"
else
    # Run the commands sequentially
    for value in "${values[@]}"; do
        command_function "$value"
    done
fi