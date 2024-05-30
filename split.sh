#!/bin/bash -e

# Function to show usage
usage() {
    echo "Usage: $0 dimension values [-j num_cores]"
    echo "  dimension    Dimension value (e.g., 3)"
    echo "  values       A comma-separated list of values (number of points in the subset) (e.g., 40,60,80,90)"
    echo "  -j num_cores (Optional) Number of cores to use in parallel (default is 1)"
    exit 1
}

# Ensure at least 2 positional arguments are provided
if [ "$#" -lt 2 ]; then
    usage
fi

# Assign positional arguments to dimension and values
dimension="$1"
shift
IFS=',' read -r -a values <<< "$1"
shift

# Default number of cores
num_cores=1

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
    # make a temp file
    tpfile=$(mktemp)
    python3 ../extract_csv.py subset_${value}.csv ${tpfile}
    ./a.out ${tpfile} 3 ${value} $(expr ${value} / 2) half_subset_${value}.txt 2>&1 | ts > log_half_subset_${value}.txt
    python3 ../match_index.py subset_${value}.csv half_subset_${value}.txt half_subset_${value}.csv other_half_subset_${value}.csv
    python3 ../extract_csv.py other_half_subset_${value}.csv other_half_subset_${value}.txt
    ./ao.out other_half_subset_${value}.txt 3 ${value} $(expr ${value} / 2) other_half_subset_${value}.txt 2>/dev/null > /dev/null
    head -n 1 half_subset_${value}.txt | python3 -c 'import sys; import re; print("half subset:", f"k={sys.argv[1]}", re.search(r"discrepancy=([\d\.]+)", input())[0])' ${value}
    head -n 1 other_half_subset_${value}.txt | python3 -c 'import sys; import re; print("other half subset:", f"k={sys.argv[1]}", re.search(r"discrepancy=([\d\.]+)", input())[0])' ${value}
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