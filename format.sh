#!/bin/bash -e

usage() {
    echo "Usage: $0 csv_file dimension values [-j num_cores]"
    echo "  csv_file     The CSV file to process (e.g., df_crit.csv)"
    echo "  dimension    Dimension value (e.g., 3)"
    echo "  values       A comma-separated list of values (number of points in the subset) (e.g., 40,60,80,90)"
    echo "  -j num_cores (Optional) Number of cores to use in parallel (default is 1)"
    exit 1
}

# Ensure at least 3 positional arguments are provided
if [ "$#" -lt 3 ]; then
    usage
fi

# Assign positional arguments
csv_file="$1"
shift
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
    python3 ../match_index.py "$csv_file" "subset_${value}.txt" "subset_${value}.csv" /dev/null 
    head -n 1 "subset_${value}.txt" | grep -oP 'discrepancy=\K[\d.]+' | awk '{print "k=" ENVIRON["value"], $0}'
}

# Export the function to be used in parallel
export -f command_function

# Run the commands either in parallel or sequentially based on num_cores
if [ "$num_cores" -gt 1 ]; then
    parallel -j "$num_cores" command_function ::: "${values[@]}"
else
    for value in "${values[@]}"; do
        command_function "$value"
    done
fi