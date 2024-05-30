#!/bin/bash -e

# Function to show usage
usage() {
    echo "Usage: $0 -v values -j num_cores"
    echo "  -v values      A comma-separated list of values (e.g., 40,60,80,90)"
    echo "  -j num_cores   Number of cores to use in parallel"
    exit 1
}

# Ensure both options are provided
if [ $# -eq 0 ]; then
    usage
fi

# Parse command line arguments
while getopts "v:j:" opt; do
    case $opt in
        v)
            IFS=',' read -r -a values <<< "$OPTARG"
            ;;
        j)
            num_cores="$OPTARG"
            ;;
        *)
            usage
            ;;
    esac
done

# Ensure required arguments are set
if [ -z "${values+x}" ] || [ -z "${num_cores+x}" ]; then
    usage
fi
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

# Export the function to be used in parallel
export -f command_function

# Run the commands in parallel with 3 cores
parallel -j ${num_cores} command_function ::: "${values[@]}"
