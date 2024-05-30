# SDSSP_Heuristics

This code is adapted from <https://github.com/frclement/SDSSP_Heuristics> for specific needs of (insert paper name). 

## Dependencies:
1. Python 3
1. gcc and associated libaries, as well as common command line utilities (should be included in common linux distros)
1. [`parallel`](https://manpages.ubuntu.com/manpages/jammy/man1/parallel.1.html) tool on linux (only if you want to run in parallel)


## Steps to reproduce:
1. Compile shift_v2nobrute.c
```bash
gcc shift_v2nobrute.c -lm -O3
```
1. Make a directory. The script will litter the directory with files and overwrite things without warning. 
```bash
mkdir run-data
```
1. cd to the directory. IMPORTANT: the scripts will overwrite files in cwd without warning.
```bash
cd run-data
```
1. put your df_crit.csv in the directory and run (for example). See the help of each files for usage
```bash
python3 ../extract_csv.py df_crit.csv df_crit.txt
../run.sh -j 10 ../a.out df_crit.txt 3 2188 30,60,80,90,100,110,120
../format.sh df_crit.csv 30,60,80,90,100,110,120 && ../split.sh -j 10 ../a.out 3 30,60,80,90,100,110,120
```

