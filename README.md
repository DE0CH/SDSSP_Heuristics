# SDSSP_Heuristics
Heuristics for the Star Discrepancy Subset Selection Problem

The four shift_ files contain the four instantiations of the heuristic.
shift_TAnobrute is the TA variant with no brute-force check.
shift_TA is the Ta variant with a brute force check.
shift_v2 is the DEM variant with a  brute-force check.
Finally, shift_v2no brute is the DEM variant without a brute force-check.

Inputs are the input pointfile, the dimension, the initial number of points n, the final number of points k and the output file.

point.c generates the 512 point Sobol set using the GNU Scientific Library, it requires the dimension as input.
