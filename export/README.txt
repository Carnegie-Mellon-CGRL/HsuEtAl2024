build_smf.py: a python script that builds an SMF input file for a given number of dimensions 

run_script_gnr.m: runs a single trial of the optimization.  objective and constraint functions are selected in this script. 
master: a bash script to run 10 trials of the optimization by running run_script_gnr in parallel - given a constraint and objective

rad_objective.m: an objective function only based on radius
dual_objective.m: an objectve function based on radius and compliance
triple_objective.m: an objective function based on radius, compliance, and inflammation

full_plot.m: generate matlab plots of radius, compliance, thickness, and inflammation

minrad_constraint.m: a constraint function keeping radius above 100% of native radius
radius_bound_constraint.m: a constraint function keeping radius between 90% and 110% of native radius
radius_min_cap_constraint.m: a constraint function keeping radius between 100% and 110% of native radius



