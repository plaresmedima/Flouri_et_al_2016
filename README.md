# Linear-2CM

Description: Simulation code for anonymised manuscript submitted to Magn Reson Med 2015, with title "Fitting the two-compartment model in DCE-MRI by linear inversion". This repository contains the IDL code to reproduce the results, figures and tables in the manuscript. 

To reproduce the figures and tables:

1) Set your IDL search path (preferences) to the folder "Linear-2CM", and include all its subfolders.
2) Change the path in the procedure SIM_PATH.pro to the folder where you want to store the results. 
3) Compile and run the procedure "ALL_FIGURES_AND_TABLES.pro". 

This will create all figures and tables at 10.000 simulations, which may take several hours. To generate a less precise set of results in a shorter time, set the parameter "nSim" in "ALL_FIGURES_AND_TABLES.pro" to a smaller number, eg. 100 or 1000.

Subfolders:

"figures_and_tables": The code for the individual figures and tables. They can all be compiled and run separately.
"measurement": The code to generate the data (AIF, exact tissue concentrations and measured concentrations)
"model_fitting": Core library with model fitting routines (LLS and NNLS) for both models 2CXM and 2CFM. This also contains the package "mpfit" by C. Markwardt, included with permission. 
"output": Default folder for storing output (initialised with final figures and tables).
