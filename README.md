# td2pLL-code
Code to reproduce the simulation results of td2pLL: An intuitive time-dose-response model.

Maintainer: Julia Duda (contact: duda@statistik.tu-dortmund.de)

To reproduce the results, follow these steps

## Step 0

Download the repository and set the working directory to where the files are.

## Step 1: Functions and packages

Open and run 00_functions/functions.R. to load required packages and source
functions needed to the simulation.

Note that you might need to install some package first before you can use them from
your library.

## Step 2: Pre-process raw data (Optional)

The raw data are stored in 01_raw_data. This is the format available for download in
the original publication from Gu et al. 2018 (doi: https://doi.org/10.1007/s00204-018-2302-0).

To pre-process the data (normalize with negative control response mean and
rescale with initial 4pLL fit upper asymptotes), 

- open an run 02_data_prep/prepare_data.R to
- generate data_refit_prep.rds

This step is optional because the pre-processed data are already available.
However, feel free to re-generate them or look at the pre-processing code.

## Step 3: Simulation preparation (Optional)

This step is optional, because the simulation preparation files are already available.
However, feel free to re-generate them or look at the code.

To do so, open and run the following R-files:

- 03_simulation_prep/01_choose_true_models.R
- 03_simulation_prep/02_choose_experimental_designs.R
- 03_simulation_prep/03_choose_sd_noise.R

For details and the genrated files, open the 03_simulation_prep/README.txt.

Note: The optimal design considerations were not of interest for the publication
of the td2pLL model of Duda et al. (2022) and can be ignored.

## Step 4: Run the simulation (Optional)

Again, this step is optional as the simulation results are available.
However, feel free to re-generate them or look at the code.

Open and run the following R-file

- 04_simulation_run/run_simulation.R

The simulation results are saved in

- 04_simulation_run/results.simulation.RData

Note that the simulation was run on a high performance computing (HPC) cluster LiDo 3 of the
TU Dortmund University with the help of the R package batchtools for convenient batchjob
organization. 

This means that if you run the code on your local computer, it might take very long.
And if you want to re-run on another HPC cluster, make sure to setup the configuration
files properly.

Note that in the code, the model M1 means strong time-dependency,
M2 is medium/small time-dependency and M3 means no time-dependency.
The naming was changed for the manuscript, where M2 means strong time dependency, M1 medium/small time-dependency and M0 means no time-dependency.





