####################################
######## 04_run_simulation #########
####################################

This directory requires:
	- ./00_functions/functions.R
	- ./03_simulation_prep/01_true_models.RData
	- ./03_simulation_prep/02_true_designs.RData
	- ./03_simulation_prep/03_lm_sd_noise.RData

Run the following script line by line. You might have batchtool installed and configured so it might run vor very long.
Consider just loading the simulation results to generate the figures.

	- run_simulation
		- requires:
			- ./00_functions/functions.R
			- ./03_simulation_prep/01_true_models.RData
			- ./03_simulation_prep/02_true_designs.RData
			- ./03_simulation_prep/03_lm_sd_noise.RData
		- generates:
			- results_simulation.RData
				- data.table / data.frame containing all the simulation results generated with batchtools().
				  The following columns are relevant for the publication:
					- algorithm 
						- anova = 2 step pipeline
						- no_pre = either single 2pLL, seperate 2pLL or td2pLL. Specified in "always"
					- always
						- joint = single 2pLL
						- sep = separate 2pLL
						- td2pLL = td2pLL (without the 2-step-pipeline)
					- model_id
						- M1 = M2 in pulication (strong time-dpendency)
						- M2 = M2 in pubcliation (small/medium time-dependency)
						- M3 = M0 in publication (no time-dependency)
					- noise_id
						- N1 =  little noise
						- N2 = medium noise
						- N3 = much noise
					- n_times_id 
						- 3 (3 different exposure durations)
						- 4 (4 different exposure durations)
	