#####################################
######## 03_simulation_prep #########
#####################################

This directory requires:
	- ./02_data_prep/data_refit_prep.RData
	- ./00_functions/functions.R

Based on the pre-processed data of Gu et al., the true models for the simulation (M1, M2 and M3),
as well as the considered designs and the added noise are calculated.

The scripts have to be run in the following order:
	- 01_choose_true_models
		- requires:
			- ./02_data_prep/data_refit_prep.RData
			- ./00_functions/functions.R
		- generates:
			- 01_true_models.RData
				- List with models M1, M2 and M3. 
				  M1 and M2 are objects of class "DERmod" and "nls",
				  M3 only contains coefficients
			- 01_true_models_data.RData
				- List with d_M1 and d_M2, real data from Gu et al. of the
					corresponding compound, model M1 and M2 were fitted to.

	- 02_choose_experimental_design
		- requires:
			- ./02_data_prep/data_refit_prep.RData
			- ./00_functions/functions.R
		
		- generates:
			- 02_true_designs.RData
				- data.frame including all designs in a long format. Columns are:
					- time: (exposure time) 1, 2, 4 or 7
					- dose
					- n	: number of replicates at this (time, dose) point
					- design: D_equ_3_n***, D_equ_4_n***, D_opt_3_n_***, D_opt_4_n***
					- n_obs : Total nuber of observations in this design

	- 03_choose_sd_noise
		- requires:
			- ./02_data_prep/data_refit_prep.RData
			- ./00_functions/functions.R

		- generates:
			- 03_plot_sd_values.RData
				- ggplot object: Scatter plots of sd values of Gu et of measurements
					per compound, Donor, exposure time and dose. Plots are stratified
					by exposure time. Linear models (y~expo+expo^2) are added, they
					are calculated seperately for each exposure time and are not the
					linear model which will be used to generate the noise sd values later.
			- 03_lm_sd_noise.RData
				- lm object: Actual noise sd model used later to generate the added noise.
			
