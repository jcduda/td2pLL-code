In this directory, the data pre-processing takes place.

The file prepare_data.R 

	-requires data from: 
		- ./01_data_raw 

	-generates: 
		-./data_prep.rds 
		-./data_refit_prep.rds

In the follwoing analyses steps, only ./02_data_prep/data_refit_prep.rds is used.

./data_refit_prep.rds includes all pre-processing of ./data_prep.rds.

The difference is, that in ./data_refit_prep.rds, an additional pre-processing
step, the refitting, is applied. Details can be found in prepare_data.R.

