
########################################################
# Choose models based on Gu et al. data for simulation #
########################################################

# For a simulation that checks on the benefit of using the 
# time-dose two-parameter log-logistic (td2pLL) model,
# we first have to decide on the "true" models that will be used in
# the data generating process

# We generally want to use three models:

# (M1): Clear dependency of EC50 on exposure time
# (M2): Little dependency of EC50 on exposute time
# (M3): No dependency of EC50 on exposure time

# The compounds for M1 and M2 were chosen subjectively by
# fitting the td2pLL to all compounds and, if it converged,
# choosing a model with a rather large \gamma (M1) and one with a
# rather small \gamma (M2).
# (\gamma represents degree of dependency of the EC50 on epxosure time).

# For (M1) and (M2) we will almost entirely (apart from scaling) use models 
# fitted to actual data of a compound of the Gu et al.data.
# For (M3), we will also base the model on the Gu et al. data but with the
# modification of setting \gamma to zero.
# Hence, using the ext2pLL model for fitting data generated from the M3 model
# is actually unnecessary, as setting \gamma to zero erradicates the depdendency
# of the EC50 on exposure time.

#########################################################

# Load the Gu et al. data: 
# Use the pre-processed data where the refitting was applied
# (cf. ./02_data_prep/prepare_data.R)
load("./02_data_prep/data_refit_prep.RData")

# Load functions
source("./00_functions/functions.R")


# Seed for identical results during fitting process
set.seed(1905)


# For (M1) use data from compound CHL
# Scale dose so that dose ranges in[0,1]

d_M1 <- data_refit_prep %>% filter(compound == "CHL") %>%
  dplyr::select(time, dose, resp) %>%
  mutate(dose = dose / max(dose))

M1 <- fit_td2pLL(data = d_M1)

# Take a look at the three dimensional model if you like:
plot(M1, xaxis_scale = "log", xaxis_title = "Concentration",
     yaxis_title = "Exposure duration", zaxis_title = "Response")

# For (M2) use data from compound ETOH

d_M2 <- data_refit_prep %>% filter(compound == "ETOH") %>%
  dplyr::select(time, dose, resp) %>%
  mutate(dose = dose / max(dose))

M2 <- fit_td2pLL(data = d_M2)

# Take a look at the three dimensional model if you like:
plot(M2, xaxis_scale = "log", xaxis_title = "Concentration",
     yaxis_title = "Exposure duration", zaxis_title = "Response")


# For (M3) use model (M2) but with gamma = 0

M3 <- c(coef(M2)["h"], coef(M2)["delta"], gamma = 0, coef(M2)["c0"])


# Take a look at the three dimensional model if you like:
plot_td2pLL(M3, dose_lim = c(1e-3, 1), time_lim = c(1, 7),
              xaxis_title = "Concentration", yaxis_title = "Exposure duration")


# For the graphics in the manuscript to be static, 
# the dynamic plots of the models M1 and M2 in the viewer were turned
# and then a snapshot was taken with the "export" option.


########
# Save #
########

true_models <- list(M1 = M1, M2 = M2, M3 = M3)
true_models_data <- list(d_M1 = d_M1, d_M2 = d_M2)

save(M1, M2, M3, file = "./03_simulation_prep/01_true_models.RData")
save(d_M1, d_M2, file = "./03_simulation_prep/01_true_models_data.RData")

