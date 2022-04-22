###################################################
# Choose experimental designs used for simulation #
###################################################

# In this script we will decide on the different experimental designs that
# will be used to generate the data in the simulation.
# Each experimental design is defined by a data.frame consisting of
# columns time, dose, n, n_obs

# The name of a design is constructed as of D_[equ/opt]_[3/4]_n[216/72]
#   - the letter D (for Design),
#   - _equ_ indicates an equidistant design on log(sqrt(10)) scale with balanced replicates
#     _opt_ indicates a pseudo bayesian, local c-optimal design for each time point
#           separately
#   - _3_ indicated that 3 exposure times are used (1, 2, and 7) as in
#         Gu et al. for the compound used for model M1 (CHL)
#     _4_ indicates that 4 exposure times are used (1, 2, 4 and 7)
#   - _n72_ / _n216_ indicated the total number of obserations

# The total number of observations, n_obs is either 
#   - 216, as  in Gu et al., or 
#   - 72,
# which is a third of the original amount of experiments and might
# be more realistic for toxicological data. 

#######################
# load required files #
#######################


source("./00_functions/functions.R")

# original data of Gu et al. used for model M1
# (Used here to extract structure of experimental design)
load("./03_simulation_prep/01_true_models_data.RData")
load("./02_data_prep/data_refit_prep.RData")

#######################################
# D_equ_3 & D_equ_4: equistant designs#
#######################################

# D_equ_3, n_obs = 216
# original design as used for M1 (compound CHL):
# n_obs = 3*6*12 = 216
D_equ_3_n216 <- d_M1 %>%
  arrange(time, dose) %>%
  group_by(time, dose) %>%
  summarize(n = n(), .groups = "drop") %>%
  dplyr::select(time, dose, n) %>%
  dplyr::mutate(design = "D_equ_3_n216", n_obs = 216)



# D_equ_3_n72
# design like D1a, but only 1/3 of original data, i.e. n_obs = 216 / 3 = 72 = 3*6*4
D_equ_3_n72 <- D_equ_3_n216 %>%
  dplyr::mutate(n = 4,
                design = "D_equ_3_n72") %>%
  dplyr::mutate(n_obs = sum(n))

# D_equ_4_n216
# close to original design but with additional expo = 4 and 
# instead reduced number of technical replicates
# n_obs = 4*6*9 = 216 

doses_low_time = D_equ_3_n216$dose[D_equ_3_n216$time == 1]
doses_high_time = D_equ_3_n216$dose[D_equ_3_n216$time == 7]

D_equ_4_n216 <- rbind.data.frame(D_equ_3_n216 %>% mutate(n = 9) %>%
                                   dplyr::select(-c("design","n_obs")),
      data.frame(time = 4, dose = doses_high_time, n = 9)) %>%
  arrange(time, dose) %>%
  mutate(design = "D_equ_4_n216", n_obs = 216)

# D_equ_4_n72
# design like D_equ_4_n216, but only 1/3 of original data,
# i.e. n_obs = 216/3=72=4*6*3

D_equ_4_n72 <- D_equ_4_n216 %>%
  dplyr::mutate(n = 3,
                design = "D_equ_4_n72") %>%
  dplyr::mutate(n_obs = sum(n))


#################################################
# D_opt_3_n*** & D_opt_4_n***: D-Optimal designs#
#################################################

# Using a locally D-optimal design with (equally) assumed shapes:
# h = 2, ED50 either 0.001, 0.01, 0.1 for all exposures.
# Work with ICAOD package andno more with DoseFinding package.

# First figure out how many support points necessary for optimal design.

# k=9 just still too large

########### 
# k = 9 ###
#
# set.seed(1905) 
# icaod_des9 <- ICAOD::robust(lx = 1e-10, ux = 1,
#                             prob = c(1/3, 1/3, 1/3),
#                             parset =  cbind(
#                               100,                 # b1
#                               0,                   # b2
#                               c(0.001, 0.01, 0.1 ),# b3
#                               2),                  # b4
#                             k = 9,
#                             iter = 2000,
#                             ICA.control = list(ncount = 150,
#                                                stop_rule = "equivalence"),
#                             fimfunc = FIM_sig_emax)
#  
#  
# my_sigEmax_equ_plot(icaod_des9) 
# k=8 and not k=9 design points clearly optimal!
##############

# k=8 correct

# k = 8 
# set.seed(1905) 
# icaod_des <- ICAOD::robust(lx = 1e-6, ux = 1,
#                            prob = c(1/3, 1/3, 1/3),
#                            parset =  cbind(
#                              100,                 # b1
#                              0,                   # b2
#                              c(0.001, 0.01, 0.1 ),# b3
#                              2),                  # b4
#                            k = 8,
#                            iter = 1500,
#                            ICA.control = list(ncount = 150,
#                                               stop_rule = "equivalence"),
#                            fimfunc = FIM_sig_emax)
#  
# save(icaod_des, file = "./03_simulation_prep/02_icaod_des.RData")
load("./03_simulation_prep/02_icaod_des.RData")

# If the line is below zero and touches zero exactly at the red points,
# the design is D-optimal - yay!

my_sigEmax_equ_plot(icaod_des)

# Interesting:
# Seems like the optimal support points are: 
# 0, ED50_guess / (1 + 2/3), ED50_guess * (1 + 2/3), ..., 1

# Comparison: Original plotting function does not allow log-scaling for x-axis:
# ICAOD:::plot.minimax(icaod_des)

# Note: k=7 not optimal (checked it)

opt_doses <- my_rndDesign(icaod_des, n = 72)$dose 

###################################
# D_opt_3_n216 and D_opt_3_n72    #
###################################


# exp. Design for original 3 exposure levels 1, 2 and 7

# discretized values using n = 72 because 72*3 = 216 observations are used in total or
# discretized values using n = 24 because 24*3 = 72 observations are used in total

# No more needed, remove when rest updated
# n_of_D_opt_3_n216 <- rndDesign(design_opt, n = 72, eps = 0.02)[ind]
# n_of_D_opt_3_n72 <- rndDesign(design_opt, n = 24, eps = 0.02)[ind]


# Create the experimental plan:
D_opt_3_n216 <- data.frame(time = rep(c(1, 2, 7), each = length(opt_doses)),
                  dose = opt_doses,
                  n = my_rndDesign(icaod_des, n = 72)$n,
                  design = "D_opt_3_n216", n_obs = 216)

# STOPPED HERE
D_opt_3_n72 <- data.frame(time = rep(c(1, 2, 7), each = length(opt_doses)),
                        dose = opt_doses,
                        n = my_rndDesign(icaod_des, n = 24)$n) %>%
              mutate(design = "D_opt_3_n72",
                     n_obs = sum(n))



# D_opt_4_n216 and D_opt_4_n72
# exp. Design for 4 exposure levels 1, 2, 4 and 7:
# discretize values using n = 54 because 4*54 = 216 observations are used in total or
# discretize values using n = 18 because 4*18 = 72 observations are used in total

# No more needed, remove when done
# n_of_D_opt_4_n216 <- rndDesign(design_opt, n = 54, eps = 0.02)[ind]
# n_of_D_opt_4_n72 <- rndDesign(design_opt, n = 18, eps = 0.02)[ind]

D_opt_4_n216 <- data.frame(time = rep(c(1, 2, 4, 7), each = length(opt_doses)),
                  dose = opt_doses,
                  n = my_rndDesign(icaod_des, n = 54)$n,
                  design = "D_opt_4_n216", n_obs = 216)

D_opt_4_n72 <- data.frame(time = rep(c(1, 2, 4, 7), each = length(opt_doses)),
                         dose = opt_doses,
                         n = my_rndDesign(icaod_des, n = 18)$n,
                         design = "D_opt_4_n72") %>%
  mutate(n_obs = sum(n))

##########
# Gather #
##########

# all designs in one long-format

# designs <- rbind.data.frame(D1a, D1b, D2a, D2b,
#                            D1a_small, D1b_small, D2a_small, D2b_small)

designs <- rbind.data.frame(D_equ_3_n72,
                            D_equ_3_n216,
                            D_equ_4_n72,
                            D_equ_4_n216,
                            D_opt_3_n72,
                            D_opt_3_n216,
                            D_opt_4_n72,
                            D_opt_4_n216) %>%
  dplyr::mutate(spacing = gsub("D_|_\\d.*","", design),
                n_times = as.numeric(gsub(".*equ_|.*opt_|_n.*", "", design)))

save(designs, file = "./03_simulation_prep/02_true_designs.RData")


# Visualize the designs to take a look:
 plot_designs(designs = designs %>% filter(n_obs == 216))



                     