################################################
# Running the simulation for the  td2pLL model #
################################################

# In this simulation, the time-dose 2pLL model
# is applied and compared to other possible dose-response
# modeling approaches that account for exposure time.

# Parameter of interest is the ED50 (dose where 50% of effect is obtained).

# Performance goal with which we compare the methods are based
# on how accurate the ED50 is calculated.

# Primary Performance Measures:
# AMAFC (Average Mean Absolute Folg Change)
# AMAcc(Average Mean Accepted)
# Note: AMAcc not used in publication.

# Secondary Performance Measures (Note: Not used for publication):
# Convergence-problems rates

# The simulation is run on LiDo3 of the TU Dortmund University
# using batchtools.

# The simulation consists of a data generation part and a method part.
# For each simulation run, data is generation with one of the
# data generation settings. Then, the generated data is analyzed with
# one of the methods.

# For each simulation scenario, 1000 replicates are calculated.



library(batchtools)
#############
# if on LiDo (set working directory accordingly)
#############

if(getwd() == "/home/smjududa") setwd("/work/smjududa/Projekte/exp_dose_resp/")



# Considerations for using batchtools: Prepare Registry

reg = makeExperimentRegistry(file.dir = "./04_simulation_run/registry",
                             packages = c("td2pLL", "dplyr", "tidyr", "scales", "stats", "plotly", "drc", "data.table"),
                             source = "./00_functions/functions.R",
                             load = c("./03_simulation_prep/01_true_models.RData",
                                    "./03_simulation_prep/02_true_designs.RData",
                                    "./03_simulation_prep/03_lm_sd_noise.RData"),
                             seed = 1)

####################
# Data Generation  #
####################

# for data generation, the factors 
# model (M1, M2, M3),
# noise (N1, N2, N3) and
# experimental design (D_equ_3_n***, D_equ_4_n***, D_opt_3_n_***, D_opt_4_n***) 
#   which considers the aspects:
#   - spacing ("equ", "opt")
#   - n_times (3, 4)
#   - n_obs (72, 216) 
# is considered. Hence, there are
# 3*3*2*2*2=72 Data generation scenarios

generate_data <- function(data, job, model_id, noise_id,
                          # expDes_id,
                          spacing_id,
                          n_times_id,
                          n_obs_id, ...){
  
  stopifnot(model_id %in% c("M1", "M2", "M3"))
  if(model_id == "M1") model <- coef(M1)
  if(model_id == "M2") model <- coef(M2)
  if(model_id == "M3") model <- M3
  
  stopifnot(noise_id %in% c("N1", "N2", "N3"))
  
  # stopifnot(expDes_id %in% c("D1a", "D1b", "D2a", "D2b"))
  stopifnot(spacing_id %in% c("equ", "opt"))
  stopifnot(n_times_id %in% c(3, 4))
  stopifnot(n_obs_id %in% c(216, 72))
  
  expDes <- designs %>% dplyr::filter(
    #design == expDes_id,
    spacing == spacing_id,
    n_times == n_times_id,
    n_obs == n_obs_id)
  
  list(stoch_data = generate_data(model = model, noise_id = noise_id, expDes = expDes),
       model_id = model_id) 
}

addProblem(name = "generate_data", data = NA, fun = generate_data, seed = 42)


# The algorithms: anova and no_pre #####

# Given a data generating scenario, the method (or algorithm, pipeline)
# is used for estimating a dose-response relationship and 
# specifically, the ED50 parameter.

# Considered method for a given, generated data set are:

# 1. "anova": Two-step procedure:
#           First, an anova pretest check if exposure time has influence
#           Second, based on pre-test, td2pLL or single "joint" dose-response
#               is applied
# "no_pre" with:
#
# 2.  "all_joint": No pre-test, always a single dose-response curve fitted
# 3.  "all_sep":   No pre-test, always separate dose-repsonse curves 
#                   (one per exposure time) fitted
# 4.  "all_td2pLL" No pre-test, always the td2pLL fitted


anova.wrapper = function(data, job, instance, alpha, ...){
   #in instance$stoch_data, the generated data should be saved
  anova_signif <- tryCatch(
    {
      td2pLL::td2pLL_anova(data = instance$stoch_data, alpha = alpha)
    },
    error = function(cond){
      return(list(MAFC = "error", 
                  MAcc = "error",
                  signif = "error",
                  is_conv = FALSE,
                  ED50s_at_times = "error"))
    }
    )
  
  if(unlist(anova_signif["signif"]) == 1){
    # if differences: fit td2pLLL
    model_fit <- tryCatch(
      {
        td2pLL::fit_td2pLL(data = instance$stoch_data)
      },
      error = function(cond){
        return(list(MAFC = "error", 
                    MAcc = "error",
                    signif = anova_signif["signif"],
                    is_conv = FALSE,
                    ED50s_at_times = "error"))
      }
      )
    
    # return with MAFC = "error", if no final fit
    if(!is.null(model_fit$MAFC)) return(model_fit)
    
    is_conv = model_fit$convInfo$isConv
    # get ED50 at all time values
    ED50s_at_times <- get_ED50s(coefs = coef(model_fit),
                                times = unique(instance$stoch_data$time)
                                )
    
  } 
  
  if(unlist(anova_signif["signif"]) == 0 | unlist(anova_signif["signif"]) == "error"){
    # if no differences or error: fit joint version (LL2.2 is used!)
    model_fit <- tryCatch(
      {
        td2pLL:::fit_joint_2pLL(data = instance$stoch_data)
      },
      error = function(cond){
        return(list(MAFC = "error",
                    MAcc = "error",
                    signif = anova_signif["signif"],
                    is_conv = FALSE,
                    ED50s_at_times = "error"))
      }
    )
    
    # return  MAFC = "error", if no final fit
    if(!is.null(model_fit$MAFC)) return(model_fit)
    
    is_conv <- model_fit$fit$convergence

    if(model_fit$fct$name[3] != "LL2.2") 
      stop("Following code assumes LL2.2, which is not at hand.")
    
    # get ED50 at all time values
    ED50s_at_times <- data.frame(time = unique(instance$stoch_data$time),
                                 ED50 = exp(coef(model_fit)[2]))
  }

  
  # compare to true ED50 values
  if(instance$model_id %in% c("M1", "M2")) true_coefs <- coef(get(instance$model_id))
  if(instance$model_id == "M3") true_coefs <- get(instance$model_id)
  
  ED50s_at_times$true_ED <- td2pLL::get_ED50s(coefs = true_coefs, times = ED50s_at_times$time)$ED50
  # Mean Absolute Fold Change
  MAFC <- mean(abs(log2(ED50s_at_times$ED50 / ED50s_at_times$true_ED))) 
  # Mean Acceptable (absolute ratio within [0.9, 1.1] )
  # MAcc <- mean(abs(ED50s_at_times$ED50 / ED50s_at_times$true_ED - 1) < 0.1)
  # Changed to using log2 ratios
  MAcc <- mean(abs(log2(ED50s_at_times$ED50 / ED50s_at_times$true_ED)) < log2(1.1))
  
  return(list(MAFC = MAFC,
              MAcc = MAcc,
              signif = anova_signif["signif"],
              is_conv = is_conv,
              ED50s_at_times = ED50s_at_times))
}

# add the anova algorithm to the experiment 
addAlgorithm(name = "anova", fun = anova.wrapper)

# define function for remaining algorithms (methods)
no_pre.wrapper = function(data, job, instance, always, ...){
  
  stopifnot(always %in% c("td2pLL", "sep", "joint"))
  
  if(always == "td2pLL"){
    model_fit <- tryCatch(
      {
        td2pLL::fit_td2pLL(instance$stoch_data)
      },
      error = function(cond){
        return(list(MAFC = "error", 
                    MAcc = "error",
                    signif = "no_pre",
                    is_conv = FALSE,
                    ED50s_at_times = "error"))
      }
    )
    # if MAFC=="error" return without MAFC calculation
    if(!is.null(model_fit$MAFC)) return(model_fit)
    
    is_conv = model_fit$convInfo$isConv
    # get ED50 at all time_values
    ED50s_at_times <- get_ED50s(coefs = coef(model_fit),
                                times = unique(instance$stoch_data$time)
    )
  }
  
  # "sep" uses LL2.2, too! 
  if(always == "sep"){
    model_fit <- tryCatch(
      {
        # Do NOT use td2pLL:::fit_sep_2pLL(), because it shares the h !
        # But here, we want completely separated curve fits.
        fit_sep_2pLL(data = instance$stoch_data)
      },
      error = function(cond){
        return(
          list(MAFC = "error", 
               MAcc = "error",
               signif = "no_pre",
               is_conv = FALSE,
               ED50s_at_times = "error")
        )
      }
    )
    if(!is.null(model_fit$MAFC)) return(model_fit)

    is_conv <- model_fit$fit$convergence
    
    
    # get ED50 related parameters of the seperate model:
    if(model_fit$fct$name != "LL2.2")
      stop("Following code assumes LL2.2")
    ED50s_pre <- coef(model_fit)[grep("e:", names(coef(model_fit)), value = T) ]
    ED50s <- c(exp(ED50s_pre[1]), exp(ED50s_pre[1] + ED50s_pre[-1]))
    ED50s_at_times <- data.frame(time = unique(instance$stoch_data$time),
                                 ED50 = ED50s, row.names = NULL)
  }
  
  # joint uses LL2.2
  if(always == "joint"){
    model_fit <- tryCatch(
      {
        td2pLL:::fit_joint_2pLL(data = instance$stoch_data)
      },
      error = function(cond){
        return(
          list(MAFC = "error",
               MAcc = "error",
               signif = "no_pre",
               is_conv = FALSE,
               ED50s_at_times = "error")
        )
      }
    )
    if(!is.null(model_fit$MAFC)) return(model_fit)
    
    is_conv <- model_fit$fit$convergence
    if(model_fit$fct$name[3] != "LL2.2")
      stop("Following code assumes LL2.2")
    ED50s_at_times <- data.frame(time = unique(instance$stoch_data$time),
                                 ED50 = exp(coef(model_fit)[2]), row.names = NULL)
    
  }
  
  # compare to true ED50 values
  if(instance$model_id %in% c("M1", "M2")) true_coefs <- coef(get(instance$model_id))
  if(instance$model_id == "M3") true_coefs <- get(instance$model_id)
  
  ED50s_at_times$true_ED <- get_ED50s(coefs = true_coefs, times = ED50s_at_times$time)$ED50

  # Mean Absolute Fold Change
  MAFC <- mean(abs(log2(ED50s_at_times$ED50 / ED50s_at_times$true_ED))) 
  # Mean Acceptable (absolute ratio within [0.9, 1.1] )
  # MAcc <- mean(abs(ED50s_at_times$ED50 / ED50s_at_times$true_ED - 1) < 0.1)
  # Changed to using log2 ratios
  MAcc <- mean(abs(log2(ED50s_at_times$ED50 / ED50s_at_times$true_ED)) < log2(1.1))
  
  return(list(MAFC = MAFC,
              MAcc = MAcc,
              signif = "no_pre",
              is_conv = is_conv,
              ED50s_at_times = ED50s_at_times))
}

# add the methods withour a pre-test to the experiment
addAlgorithm(name = "no_pre", fun = no_pre.wrapper)

##################
# Problem design #
##################

# data generation (problem design for batchtools) is factorial
pdes = list(generate_data = CJ(model_id = c("M1", "M2", "M3"),
                               noise_id = c("N1", "N2", "N3"),
                               # expDes_id = c("D1a", "D1b", "D2a", "D2b"),
                               spacing_id = c("equ", "opt"),
                               n_times_id = c(3, 4),
                               n_obs_id = c(216, 72)
)
)

####################
# algorithm design #
#################### 

ades = list(
  anova = data.table(alpha = 0.05),
  no_pre = data.table(always = c("td2pLL", "sep", "joint"))
)


# change repls to 1 for a minimal simulation with 72*4=288 small jobs

addExperiments(pdes, ades, repls = 1000)

# take a look:
summarizeExperiments(by = c("problem", "algorithm",
                            "model_id",
                            "always", "noise_id", "spacing_id", "n_times_id", "n_obs_id"))

######################################
# Optional: testing before submission
######################################

# model_id = "M1"; noise_id= "N1"; spacing_id = "equ"; n_times_id = 3; n_obs_id = 216
# easy scenario: Clear time-dependency, little noise, large n
id = head(findExperiments(algo.name = "anova",
                          prob.pars = (model_id == "M1" &
                                         noise_id == "N1" &
                                         spacing_id == "equ" &
                                         n_times_id == 3 &  
                                         n_obs_id == 216)),
          1)

testJob(id= id)

# more difficult: Much noise, small n, no time-dependency
# model_id = "M3"; noise_id= "N3"; spacing_id = "equ"; n_times_id = 4; n_obs_id = 72
id = head(findExperiments(algo.name = "anova",
                          prob.pars = (model_id == "M3" &
                                         noise_id == "N3" &
                                         spacing_id == "equ" &
                                         n_times_id == 4 &
                                         n_obs_id == 72)),
          1)

testJob(id= id)

# Separate curves, large noise
# model_id = "M3"; noise_id= "N3"; spacing_id = "equ"; n_times_id = 3; n_obs_id = 72
id = head(findExperiments(algo.name = "no_pre",
                          algo.pars = (always == "sep"),
                          prob.pars = (model_id == "M3" &
                                         noise_id == "N3" &
                                         spacing_id == "opt" &
                                         n_times_id == 3 &
                                         n_obs_id == 72)),
          1)
testJob(id = id)

# separate curves, little noise
# model_id = "M1"; noise_id= "N1"; spacing_id = "opt"; n_times_id = 3; n_obs_id = 216
id = head(findExperiments(algo.name = "no_pre",
                          algo.pars = (always == "sep"),
                          prob.pars = (model_id == "M1" &
                                         noise_id == "N1" &
                                         spacing_id == "opt" &
                                         n_times_id == 3 &
                                         n_obs_id == 216)),
          1)
testJob(id = id)

# Join curves, strong time-dep., much noise, small n
id = head(findExperiments(algo.name = "no_pre",
                          algo.pars = (always == "joint"),
                          prob.pars = (model_id == "M1" &
                                         noise_id == "N3" &
                                         spacing_id == "equ" &
                                         n_times_id == 4 &
                                         n_obs_id == 72)), 1)
testJob(id = id)

# Difficult: Always td2pLL, no time-dep, small n, large noise

id = head(findExperiments(algo.name = "no_pre",
                          algo.pars = (always == "td2pLL"),
                          prob.pars = (model_id == "M3" &
                                         noise_id == "N3" &
                                         spacing_id == "equ" &
                                         n_times_id == 4 &
                                         n_obs_id == 72)), 1)
testJob(id = id)

# In case there are errors: Find the corresponding simulation scenario's 
#   factor levels

# id = head(findErrors(), 1)
# getJobPars(id = head(findErrors(), 1))

#####################################
# Chunk and submit
######################################

# Put the many jobs into chunks to speed up the computing

all_jobs <- findJobs()

# for minimal simulation, remove following row that chunks the jobs
all_jobs[,chunk:=chunk(all_jobs$job.id, n.chunks = 288)] # 54*3


submitJobs(all_jobs, resources = list(walltime = 2500L, memory = 2048))

# check if jobs are done
getStatus()


################
# Save results
################

# In case you return to the session on a HPC cluster, you may need something like this
# loadRegistry(file.dir = "./04_simulation_run/registry/",
# writeable = TRUE, conf.file = "/work/smjududa/.batchtools.conf.R")

# gather results as are (big list)
results_simulation_raw <- reduceResultsList()
save(results_simulation_raw, file = "./04_simulation_run/results_simulation_raw.RData")

# generate more handy format (data.frame) for results
add_res <- lapply(results_simulation_raw, function(x) unlist(x[1:4])) %>%
  do.call(rbind, .) %>% as.data.frame()
colnames(add_res) <- c("MAFC", "MAcc", "signif", "is_conv")

# save handy results
results_simulation <- cbind.data.frame(unwrap(getJobPars()), add_res)
save(results_simulation, file = "./04_simulation_run/results_simulation.RData")

