#######################################
# functions for the td2pLL simulation #
#######################################

####################
# Required packages
####################

library(dplyr)
library(DoseFinding)
library(ICAOD)
library(tidyr)
library(plotly)
library(drc)

# To install the td2pLL package:
library(devtools)
# devtools::install_github("jcduda/td2pLL")
library(td2pLL)


#####################
# Required functions
#####################

# This function is not included in the td2pLL package. It is tailored
# to the simulations study

# Function fit_sep_2pLL
#
# Fit separate 2pLL curves to a time-dose-response data set where the
# upper asymptote is fixed to 100 and the lower to 0.
# The functions uses the drc::drm() function.
# As numerical algorithm, the function first uses Nelder-Mead.
# If this fails, it uses BFGS.
#
# Input:
# data - data.frame with numerical variables dose, time, resp.
#
# Output: A drm-object generated with drc::drm()



fit_sep_2pLL <- function(data){
  tryCatch(
    {
      drm(resp ~ dose,
          curveid = time,
          data = data %>% dplyr::mutate(time = factor(time)),
          control = drmc(method = "Nelder-Mead"),
          fct = LL2.2(upper = 100),
          pmodels = list(~time, # h
                         ~time) # EC50
      )
    },
    error = function(cond){
      drm(resp ~ dose,
          curveid = time,
          data = data %>% dplyr::mutate(time = factor(time)),
          fct = LL2.2(upper = 100),
          pmodels = list(~time, # h
                         ~time) # EC50
      )
    }
  )
}



#######################
# Data generation 
#######################



# function: generate_data
#           generates data of an extended 3pLL model for specified noise and experimental design (expDes)
# Input:
# - model:    named, numeric character with model parameters h, delta, gamma and c0
# - noise:    character "N1", "N2, or "N3" for low, baseline or high noise level respectively
#             Details: For the baseline N2, we use a linear model that uses dose and time as covariates,
#             which has to be tranformed via my_pseudo_log first, to define the variance of the mean-zero, normal distributed
#             noise at each dose-expo point of expDes. For N1, the variance is halfed. For N3, it is multiplied by 1.5.
# - expDes:   data.frame with columns expo, dose and n where n specifies how many replicates are measured at given expo-dose point
#
# Output:
# - res:      data.frame containing generated responses. Columns are expo, dose, resp

generate_data <- function(model, noise_id, expDes){
  # grid for time, dose, n, h, delta, gamma, c0
  inputs <- cbind.data.frame(expDes %>% dplyr::select(-design),
                             t(as.data.frame(model)) %>%
                               `rownames<-`(NULL))
  
  # calculate mean_resp (true mean) at each time  dose point
  inputs$mean_resp <- apply(inputs, 1, function(x){
    td2pLL::td2pLL(time = as.numeric(x["time"]),
                   dose = as.numeric(x["dose"]),
                   h = as.numeric(x["h"]),
                   gamma = as.numeric(x["gamma"]),
                   c0 = as.numeric(x["c0"]),
                   delta = as.numeric(x["delta"]))
  } 
  )
  
  # Calculate noise sd first
  inputs$noise_sd <- predict(lm_sd_noise,
                             newdata = inputs[, c("time", "dose")] %>%
                               dplyr::mutate(dose = my_pseudo_log(dose)))
  
  # check via noise parameter, if sd will be halfed (N1), or doubled (N3)
  
  stopifnot(noise_id %in% c("N1", "N2", "N3"))
  
  if(noise_id == "N1"){
    inputs$noise_sd <- inputs$noise_sd * 0.5
  }
  
  if(noise_id == "N3"){
    inputs$noise_sd <- inputs$noise_sd * 2
  }
  
  # if noise_id == "N2": Do nothing, this is the baseline noise sd.
  
  
  # generate random noise values
  help_add_noise <- apply(inputs[, c("time", "dose", "n", "noise_sd")], 1, function(x){
    suppressWarnings(
      cbind.data.frame(time = x["time"], dose = x["dose"],
                       noise_value = rnorm(n = x["n"], sd = x["noise_sd"]))
    )
  }) %>% do.call(rbind.data.frame, .)
  
  # put together
  res_pre <- left_join(inputs, help_add_noise, by = c("time", "dose")) %>%
    dplyr::mutate(resp = mean_resp + noise_value) %>%
    dplyr::select(time, dose, resp)
  
  # apply regular pre-processing
  
  # First, divide by mean response at dose 0
  
  mean_resp_0 <- res_pre %>% dplyr::filter(dose == 0) %>% dplyr::select(resp) %>% unlist %>% mean
  res_pre <- dplyr::mutate(res_pre, resp = (resp / mean_resp_0) * 100 )
  
  # refit procedure: seperately fit 4pLL dose-response curves at each exposure time.
  # Divide by resulting left (upper) asymptote
  all_times <- unique(res_pre$time)
  
  # get left asymptote (e0) of each exposure time
  help_left_asymp <- sapply(all_times, function(curr_time){
    c(time = curr_time,
      fitMod(dose, resp,
             data = res_pre %>% dplyr::filter(time == curr_time) %>% dplyr::select(dose, resp),
             model = "sigEmax") %>%
        coef() %>% 
        .["e0"]
    )
  }) %>% t() %>% as.data.frame()
  
  # divide by left asymptote stratified by exposure time
  res <- left_join(res_pre, help_left_asymp, by = "time") %>%
    dplyr::mutate(resp = (resp / e0) * 100) %>%
    dplyr::select(time, dose, resp)
  
  
  return(res)
  
}


#############################################################################
# Helper functions for plotting results in 05_results/01_generate_figures.R
#############################################################################

# function: get_2pll_at_t (a small helper function)
#
# From a td2pLL model specified through "coef",
# get the response values at "dose" and a fixed time "time"

get_2pll_at_t <- function(dose, time, coefs){
  ED50 <- get_ED50s(coefs = coefs, times = time)$ED50
  res <- sigEmax(dose, e0 = 100, eMax = -100, ed50 = ED50, h = coefs["h"])
  return(unname(res))
}


# function: my_plot_designs (small helper function)
#
# Visualizes the time-dose design with a ggplot.
#
#
# Input: design (data.frame) which looks like this:
#
#   time   dose     n design       n_obs spacing n_times              pt_size
#   1     1 0         12 D_equ_3_n216   216 equ     3 exposure durations       5
#   2     1 0.01      12 D_equ_3_n216   216 equ     3 exposure durations       5
#   3     1 0.0316    12 D_equ_3_n216   216 equ     3 exposure durations       5
#   4     1 0.1       12 D_equ_3_n216   216 equ     3 exposure durations       5
#                               ...

my_plot_designs <- function(designs){
  
  
  l_br <- floor(log10(min(designs$dose[designs$dose > 0])))
  up_br <- round(log10(max(designs$dose)))
  
  break_doses <- c(0, 10^(l_br:up_br))
  
  
  ggplot(data = designs, aes(x = dose, y = time, size = pt_size)) +
    geom_point(alpha = 0.7, size = designs$pt_size, color = "steelblue") +
    geom_text(aes(label=n), size = 3, hjust=+0.55, vjust=0.4) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = break_doses[2]/2, base = sqrt(10)),
                       breaks = break_doses,
                       labels = round(break_doses, 4)) +
    scale_y_continuous(breaks = c(1, 2, 4, 7)) +
    labs(y = "Exposure duration", x = "Concentration") +
    guides(size = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(facets = "n_times")
  
}


# function plot_td2pLL_2dim 
#
# Plots a td2pLL_mod object on a regular plot with only a dose-axis and a response-axis.
# The time-dimension is represented by seperately fitting the marginal dose-response curves at
# times 1, 2 and 7.
# This function is not general and only tailored to the situation of the simulation results.

plot_td2pLL_2dim <- function(td2pLL_mod, sigma_scal = sqrt(10)){
  
  params <- coef(td2pLL_mod)
  data <- td2pLL_mod$orig_data %>%
    group_by(time, dose) %>%
    summarize(resp = mean(resp)) %>%
    mutate(Time = as.factor(time))
  break_doses <- data$dose %>% unique
  ED50s <- get_ED50s(coefs = params, times = c(1, 2, 7))[,2]
  
  my_2pLL_1 <- function(x){
    100 - 100 * x^params["h"] / (ED50s[1]^params["h"] + x^params["h"])
  }
  
  my_2pLL_2 <- function(x){
    100 - 100 * x^params["h"] / (ED50s[2]^params["h"] + x^params["h"])
  }
  
  my_2pLL_3 <- function(x){
    100 - 100 * x^params["h"] / (ED50s[3]^params["h"] + x^params["h"])
  }
  
  p <- ggplot(data = data,
              aes(x = dose, y = resp, shape = Time, linetype = Time)) +
    geom_point() +
    scale_shape_manual(values = c(1, 2, 3)) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = break_doses[2]/sigma_scal,
                                                        base = sqrt(10)),
                       breaks = break_doses,
                       labels = round(break_doses, 4)) +
    stat_function(fun = my_2pLL_1) +
    stat_function(fun = my_2pLL_2, linetype = "dashed") +
    stat_function(fun = my_2pLL_3, linetype = "dotted") +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2, 3),
                                                    linetype = c(1, 2, 3)))) +
    labs(x = "Concentration", y = "Response") +
    theme(legend.position = c(0.9, 0.9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  return(p)
  
}


# function plot_2pLL
#
# Plots a 2pLL model fitted with drc::drm(., curveid = time) in a customized way
# for visualizing the results.
# This function is not general but customized to for the figures of the manuscript.

plot_2pLL <- function(LL2_fit, sigma_scal = sqrt(10)){
  params <- coef(LL2_fit)
  
  data <- LL2_fit$origData %>%
    group_by(time, dose) %>%
    summarize(resp = mean(resp)) %>%
    mutate(Time = as.factor(time))
  
  break_doses <- data$dose %>% unique
  
  ED50_1 <- exp(params["e:(Intercept)"])
  ED50_2 <- exp(params["e:(Intercept)"] + params["e:time2"])
  ED50_3 <- exp(params["e:(Intercept)"] + params["e:time7"])
  
  h_1 <- params["b:(Intercept)"]
  h_2 <- params["b:(Intercept)"] + params["b:time2"]
  h_3 <- params["b:(Intercept)"] + params["b:time7"]
  
  my_2pLL_1 <- function(x){
    100 - 100 * x^h_1 / 
      (ED50_1^h_1 + x^h_1)
  }
  
  my_2pLL_2 <- function(x){
    100 - 100 * x^h_2 / 
      (ED50_2^h_2 + x^h_2)
  }
  
  my_2pLL_3 <- function(x){
    100 - 100 * x^h_3 / 
      (ED50_3^h_3 + x^h_3)
  }
  
  p <- ggplot(data = data,
              aes(x = dose, y = resp, shape = Time, linetype = Time)) +
    geom_point() +
    scale_shape_manual(values = c(1, 2, 3)) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = break_doses[2]/sigma_scal,
                                                        base = sqrt(10)),
                       breaks = break_doses,
                       labels = round(break_doses, 4)) +
    stat_function(fun = my_2pLL_1) +
    stat_function(fun = my_2pLL_2, linetype = "dashed") +
    stat_function(fun = my_2pLL_3, linetype = "dotted") +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2, 3),
                                                    linetype = c(1, 2, 3)))) +
    labs(x = "Concentration", y = "Response") +
    theme(legend.position = c(0.9, 0.9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  return(p)
  
  
}
