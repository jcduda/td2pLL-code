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



# function: plot_designs
#           plots the considered experimental designs as a ggplot. 
#
# Input:
# - designs:  data.frame with columns expo, dose, n and design
#
# Details:
#   Right now, inceased doses in a log scale with base sqrt(10) is assumed
#   for plotting the dose-axes accordingly.

plot_designs <- function(designs){

    
    l_br <- floor(log10(min(designs$dose[designs$dose > 0])))
    up_br <- round(log10(max(designs$dose)))
    
    break_doses <- c(0, 10^(l_br:up_br))
    
    
    ggplot(data = designs, aes(x = dose, y = time, size = factor(n))) +
      geom_point(alpha = 0.7, color = "steelblue") +
      geom_text(aes(label=n), size = 3, hjust=+0.5, vjust=0.5) +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = break_doses[2]/2, base = sqrt(10)),
                         breaks = break_doses,
                         labels = round(break_doses, 4)) +
      labs(y = "exposure time") +
      guides(size = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(facets = "design", labeller = "label_both")
  
}



# function: my_pseudo_log
#           Creates a pseudo_log of the doses, so that dose=0
#           does not hrow an error
# Input:
# - x:  positive numeric
my_pseudo_log <- function(x) asinh(x/(2 * 0.0001))/log(sqrt(10))



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

