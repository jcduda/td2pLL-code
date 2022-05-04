

####################################
####################################
#                                  #
# pre-process data of Gu et al.    #
#                                  #
####################################
####################################

library(readxl)
library(dplyr)
library(tidyr)
library(DoseFinding)
library(ggplot2) 

############################

# find all compound name using the filenames of the "Day1"-Experiments:
compounds <- list.files(path = "./01_data_raw/Paket13Day1/") %>%
  gsub("PHH_", "", .) %>%
  gsub("_.*", "", .)

# find all directory names for the exposure times
time_dirs <- list.files("./01_data_raw")

# considered exposure times (here in days)
times <- c(1, 2, 7)



############################


# function: prep_single_compound
#
#   Generates a long format data.frame to store all data of one compound 
#
# Input:
#   - compound: character for the compound name
#   - time_dirs: directory names of the exposure time files (one per exposure time)
#   - times:    numeric vector of actually used exposure times. Here: c(1, 2, 7)
#
# Output:
#   - data.frame (long format) with all data for one compound with columns:
#       - compound
#       - Donor
#       - time
#       - dose
#       - resp (= raw_resp / control_mean)
#       - sampleID
#       - Donor_name
#       - control_mean (mean resp value for given time and Donor)
#       - raw_resp
#       - by_control
#           - 0: only one control was used to normalize all raw_resp
#           - 1: there were control_1 and control_2, this value was normalized
#                 with control_1
#           - 2: there were control_1 and control_2, this value was normalized
#                 with control_2

prep_single_compound <- function(compound, time_dirs, times){
  
  print(compound)
  
  # get the three paths for the compound belonging to the three exposure time levels
  for(i in 1:length(times)){
    assign(paste0("d",i,"_path"),
           grep(compound, list.files(paste0("./01_data_raw/", time_dirs[i]),
                                     full.names = T), value = T)
    )
  }
  
  # function: adjust_control
  #
  #   Pre-processes wide-format tmp_dat: generates response that is divided by 
  #   mean control response of this donor.
  #
  # Input:
  #   - tmp_dat: data.frame (wide format) of current compound with columns
  #     "SampleID", "Donor", then either 
  #        - "control" followed by some concentration levels or 
  #        - "control_1", some concentration levels, "control_2", some 
  #           concentration levels, then
  #     "Donor_name", "time"
  #
  #   - nof_control: 1 or 2 (numeric), the number of control-columns
  #
  # Output:
  #   - data.frame (long format) for this compound with columns
  #       - SampleID
  #       - Donor
  #       - Donor_name
  #       - time
  #       - control_mean (mean raw_resp at dose 0 for given time and Donor)
  #       - dose
  #       - raw_resp (raw response)
  #       - resp (response, i.e. raw_resp / control_mean)
  
  adjust_control <- function(tmp_dat, nof_control){
    # only one control
    if(nof_control == 1){
      stopifnot(sum(grepl("control", colnames(tmp_dat))) == 1)
      tmp_dat <- tmp_dat %>% 
        group_by(time, Donor) %>%
        mutate(control_mean = mean(control)) %>%
        ungroup() %>%
        rename("0mM" = "control") %>%
        pivot_longer(cols = contains("mM"), names_to = "dose",
                     values_to = "raw_resp") %>%
        mutate(dose = as.numeric(gsub("mM|M", "", dose))) %>%
        mutate(resp = (raw_resp / control_mean) * 100) 
      return(tmp_dat)
    }
    
    # when two controls are used for different dose levels
    if(nof_control == 2){
      # everything that used control_1
      part1 <- tmp_dat %>% dplyr::select(SampleID:control_2, Donor_name:time) %>%
        dplyr::select(-control_2) %>%
        rename("control" = "control_1") 
      part1 <- adjust_control(tmp_dat = part1, nof_control = 1) %>%
        mutate(by_control = 1)
      
      # everything that used control_2
      part2 <- tmp_dat %>% dplyr::select(SampleID:Donor, control_2:time) %>%
        rename("control" = "control_2")
      part2 <- adjust_control(tmp_dat = part2, nof_control = 1) %>%
        mutate(by_control = 2)
      
      # put together and use by_control as a flag
      tmp_dat <- rbind.data.frame(part1, part2)
      return(tmp_dat)
    }
  }
  
  # get the three files for the substance and do the processing into a
  # long format with above function
  for(i in 1:length(times)){
    tmp_dat <- suppressMessages(read_excel(path = get(paste0("d",i,"_path"))) %>%
      dplyr::select_if(~sum(!is.na(.)) > 0) %>%
      mutate(time = times[i]))
    
    
    if("control_1" %in% colnames(tmp_dat) & !("control_2" %in% colnames(tmp_dat))){
      tmp_dat <- tmp_dat %>% rename("control" = "control_1")
    }
    
    which_col <- which(colnames(tmp_dat) == "time") - 1 # contains Donor name
    
    colnames(tmp_dat)[which_col] <- "Donor_name"
    
    # assign(paste0("data.time.",expos[i]), tmp_dat)
    
    # typo in VITC
    if("0mM" %in% colnames(tmp_dat)){
      tmp_dat <- rename(tmp_dat, "control" = "0mM")
    }
    
    # check if there is only one control or control_1 and control_2,
    # then adjust for control accordingly
    nof_control <- sum(grepl("control", colnames(tmp_dat)))
    if(!(nof_control %in% c(1,2))){
      stop(paste0("Compound ", compound," not used, because not 1 or 2 controls. \n"))
    }
    
    # typo in LEV
    if("0.316M" %in% colnames(tmp_dat)){
      tmp_dat <- rename(tmp_dat, "0.316mM" = "0.316M")
    }
    
    # typo in GLC
    if("316µM" %in% colnames(tmp_dat)){
      tmp_dat <- rename(tmp_dat, "0.316mM" = "316µM")
    }

    # bigger processing step
    prep_data <- adjust_control(tmp_dat = tmp_dat, nof_control = nof_control) 
    if(!("by_control" %in% colnames(prep_data))) prep_data <- prep_data %>% mutate(by_control = 0)
    
    assign(paste0("prep_data_", i), prep_data)
  }
  
  # some final arrangements
  res <- rbind.data.frame(prep_data_1, prep_data_2, prep_data_3) %>%   
    mutate_at(vars(SampleID, Donor, Donor_name, by_control), as.factor) %>% 
    arrange(Donor, time, dose) %>%
    mutate(compound = compound) %>%
    dplyr::select(compound, Donor, time, dose, resp, everything())
  
  return(as.data.frame(res))
  
}


#####################################################################
#####################################################################

# Apply above pre-processing function to all compounds
# to generate a neat, long format data.frame


data_prep <- lapply(compounds, function(compound){
  
  tryCatch({
    prep_single_compound(compound, time_dirs = time_dirs, times = times)
  },
  error = function(cond){
    message(cond)
  })
}
) %>% 
   do.call(rbind.data.frame, .) 

#######################################################################

# Rows with (originally) missing values will be deleted

# Take a look at which rows are of concern,if desired:
# data_prep %>% filter(is.na(resp)) 

# Delete rows with missing resp
data_prep <- data_prep %>% filter(!is.na(resp)) 

save(data_prep, file = "./02_data_prep/data_prep.RData")


######################################################################
######################################################################

# Additional pre-processing: refitting-algorithm 

# The alreasy pre-processed data (put into long-format) will go thourgh
# a secondary, optional pre-processing step, yielding to the
# "data_refit_prep" file.
#
# For the entire analysis, the pre-processed data which includes this step
# is always used. 
#
# The above generated data_prep is for possible, individual analysis purposes.
#
# The refitting algorithm does the follwing:
#
# For a given compound, exposure time and donor, a dose-response curve (4pLL) 
# is fitted.
# The corresponding data is divided through the resulting left (upper) asymptote
# and again multiplied by 100.
# This way, the data is expected to better follow the assumption of having a 
# left asymptote at 100 [percent].

data_refit_prep <- data_prep

data_refit_prep$left_asymp <- NA
data_refit_prep$resp_refit <- NA  


for(curr_compound in compounds){
  tmp_compound <- data_refit_prep %>% filter(compound == curr_compound)
  times <- tmp_compound$time %>% unique
  
  for(curr_time in times){
    tmp_compound_time <- tmp_compound %>% filter(time == curr_time)
    donors <- tmp_compound_time$Donor %>% unique
    
    for(curr_donor in donors){
      # fit 4pLL (aka sigEmax) model using DoseFinding package
      suppressMessages(
      left_asymp <- fitMod(dose, resp,
                           data = tmp_compound_time %>%
                             filter(Donor == curr_donor) %>%
                             dplyr::select(dose, resp),
                           model = "sigEmax") %>%
        coef() %>% 
        .["e0"]
      )
      
      data_refit_prep[data_refit_prep$compound == curr_compound &
                data_refit_prep$time == curr_time &
                data_refit_prep$Donor == curr_donor, "left_asymp"] <- left_asymp
      
    }
  }
}

# add the new, adjusted repsonse: (old response / left_asymptote) * 100

data_refit_prep$resp_refit <- (data_refit_prep$resp / data_refit_prep$left_asymp) * 100

# Visually check if the control response level is now closer to 100


# data_refit_prep %>% filter(dose == 0) %>%
#   select(resp, resp_refit) %>%
#   pivot_longer(cols = c(1, 2)) %>%
#   ggplot(., aes(x= value)) +
#   geom_histogram() +
#   facet_wrap(facets = "name")

# Actually does not look like it helped.

# overwrite
data_refit_prep$resp <- data_refit_prep$resp_refit

save(data_refit_prep, file="./02_data_prep/data_refit_prep.RData")











