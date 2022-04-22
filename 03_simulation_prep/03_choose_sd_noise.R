#################################
# Choose noise levels for doses #
#################################

# When generating data using the models (M1) - (M3) and the
# experimental designs (P1) - (P3b), it remains to decide
# on the noise that will be added to the data.


# Assumption: We assume normal, mean-zero noise. Hence,  only
# the standard veiation sd for the noise must be specified.
# This is done here, based on the data of Gu et al.

library(ggplot2)
library(dplyr)

source("./00_functions/functions.R")
load("./02_data_prep/data_refit_prep.RData")

data <- data_refit_prep %>% dplyr::select(compound, Donor, dose, time, resp)

# By compound, Donor, time and dose: Calculate sd

data <- data %>% dplyr::group_by(by = compound, Donor) %>%
  mutate(maxD = max(dose)) %>%
  mutate(dose = dose / maxD) %>%
  ungroup %>% 
  group_by(compound, Donor, time, dose) %>%
  mutate(sd = sqrt(var(resp)))

# Plot the sd values:
# Within each exposre time seperately, fit a linear model with
# stat_smooth

plot_sd_values <- ggplot(data, aes(x = round(dose, 4), y = sd)) +
  geom_jitter(width = 0.2, alpha = 0.05) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  theme_bw() +
  labs(x = expression(dose~x[j]), y = expression(sigma(t[i],x[j]))) +
  coord_cartesian(ylim = c(0, 65)) +
  # facet_wrap(facets = "time", labeller = "label_both") +
  facet_wrap(facets = "time", labeller = label_bquote(t[i] == .(time))) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.0001, base = sqrt(10)),
                     breaks = sort(unique(round(data$dose, 3))),
                     labels = sort(unique(round(data$dose, 3)))
  )

# plot_sd_values

saveRDS(plot_sd_values, file = "./03_simulation_prep/03_plot_sd_values.rds")
save(plot_sd_values, file = "./03_simulation_prep/03_plot_sd_values.RData")
# apparently, the noise rises at first and then declines again,
# especially at dose = 1
# the overall mean of sd seems ot increase with time increasing

##############################
# create linear model for sd #
##############################


lm_sd_noise <- lm(sd ~ dose +
                    I(dose^2)+
                     I(dose^2*time) + 
                    time,
                  data = data %>% mutate(dose = my_pseudo_log(dose)))

# summary(lm_sd_noise)

saveRDS(lm_sd_noise, file = "./03_simulation_prep/03_lm_sd_noise.rds")
save(lm_sd_noise, file = "./03_simulation_prep/03_lm_sd_noise.RData")
###################
# Look at the model
###################

seq_data <- expand.grid(time = c(1, 2, 7),
                        dose = c(0, exp(seq(log(1e-4), log(1), length.out = 100)))
)
seq_data<- seq_data %>% mutate(sd_pred = predict(lm_sd_noise,  newdata = seq_data %>%
                                                   mutate(dose = my_pseudo_log(dose))))


# Take a look:

ggplot(data = seq_data, aes(x = my_pseudo_log(dose), y = sd_pred)) +
  geom_line(color = "steelblue", size = 2) +
  facet_wrap(facets = "time", labeller = "label_both") +
  labs(x = "(pseudo) log(dose)", y = "predicted sd") +
  coord_cartesian(ylim = c(0, 65)) +
  theme_bw()


####################

# For the data generation, we will consider three factors of noise:
# N1, N2, N3, where 
#
# N1 is half of what the sd-model gives (clear data)
# N2 is what the sd model gives (baseline)
# N3 is 2 times what the sd model gives (noisy data)



