
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggplot2); theme_set(theme_bw())
library(td2pLL)
library(DoseFinding)

setwd("C:/Users/duda/Projekte/exp_dose_resp_code")

source("./00_functions/functions.R")
load("./02_data_prep/data_refit_prep.RData")
load("./03_simulation_prep/02_true_designs.RData")
load("./03_simulation_prep/03_plot_sd_values.RData")
load("./03_simulation_prep/03_lm_sd_noise.RData")
load("./04_simulation_run/results_simulation.RData")






# TODO:
# Clean up: figure by figure
#           - generate ennumerated figure plots
# Clean up 00_functions/functions. Eventually add individual functions used here
# to keep this script cleaner.
#
# Clean up 02_data_prer: Do we anywhere need .rds objects? Maybe clean up the README there, too.


############
# Figure 1 #
############

# Note: The figures appear in the Viewer and were manually exported and saved.

# 1a
plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 0.1, gamma = 2, c0 = 0.01),
            dose_lim = c(0.001, 1), time_lim = c(1, 5), 
            xaxis_title = "Concentration", yaxis_title = "Exposure duration")
# 1b
plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 0.1, gamma = 2, c0 = 0.1),
            dose_lim = c(0.001, 1), time_lim = c(1, 5), 
            xaxis_title = "Concentration", yaxis_title = "Exposure duration")
# 1c
plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 0.1, gamma = 4, c0 = 0.01),
            dose_lim = c(0.001, 1), time_lim = c(1, 5), 
            xaxis_title = "Concentration", yaxis_title = "Exposure duration")

# 1d
plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 0.1, gamma = 4, c0 = 0.1),
            dose_lim = c(0.001, 1), time_lim = c(1, 5), 
            xaxis_title = "Concentration", yaxis_title = "Exposure duration")


#############
## Figure 2
#############


# function: get_2pll_at_t
#
# --> a small helper function
#
# From a td2pLL model specified through "coef",
# get the response values at "dose" and a fixed time "time"

get_2pll_at_t <- function(dose, time, coefs){
  ED50 <- get_ED50s(coefs = coefs, times = time)$ED50
  res <- sigEmax(dose, e0 = 100, eMax = -100, ed50 = ED50, h = coefs["h"])
  return(unname(res))
}


ggplot(data = NULL) +
  
  stat_function(fun = get_2pll_at_t, args =
                  list(coefs = c(h = 2, delta = 0.1, gamma= 4, c0 = 0.01), time = 1),
                aes(colour = "t = 1")) +
  
  stat_function(fun = get_2pll_at_t, args =
                  list(coefs = c(h = 2, delta = 0.1, gamma= 4, c0 = 0.01), time = 2),
                aes(colour = "t = 2")) +
  
  stat_function(fun = get_2pll_at_t, args =
                  list(coefs = c(h = 2, delta = 0.1, gamma= 4, c0 = 0.01), time = 3),
                aes(colour = "t = 3")) +
  
  stat_function(fun = get_2pll_at_t, args =
                  list(coefs = c(h = 2, delta = 0.1, gamma= 4, c0 = 0.01), time = 100),
                aes(colour = "t = 100")) +
  
  geom_hline(yintercept = 50, linetype = "dashed")+
  geom_point(data = data.frame(dose = get_ED50s(coefs = c(h = 2, delta = 0.1, gamma= 4, c0 = 0.01),
                                                times = c(1, 2, 3, 100))$ED50,
                               resp = 50), 
             mapping = aes(x = dose, y = resp), size = 2, color = "black") +
  
  annotate("text", x = 0.01, y = 0, label = expression(C[0]==0.01), size = 5) +
  
  geom_segment(aes(x = 0.01, y = 37.5, xend = 0.11, yend = 37.5),
               arrow = arrow(length = unit(0.3, "cm"), ends = "both")) +
  geom_segment(aes(x = 0.01, y = 37.5, xend = 0.01, yend = 50),
               linetype = "dashed") +
  geom_segment(aes(x = 0.11, y = 37.5, xend = 0.11, yend = 50),
               linetype = "dashed") +
  
  annotate("text", x = 0.04, y = 42, label = expression(Delta==0.1), size = 5) +
  
  coord_cartesian(xlim = c(0.001, 1), ylim = c(0, 100)) +
  scale_x_continuous(trans = "log10", limits = c(0.0001, 1)) +
  
  labs(x = "Concentration", y = "Response", colour = "Exp. dur.") +
  scale_color_discrete(breaks=c("t = 1","t = 2","t = 3","t = 100")) +
  theme(legend.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.position = c(0.85, 0.75),
        legend.text = element_text(size = rel(1.3), margin = margin(t = 1)),
        plot.margin = margin(0, 10, 0, 1)) +
  guides(colour = guide_legend(override.aes = list(size = unit(2, "cm"))))

ggsave(filename = "./05_results/fig_2_explain_model_2.tiff",
       width = 5, height = 4)
ggsave(filename = "./05_results/fig_2_explain_model_2.pdf",
       width = 5, height = 4)

###############################################################################
###############################################################################

# STOPPED HERE

# Export manually into jpegs


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

my_plot_designs(designs = designs %>% filter(n_obs == 216,
                                             spacing == "equ") %>%
                  mutate(pt_size = ifelse(n == 12, 5, 4),
                         n_times = paste0(n_times," exposure durations")))

ggsave("fig_expDes_sim.pdf", width = 5, height = 2.75, path = "publication_progress/Figures/")


load("./03_simulation_prep/03_plot_sd_values.RData")
load("./03_simulation_prep/03_lm_sd_noise.RData")
plot_sd_values +
  labs(x = expression(Concentration~x[j]))

ggsave("fig_sd_sim.pdf", width = 7.5, height = 2.75, path = "publication_progress/Figures/")


################################################################################
################################################################################


data_subset <- data_refit_prep %>% filter(compound == "ASP") %>% 
  #filter(Donor == "Don_3") %>%
  dplyr::select(time, dose, resp)

fit <- fit_td2pLL(data = data_subset)
plot(fit, add_data = data_subset, xaxis_title = "Concentration",
     yaxis_title = "Exposure duration", zaxis_title = "Response")



data(cytotox)

all_compounds <- cytotox$compound %>% unique
app_res <- data.frame(compound = all_compounds,
                      td2pLL_rse = NA,
                      sep_rse = NA,
                      joint_rse = NA,
                      anova_signif = NA)

model_list <- myList <- vector("list", length(all_compounds))
names(model_list) <- all_compounds

pdf("publication_progress/Figures/fig_all_compounds.pdf", width = 10, height = 12)
par(mfrow = c(6, 5), mar = c(3, 3, 1, 0.5), mgp = c(1.75, 0.5, 0))
i <- 1
set.seed(1905)
for(i in 1:nrow(app_res)){
  print(i)
  curr_comp <- app_res$compound[i]
  
  data_subset <- cytotox[cytotox$compound == curr_comp, c("expo", "dose", "resp")]
  colnames(data_subset)[1] <- "time"
  
  # td2pLL
  m_fit_td2pLL <- fit_td2pLL(data = data_subset)
  td2pLL_predict <- predict(m_fit_td2pLL, newdata = data_subset)
  td2pLL_rse <- sqrt(sum((td2pLL_predict - data_subset$resp)^2) /
                       (nrow(data_subset) - length(coef(m_fit_td2pLL))))
  
  # seperate
  m_fit_sep <- tryCatch(
    {fit_sep_2pLL(data_subset)},
    error = function(cond) return(NA))
  
  if(length(m_fit_sep) == 1){
    sep_rse <- NA 
  } else {
    sep_predict <- predict(m_fit_sep, newdata = data_subset)
    sep_rse <- sqrt(sum((sep_predict - data_subset$resp)^2) /
                      (nrow(data_subset) - length(coef(m_fit_sep))))
  }
  plot(m_fit_sep, main = all_compounds[i], xlab = "Concentration [mg]", ylab = "Viability [%]",
       broken = TRUE, 
       legend = ifelse(i==1, TRUE, FALSE))
  
  # joint
  m_fit_joint <- td2pLL:::fit_joint_2pLL(data_subset)
  joint_predict <- predict(m_fit_joint, newdata = data_subset)
  joint_rse <- sqrt(sum((joint_predict - data_subset$resp)^2) /
                      (nrow(data_subset) - length(coef(m_fit_joint))))
  
  
  
  app_res$td2pLL_rse[i] <- td2pLL_rse
  app_res$sep_rse[i] <- sep_rse
  app_res$joint_rse[i] <- joint_rse
  app_res$anova_signif[i] <- as.numeric(td2pLL_anova(data = data_subset, alpha = 0.05)$signif)
  model_list[[i]] <- list(m_fit_td2pLL = m_fit_td2pLL,
                          m_fit_sep = m_fit_sep)
}
dev.off()

save(model_list, file = "./publication_progress/Figures/model_list.RData")

# td2pLL vs sep
app_res %>%
  # Remove compounds where computational issues for separate 2pLL fitting arised
  filter(!compound %in% c("BUSF", "MePA", "RIF", "FAM")) %>%
  ggplot(aes(x = td2pLL_rse, y = sep_rse,
             #col = as.factor(anova_signif),
             label = compound)) +
  geom_point(size = 3) +
  geom_text_repel(alpha = 0.8, color = "grey", max.overlaps = 100) +
  geom_abline(intercept = 0, slope = 1) +
  labs(color = "ANOVA signif.", x = "New td2pLL model: Deviation model fit to data",
       y = "Separate 2pLL models: Deviation model fit to data") +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,5,0),
        legend.box.margin=margin(-10,-10,-10,-10))

ggsave("fig_td2pLL_app.pdf", width = 4, height = 4,  path = "publication_progress/Figures/")


################################################
# Compare td2pLL and separate fit
################################################

# Choose CHL
td2pLL_mod <- model_list$PPL$m_fit_td2pLL
plot_td2pLL_2dim(td2pLL_mod)
# plot a td2pLL_fit in a two-dimensional manner

#############################################################
# plot td2pLL in normal, two-dimensional plot
#############################################################

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

save(plot_td2pLL_2dim, file = "./publication_progress/Figures/fct_plot_td2pLL_2dim.RData")

compound <- "PPL"
LL2_fit <- model_list[[compound]]$m_fit_sep
plot_2pLL(LL2_fit)
# plot a separated (with curveid = time) drm model fit of the LL2.2

#############################################################
# Custom plot 2pLL in normal, two-dimensional ggplot
#############################################################

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

save(plot_2pLL, file = "./publication_progress/Figures/fct_plot_2pLL.RData")


compound <- "CHL"

for(compound in all_compounds) {
  print(plot_2pLL(model_list[[compound]]$m_fit_sep)  +
    labs(title = paste0(compound, " (Sep. fit)")))
  print(plot_td2pLL_2dim(model_list[[compound]]$m_fit_td2pLL) +
    labs(title = paste0(compound, " (td2pLL. fit)")))
}

# Selections:
# PPL, GLC: Separate better
# LAB: td2pLL looks more stable
# FAM: No computational problems in td2pLL
# DFN: td2pLL more stable as ED50 must decrease with time 

require(gridExtra)

## PPL ##
p11 <- plot_2pLL(model_list[["PPL"]]$m_fit_sep, sigma_scal = 100)  +
  labs(title = "PPL") +
  theme(legend.position = c(0.2, 0.4),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

p12 <- plot_td2pLL_2dim(model_list[["PPL"]]$m_fit_td2pLL, sigma_scal = 100) +
  labs(title = "") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

## GLC ##
p21 <- plot_2pLL(model_list[["GLC"]]$m_fit_sep, sigma_scal = 10)  +
  labs(title = "GLC") +
  theme(legend.position = "none",
        axis.title.x = element_blank())

p22 <- plot_td2pLL_2dim(model_list[["GLC"]]$m_fit_td2pLL, sigma_scal = 10) +
  labs(title = "") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

## LAB ##
p31 <- plot_2pLL(model_list[["LAB"]]$m_fit_sep, sigma_scal = 10)  +
  labs(title = "LAB") +
  theme(legend.position = "none",
        axis.title.x = element_blank())

p32 <- plot_td2pLL_2dim(model_list[["LAB"]]$m_fit_td2pLL, sigma_scal = 10) +
  labs(title = "") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

## DFN ##
p41 <- plot_2pLL(model_list[["DFN"]]$m_fit_sep, sigma_scal = 10)  +
  labs(title = "DFN") +
  theme(legend.position = "DFN")

p42 <- plot_td2pLL_2dim(model_list[["DFN"]]$m_fit_td2pLL, sigma_scal = 10) +
  labs(title = "") +
  theme(legend.position = "none",
        axis.title.y = element_blank())


grid.arrange(p11, p12, p21, p22,p31, p32, p41, p42, ncol=2)

arranged_g <- arrangeGrob(p11, p12, p21, p22,p31, p32, p41, p42, ncol=2) #generates g


ggsave(file = "fig_ex_curves.pdf", plot = arranged_g, width = 6, height = 10,  path = "publication_progress/Figures/")

##############################################################################################
##############################################################################################




results_simulation <-
  results_simulation %>%
  mutate(MAcc = as.character(MAcc)) %>%
  mutate(MAcc = as.numeric(ifelse(MAFC == "error", "0", MAcc))) %>%
  mutate(method = ifelse(is.na(always), algorithm,
                         paste0("all_", always))) %>%
  mutate(is_conv = as.character(is_conv))%>%
  mutate(is_conv = case_when(
    is_conv %in% c("FALSE", "0") ~0,
    TRUE ~ 1
  ))

res_conv <- results_simulation %>%
  filter(is_conv == 1, spacing_id == "equ") %>%
  mutate(method = case_when(
    method == "all_joint" ~ "Single 2pLL",
    method == "all_sep" ~ "Sep. 2pLL",
    method == "all_td2pLL" ~ "Always td2pLL",
    method == "anova" ~ "Two-step"
  ),
  noise_id = case_when(
    noise_id == "N1" ~ "Little (N1)",
    noise_id == "N2" ~ "Medium (N2)",
    noise_id == "N3" ~ "Large (N3)"
  ),
  model_id = case_when(
    model_id == "M1" ~ "Large (M2)",
    model_id == "M2" ~ "Small (M1)",
    model_id == "M3" ~ "None (M0)"
  )) %>%
  mutate(MAFC = as.numeric(as.character(MAFC)),
         noise_id = factor(noise_id, levels = c("Little (N1)", "Medium (N2)", "Large (N3)")),
         model_id = factor(model_id, levels = c("None (M0)", "Small (M1)", "Large (M2)")))
res_conv %>%
  ggplot(., aes(x = method, y = MAFC, fill = as.character(n_obs_id))) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  facet_grid(rows = vars(model_id), cols = vars(noise_id)) +
  labs(x = "Method",
       y = "Estimated EC50 to true EC50 (Absolute log2 deviation)",
       title = "Simulated background noise level",
       tag= "Simulated effect of exposure duration on EC50",
       fill = expression(n[obs])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        plot.tag=element_text(angle=-90),
        plot.tag.position= c(1.03, 0.6),
        plot.margin = margin(0.2,1,0,0.2, "cm"),
        legend.position = "bottom") +
coord_cartesian(ylim = c(0, 1))

ggsave("fig_sim_res_AMAFC.pdf", width = 6, height = 6,  path = "publication_progress/Figures/")
