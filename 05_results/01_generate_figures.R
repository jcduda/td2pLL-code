
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggplot2); theme_set(theme_bw())
library(td2pLL)
library(DoseFinding)
library(gridExtra)

setwd("C:/Users/duda/Projekte/exp_dose_resp_code")

source("./00_functions/functions.R")
load("./02_data_prep/data_refit_prep.RData")
load("./03_simulation_prep//01_true_models.RData")
load("./03_simulation_prep/02_true_designs.RData")
load("./03_simulation_prep/03_plot_sd_values.RData")
load("./03_simulation_prep/03_lm_sd_noise.RData")
load("./04_simulation_run/results_simulation.RData")

data(cytotox) # loaded through td2pLL package

# TODO:
# Generate the figure that was proposed by the reviewer and re-name the figure numbers.



############
# Figure 1 #
############

# Note: The figures appear in the Viewer after a second and were manually exported and saved.
# Within the viewer, you can move the plot interactively.

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



#############
## Figure 3
#############

# Note: For the manuscript, the models where renames as follows:
# 2a: No time-dependency:       M3 (Code) --> M0 (Paper)
# 2b: Weak time-dependency:     M2 (Code) --> M1 (Paper)
# 2c: Strong time-dependency:   M1 (Code) --> M2 (Paper)
#
# The graphics appear in the Viewer and were exported and cropped manually.

# 3a (No time-dependency)
plot_td2pLL(M3, dose_lim = c(1e-3, 1), time_lim = c(1, 7),
            xaxis_title = "Concentration", yaxis_title = "Exposure duration")

# 3b (Weak time-dependency)
plot(M2, xaxis_scale = "log", xaxis_title = "Concentration",
     yaxis_title = "Exposure duration", zaxis_title = "Response")

# 3c (Strong time-dependency)
plot(M1, xaxis_scale = "log", xaxis_title = "Concentration",
     yaxis_title = "Exposure duration", zaxis_title = "Response")


#############
## Figure 4
#############

my_plot_designs(designs = designs %>% filter(n_obs == 216,
                                             spacing == "equ") %>%
                  mutate(pt_size = ifelse(n == 12, 5, 4),
                         n_times = paste0(n_times," exposure durations")))

ggsave("./05_results/fig_4_expDes_sim.pdf", width = 5, height = 2.75)


#############
## Figure 5
#############

plot_sd_values +
  labs(x = expression(Concentration~x[j]))

ggsave("./05_results/fig_5_sd_sim.pdf", width = 7.5, height = 2.75)


##################
## Figure 7 and 6
##################

# Figure 7

all_compounds <- cytotox$compound %>% unique
app_res <- data.frame(compound = all_compounds,
                      td2pLL_rse = NA,
                      sep_rse = NA,
                      joint_rse = NA,
                      anova_signif = NA)

model_list <- myList <- vector("list", length(all_compounds))
names(model_list) <- all_compounds

pdf("./05_results/fig_7_all_compounds.pdf", width = 10, height = 12)
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
  
  # separate
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


# Figure 6

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

ggsave("./05_results/fig_6_td2pLL_app.pdf", width = 4, height = 4)


############
## Figure 8
############

# Selections:
# PPL, GLC: Separate better
# LAB: td2pLL looks more stable
# FAM: No computational problems in td2pLL
# DFN: td2pLL more stable as ED50 should decrease with time 



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

ggsave(file = "./05_results/fig_8_ex_curves.pdf", plot = arranged_g, width = 6, height = 10)


############
## Figure 9
############

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
  # For the manuscript, consider only equally spaced scenarios.
  # Scenarios with optimal designs were not scope of the manuscript.
  #
  # Also, "is_conv ==1" means that only cases where the model fit converged are
  # considered. 
  # This is negligible as only 27 of 144000 rows have is_conv == 0, i.e. non.convergence.
  filter(is_conv == 1, spacing_id == "equ") %>%
  mutate(method = case_when(
    # method names were slightly renamed for the manuscript!
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
    # Model names were renamed for the manuscript!
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

ggsave("./05_results/fig_9_sim_res_AMAFC.pdf", width = 6, height = 6)
