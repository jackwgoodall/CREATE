#### POWER CALCULATIONS FOR CREATE FELLOWSHIP APPLICATION #### 

# This detailed the power calculations for Staph aureus and Group A Streptococcus

# The number of viral events is already known from our dataset (18% of all swabs were taken when the participant had a viral infection in the last 14 days)

# The intercept is calculated from this - log(baseline_prob / (1 - baseline_prob)

# This is then corrected to the baseline rate (i.e. the value it would be then there was no
# viral infection, assuming the effect size we predict is correct)

# The random effect stuctures are also pre-determined and their variance has been extracted from the basic model of Strep pneumoniae

library(simr) # Critical that this is loaded before the tidyverse
library(lme4)
# lastResult()$errors finds errors 

# Create function to estimate baseline prevalence
estimate_baseline_prevalence <- function(p_obs, q, OR, tol = 1e-7) {
  fn <- function(p0) {
    p1 <- (OR * p0) / (1 - p0 + OR * p0)
    p_hat <- p0 * (1 - q) + p1 * q
    return(p_hat - p_obs)
  }
  
  # Use uniroot to solve for p0
  result <- uniroot(fn, interval = c(1e-6, 0.999), tol = tol)
  return(result$root)
}

set.seed(12345)

sim_data <- STPN %>%
  select(pid, household, any_within_14_days) %>%
  filter(if_all(everything(), ~ !is.na(.))) 

# Set all prevelances
STPN_prev = 0.2
GAS_low_prev = 0.046
GAS_prev = 0.08
STAU_prev = 0.264

#save(sim_data, file = "data/objects/sim_data.Rdata")

## STPN as a test
# Calculate STPN baseline
STPN_BL <- estimate_baseline_prevalence(p_obs = 0.046, q = 0.18, OR = 2.25)

STPN_BL <- 0.046

# Specify random effects variance
random_effects <- list("pid:household" = 0.1331, "household" = 0.2288)

# Create the model (STPN - with known output just to check this is all working!)
STPN_model <- makeGlmer(
  STPN ~ any_within_14_days + (1 | household/pid),
  fixef = c(log(STPN_BL / (1 - STPN_BL)), log(2.25)),
  VarCorr = random_effects,
  family = "binomial",
  data = sim_data
)

# Run power analysis with effect size of 0 (null hypothesis)
STPN_power_result <- powerSim(STPN_model, fixed("any_within_14_daysTRUE"), nsim = 10)

print(STPN_power_result)


###### OPTIMAL EFFECT SIZE CALCULATIONS AND GRAPHS ------

effect_sizes <- seq(1.1, 2, by = 0.1) # Odds ratios from 1 (null) to 3

### STAU -----

STAU_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())

for (effect in effect_sizes) {
  
  # calculate the baseline
  STAU_BL <- estimate_baseline_prevalence(p_obs = STAU_prev, q = 0.18, OR = effect)
  
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    STAU ~ any_within_14_days + (1 | household/pid),
    fixef = c(log(STAU_BL / (1 - STAU_BL)), log(effect)),
    VarCorr = random_effects,
    family = "binomial",
    data = sim_data
  )
  
  # Run power simulation
  power_sim <- powerSim(current_model, fixed("any_within_14_daysTRUE"), nsim = 100) 
  
  # Extract power from summary
  power_summary <- summary(power_sim)
  power_mean <- as.numeric(power_summary$mean) * 100 # Convert percentage to proportion
  power_lower <- as.numeric(power_summary$lower) * 100 
  power_upper <- as.numeric(power_summary$upper) * 100 
  
  STAU_power_results <- rbind(
    STAU_power_results,
    data.frame(effect_size = effect, 
               power_mean = power_mean,
               power_upper = power_upper,
               power_lower = power_lower
    ))
}

## GAS ----

GAS_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())

for (effect in effect_sizes) {
  
  # calculate the baseline
  GAS_BL <- estimate_baseline_prevalence(p_obs = GAS_prev, q = 0.18, OR = effect)
  
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    GAS ~ any_within_14_days + (1 | household/pid),
    fixef = c(log(GAS_BL / (1 - GAS_BL)), log(effect)),
    VarCorr = random_effects,
    family = "binomial",
    data = sim_data
  )
  
  # Run power simulation
  power_sim <- powerSim(current_model, fixed("any_within_14_daysTRUE"), nsim = 100) 
  
  # Extract power from summary
  power_summary <- summary(power_sim)
  power_mean <- as.numeric(power_summary$mean) * 100 # Convert percentage to proportion
  power_lower <- as.numeric(power_summary$lower) * 100 
  power_upper <- as.numeric(power_summary$upper) * 100 
  
  GAS_power_results <- rbind(
    GAS_power_results,
    data.frame(effect_size = effect, 
               power_mean = power_mean,
               power_upper = power_upper,
               power_lower = power_lower
    ))
}


## GAS low ----
GAS_low_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())


for (effect in effect_sizes) {
  # calculate the baseline
  GAS_low_BL <- estimate_baseline_prevalence(p_obs = GAS_low_prev, q = 0.18, OR = effect)
  
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    GAS_low ~ any_within_14_days + (1 | household/pid),
    fixef = c(log(GAS_low_BL / (1 - GAS_low_BL)), log(effect)),
    VarCorr = random_effects,
    family = "binomial",
    data = sim_data
  )
  
  # Run power simulation
  power_sim <- powerSim(current_model, fixed("any_within_14_daysTRUE"), nsim = 100) 
  
  # Extract power from summary
  power_summary <- summary(power_sim)
  power_mean <- as.numeric(power_summary$mean) * 100 # Convert percentage to proportion
  power_lower <- as.numeric(power_summary$lower) * 100 
  power_upper <- as.numeric(power_summary$upper) * 100 
  
  GAS_low_power_results <- rbind(
    GAS_low_power_results,
    data.frame(effect_size = effect, 
               power_mean = power_mean,
               power_upper = power_upper,
               power_lower = power_lower
    ))
}


#### Combine
combined_power <- rbind(cbind(
  GAS_low_power_results, 
  group = rep("Group A Strep low (4.6%)", nrow(GAS_low_power_results))
),
cbind(
  GAS_power_results, 
  group = rep("Group A Strep expected (8%)", nrow(GAS_power_results))
),
cbind(
  STAU_power_results, 
  group = rep("Staph aureus (26.4%)", nrow(STAU_power_results))
))


## Plot 
library(ggplot2)
ggplot(combined_power, aes(x = effect_size, y = power_mean, colour = group)) + 
  geom_point() + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = power_lower, ymax = power_upper), width = 0, alpha = 0.5, size = 1.5) + 
  geom_hline(yintercept = 80, linetype = 2) + 
  labs(x = "Effect size",
       y = "Predicted power",
       title = "Power estimations",
       colour = "Organism") + 
  scale_colour_manual(values = c(
    "Group A Strep expected (8%)" = "#0099FF", 
    "Staph aureus (26.4%)" = "#ff7f0e", 
    "Group A Strep low (4.6%)" = "#000066"  
  )) + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 14))






