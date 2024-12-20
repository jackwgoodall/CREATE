#### POWER CALCULATIONS FOR CREATE FELLOWSHIP APPLICATION #### 

# This detailed the power calculations for Staph aureus and Group A Streptococcus
# The number of viral events is already known from our dataset (18% of all swabs were taken when the participant had a viral infection in the last 14 days)
# The intercept is calculated from this - log(baseline_prob / (1 - baseline_prob)
# The random effect stuctures are also pre-determined and their variance has been extracted from the basic model of Strep pneumoniae

library(simr) # Critical that this is loaded before the tidyverse
library(lme4)
# lastResult()$errors finds errors 

set.seed(12345)

sim_data <- STPN %>%
  select(pid, household, any_within_14_days) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  mutate(STPN = rbinom(nrow(.), 1, 0.2),
         GAS_low = rbinom(nrow(.), 1, 0.046),
         GAS = rbinom(nrow(.), 1, 0.08),
         STAU = rbinom(nrow(.), 1, 0.264))

#save(sim_data, file = "data/objects/sim_data.Rdata")

# load synthetic data from the main project
# load("~/Documents/R projects/Frequentist_STPN/data/objects/sim_data.Rdata")

sim_data$STPN <- factor(sim_data$STPN)
sim_data$any_within_14_days <- factor(sim_data$any_within_14_days)
sim_data$STAU <- factor(sim_data$STAU)
sim_data$GAS <- factor(sim_data$GAS)
sim_data$GAS_low <- factor(sim_data$GAS_low)

## STPN as a test
# Calculate the intercept
STPN_p <- sum(sim_data$STPN == 1) / sum(!is.na(sim_data$STPN)) # probability of STPN positive at baseline

STPN_intercept <- log(STPN_p/(1-STPN_p))

# Specify random effects variance
random_effects <- list("pid:household" = 0.1331, "household" = 0.2288)

# Create the model (STPN - with known output just to check this is all working!)
STPN_model <- makeGlmer(
  STPN ~ any_within_14_days + (1 | household/pid),
  fixef = c(STPN_intercept, log(2.25)),
  VarCorr = random_effects,
  family = "binomial",
  data = sim_data
)

# Run power analysis with effect size of 0 (null hypothesis)
STPN_power_result <- powerSim(STPN_model, fixed("any_within_14_daysTRUE"), nsim = 3)

print(STPN_power_result)

### STAU
# Calculate the intercept
STAU_p <- sum(sim_data$STAU == 1) / sum(!is.na(sim_data$STAU)) # probability of STPN positive at baseline

STAU_intercept <- log(STAU_p/(1-STAU_p))

# Create the model (STPN - with known output just to check this is all working!)
STAU_model <- makeGlmer(
  STAU ~ any_within_14_days + (1 | household/pid),
  fixef = c(STAU_intercept, log(2.25)),
  VarCorr = random_effects,
  family = "binomial",
  data = sim_data
)

# Run power analysis with effect size of 0 (null hypothesis)
STAU_power_result <- powerSim(STAU_model, fixed("any_within_14_daysTRUE"), nsim = 3)

print(STAU_power_result)

### GAS 
GAS_p <- sum(sim_data$GAS == 1) / sum(!is.na(sim_data$GAS)) # probability of STPN positive at baseline

GAS_intercept <- log(GAS_p/(1-GAS_p))

# Create the model (STPN - with known output just to check this is all working!)
GAS_model <- makeGlmer(
  GAS ~ any_within_14_days + (1 | household/pid),
  fixef = c(GAS_intercept, log(2.25)),
  VarCorr = random_effects,
  family = "binomial",
  data = sim_data
)

# Run power analysis with effect size of 0 (null hypothesis)
GAS_power_result <- powerSim(GAS_model, fixed("any_within_14_daysTRUE"), nsim = 3)

print(GAS_power_result)

### GAS_low 
GAS_low_p <- sum(sim_data$GAS_low == 1) / sum(!is.na(sim_data$GAS_low)) 
GAS_low_intercept <- log(GAS_low_p/(1-GAS_low_p))

# Create the model (STPN - with known output just to check this is all working!)
GAS_low_model <- makeGlmer(
  GAS_low ~ any_within_14_days + (1 | household/pid),
  fixef = c(GAS_low_intercept, log(2.25)),
  VarCorr = random_effects,
  family = "binomial",
  data = sim_data
)

# Run power analysis with effect size of 0 (null hypothesis)
GAS_low_power_result <- powerSim(GAS_low_model, fixed("any_within_14_daysTRUE"), nsim = 3)

print(GAS_low_power_result)



###### OPTIMAL EFFECT SIZE CALCULATIONS AND GRAPHS #######

# optimal power
effect_sizes <- seq(1.1, 2, by = 0.1) # Odds ratios from 1 (null) to 3

STAU_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())


for (effect in effect_sizes) {
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    STAU ~ any_within_14_days + (1 | household/pid),
    fixef = c(STAU_intercept, log_effect),
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


GAS_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())


for (effect in effect_sizes) {
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    GAS ~ any_within_14_days + (1 | household/pid),
    fixef = c(GAS_intercept, log_effect),
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

## STAU low - optimiser
STAU_low_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())


for (effect in effect_sizes) {
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    STAU_low ~ any_within_14_days + (1 | household/pid),
    fixef = c(STAU_low_intercept, log_effect),
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
  
  STAU_low_power_results <- rbind(
    STAU_low_power_results,
    data.frame(effect_size = effect, 
               power_mean = power_mean,
               power_upper = power_upper,
               power_lower = power_lower
    ))
}

## GAS low
GAS_low_power_results <- data.frame(
  effect_size = numeric(), 
  power_mean = numeric(),
  power_upper = numeric(),
  power_lower = numeric())


for (effect in effect_sizes) {
  # Set the fixed effect for the current odds ratio
  log_effect <- log(effect)
  current_model <- makeGlmer(
    GAS_low ~ any_within_14_days + (1 | household/pid),
    fixef = c(GAS_low_intercept, log_effect),
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
  geom_line() +
  geom_errorbar(aes(ymin = power_lower, ymax = power_upper), width = 0, alpha = 0.5, size = 0.8) + 
  geom_hline(yintercept = 80, linetype = 2) + 
  labs(x = "Effect size",
       y = "Predicted power",
       title = "Power estimations",
       colour = "Organism") + 
  scale_colour_manual(values = c(
    "Group A Strep expected (8%)" = "#0099FF", 
    "Staph aureus (26.4%)" = "#ff7f0e", 
    "Group A Strep low (4.6%)" = "#000066"  
  ))


