### This script contains model for English data

suppressPackageStartupMessages({
  library(brms)
  library(tidyr)
  library(ggplot2)
  library(tidybayes)
  library(bayesplot)
  library(dplyr)
  library(projpred)
  library(assertthat)
})


melodic_df = read.csv('data/model_data.csv')
seed = 469852

english_df = melodic_df %>% 
  filter(society == "English")

# This dataset should contain 254 rows
assert_that(nrow(english_df) == 252, msg = "Row numbers are wrong")

#### Null model ####
fit.1 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ 1 + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/null_english")
fit.1 = add_criterion(fit.1, "loo", moment_match = TRUE,
                      file = "results/null_english")


#### Note Frequency ####
fit.2 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ frequency1:frequency2 + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/notefrequency_english")
fit.2 = add_criterion(fit.2, "loo", moment_match = TRUE,
                      file = "results/notefrequency_english")

#### Baseline (Note & Function Frequency) ####
fit.3 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/notefunctionsfrequency_english")
fit.3 = add_criterion(fit.3, "loo", moment_match = TRUE,
                      file = "results/notefunctionsfrequency_english")

#### Distance model - strong function ####
english_strong = english_df %>% 
  dplyr::filter(functional_change == "s")

fit.4.1 <-
  brm(data = english_strong, family = poisson,
      substitution_count ~ 
        semitonal_distance + 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/strongdistancefrequency_english")


#### Distance model - weak function ####
english_weak = english_df %>% 
  dplyr::filter(functional_change == "w")

fit.4.2 <-
  brm(data = english_weak, family = poisson,
      substitution_count ~ 
        semitonal_distance + 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/weakdistancefrequency_english")


#### Functional model ####
fit.4.3 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ functional_change + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functional_english")
fit.4.3 = add_criterion(fit.4.3, "loo", moment_match = TRUE,
                        file = "results/functional_english")

#### Function + Distance ####
fit.5 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ functional_change + 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functionplusdistance_english")
fit.5 = add_criterion(fit.5, "loo", moment_match = TRUE,
                      file = "results/functionplusdistance_english")

#### Function * Distance ####
fit.6 <-
  brm(data = english_df, family = poisson,
      substitution_count ~ functional_change * 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functiontimesdistance_english")
fit.6 = add_criterion(fit.6, "loo", moment_match = TRUE,
                      file = "results/functiontimesdistance_english")
