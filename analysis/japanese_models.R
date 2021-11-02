### This script contains model for Japanese data

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

# add note pairs
melodic_df$notepair = paste0(melodic_df$note1, melodic_df$note2)

japanese_df = melodic_df %>% 
  filter(society == "Japanese")

# This dataset should contain 254 rows
assert_that(nrow(japanese_df) == 252, msg = "Row numbers are wrong")

#### Null model ####
fit.1 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ 1 + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/null_japanese")
fit.1 = add_criterion(fit.1, "loo", moment_match = TRUE, 
                      file = "results/null_japanese")

#### Note Frequency ####
fit.2 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ frequency1:frequency2 + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/notefrequency_japanese")
fit.2 = add_criterion(fit.2, "loo", moment_match = TRUE,
                      file = "results/notefrequency_japanese")


#### Baseline (Note & Function Frequency) ####
fit.3 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/notefunctionsfrequency_japanese")
fit.3 = add_criterion(fit.3, "loo", moment_match = TRUE,
                      file = "results/notefunctionsfrequency_japanese")


#### Distance model - strong function ####
japanese_strong = japanese_df %>% 
  dplyr::filter(functional_change == "s")

fit.4.1 <-
  brm(data = japanese_strong, family = poisson,
      substitution_count ~ 
        semitonal_distance + 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/strongdistancefrequency_japanese")

#### Distance model - weak function ####
japanese_weak = japanese_df %>% 
  dplyr::filter(functional_change == "w")

fit.4.2 <-
  brm(data = japanese_weak, family = poisson,
      substitution_count ~ 
        semitonal_distance + 
        frequency1:frequency2:functional_total + 
        (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/weakdistancefrequency_japanese")

#### Functional model ####
fit.4.3 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ functional_change + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functional_japanese")
fit.4.3 = add_criterion(fit.4.3, "loo", moment_match = TRUE,
                        file = "results/functional_japanese")

#### Function + Distance ####
fit.5 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ functional_change + 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functionplusdistance_japanese")
fit.5 = add_criterion(fit.5, "loo", moment_match = TRUE,
                      file = "results/functionplusdistance_japanese")

# bayes_R2(fit.5)
# exp(fixef(fit.5))


#### Function * Distance ####
fit.6 <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ functional_change * 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/functiontimesdistance_japanese")
fit.6 = add_criterion(fit.6, "loo", moment_match = TRUE,
                      file = "results/functiontimesdistance_japanese")
