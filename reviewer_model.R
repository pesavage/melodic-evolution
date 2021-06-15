library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(dplyr)

melodic_df = read.csv('results/reviewer_modeldata.csv')
seed = 88877

# Extra variables
melodic_df$frequency_interaction = melodic_df$frequency_1 * melodic_df$frequency_2
melodic_df$change_proportion = melodic_df$change_frequency / melodic_df$total_changes
melodic_df$log_frequencyinteraction = log(melodic_df$frequency_1) * log(melodic_df$frequency_2)

# center predictor variables
melodic_df$c_semitonaldistance = melodic_df$semitonal_distance - mean(melodic_df$semitonal_distance)
melodic_df$c_frequency1 = melodic_df$frequency_1 - mean(melodic_df$frequency_1)
melodic_df$c_frequency2 = melodic_df$frequency_2 - mean(melodic_df$frequency_2)
melodic_df$c_frequencyinteraction = melodic_df$frequency_interaction - mean(melodic_df$frequency_interaction)

#### Summary ####

melodic_df %>% 
  group_by(society) %>% 
  summarise(mean(semitonal_distance), 
            mean(frequency_1), 
            mean(frequency_2), 
            mean(frequency_interaction))

tapply(melodic_df$change_proportion, melodic_df$society, mean)
tapply(melodic_df$change_proportion, melodic_df$society, median)
tapply(melodic_df$change_proportion, melodic_df$society, sd)

#### Plots of variables ####
melodic_long = pivot_longer(melodic_df,
                            cols = c("change_frequency", "semitonal_distance",
                                     "frequency_1", "frequency_2", "minimum_frequency"))

ggplot(melodic_long, aes(x = value, group = society)) +
  geom_density() +
  facet_wrap(~name, scales = "free")


#### Models ####

fit.1 <-
  brm(data = melodic_df, family = binomial,
      change_frequency | trials(total_changes) ~ 1,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/intercept_model")

fit.1.1 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ 1,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/intercept_negbin")


fixef(fit.1.1)[1] %>% exp()

# pp_check(fit.1, nsamples = 500)
# pp_check(fit.1.1, nsamples = 500)

fit.2 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ 1 + society,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE, file = "results/society_negbinom")

fixef(fit.2) %>%
  exp()

# pp_check(fit.2, nsamples = 500)

plot(conditional_effects(fit.2), ask = FALSE)

melodic_df %>% 
  mutate(change_prop = change_frequency / total_changes) %>% 
  group_by(society) %>% 
  summarise(mean(change_frequency))

fit.3 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ semitonal_distance + (1|society),
      c(prior(normal(0, 4), class = Intercept),
        prior(cauchy(0, 2), class = b),
        prior(cauchy(0, 0.5), class = sd)),
      seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE,
      file = "results/semitonaldistance_negbin")


summary(fit.3)

fixef(fit.3) %>%
  exp()

pp_check(fit.3, nsamples = 500)

#### Frequency ####

fit.4.1 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ log(frequency_1) + (1|society),
      c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 10), class = b),
        prior(cauchy(0, 5), class = sd)),
      seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE, file = "results/frequency1_negbin",
      init_r = 0)

fit.4.2 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ log(frequency_2) + (1|society),
      c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 10), class = b),
        prior(cauchy(0, 5), class = sd)),
      seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE, file = "results/frequency2_negbin",
      init_r = 0)

summary(fit.4.2)

fit.4.3 = brm(data = melodic_df, family = negbinomial(),
                change_frequency ~ log(frequency_1) * log(frequency_2) + (1|society),
                c(prior(normal(0, 4), class = Intercept),
                  prior(cauchy(0, 2), class = b),
                  prior(cauchy(0, 0.5), class = sd)),
                seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
                control = list(max_treedepth = 15, adapt_delta = 0.99), 
                sample_prior = TRUE, file = "results/frequencyint_negbin",
                init_r = 0, save_pars = save_pars(all = TRUE))

fit.4.1 = add_criterion(fit.4.1, "loo")
fit.4.2 = add_criterion(fit.4.2, "loo")
fit.4.3 = add_criterion(fit.4.3, "loo", moment_match = TRUE)

loo_compare(fit.4.1, fit.4.2, fit.4.3)

pp_check(fit.4.1, nsamples = 500)
pp_check(fit.4.2, nsamples = 500)
pp_check(fit.4.3, nsamples = 500)

fit.4.4 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ log_frequencyinteraction + (1|society),
              c(prior(cauchy(0, 0.5), class = Intercept),
                prior(cauchy(0, 0.5), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 20000, warmup = 15000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/frequencyint2_negbin")

fit.4.4 = add_criterion(fit.4.4, "loo")
loo_compare(fit.4.1, fit.4.2, fit.4.3, fit.4.4)


fit.4.5 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ minimum_frequency + (1|society),
              c(prior(cauchy(0, 4), class = Intercept),
                prior(cauchy(0, 2), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/minfrequency_negbin")


fit.4.5 = add_criterion(fit.4.5, "loo")
loo_compare(fit.4.1, fit.4.2, fit.4.3, fit.4.4, fit.4.5)

#### Full Model ####

fit.5 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ semitonal_distance + log(frequency_1) * log(frequency_2) + (1|society),
              c(prior(normal(0, 4), class = Intercept),
                prior(cauchy(0, 2), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 20000, warmup = 15000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/full_negbin")

fit.5.1 = brm(data = melodic_df, family = negbinomial(),
            change_frequency ~ semitonal_distance + log_frequencyinteraction + (1|society),
            c(prior(cauchy(0, 0.5), class = Intercept),
              prior(cauchy(0, 0.5), class = b),
              prior(cauchy(0, 0.5), class = sd)),
            seed = 10, iter = 20000, warmup = 15000, chains = 2, cores = 2,
            control = list(max_treedepth = 15, adapt_delta = 0.99), 
            sample_prior = TRUE, file = "results/full2_negbin")

summary(fit.5)
summary(fit.5.1)


pp_check(fit.5, nsamples = 500)
