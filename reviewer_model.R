library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(dplyr)

melodic_df = read.csv('results/reviewer_modeldata.csv')
seed = 88878

# Extra variables
melodic_df$frequency_interaction = melodic_df$frequency_1 * melodic_df$frequency_2
melodic_df$change_proportion = melodic_df$change_frequency / melodic_df$total_changes
melodic_df$log_frequencyinteraction = log(melodic_df$frequency_1) * log(melodic_df$frequency_2)

# center predictor variables
melodic_df$c_semitonaldistance = melodic_df$semitonal_distance - mean(melodic_df$semitonal_distance)
melodic_df$c_logfrequency1 = log(melodic_df$frequency_1) - mean(log(melodic_df$frequency_1))
melodic_df$c_logfrequency2 = log(melodic_df$frequency_2) - mean(log(melodic_df$frequency_2))
melodic_df$c_frequencyinteraction = melodic_df$frequency_interaction - mean(melodic_df$frequency_interaction)

# standardize variables
melodic_df$std_semitonaldistance = melodic_df$semitonal_distance / max(melodic_df$semitonal_distance)
melodic_df$std_frequency1 =  melodic_df$frequency_1 / max(melodic_df$frequency_1)
melodic_df$std_frequency2 =  melodic_df$frequency_2 / max(melodic_df$frequency_2)


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
      change_frequency | trials(total_notes) ~ 1,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/intercept_bin")

fit.1.1 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ 1,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/intercept_negbin")

fit.1.1 = add_criterion(fit.1.1, "loo")

fixef(fit.1)[1] %>% exp()

# pp_check(fit.1, nsamples = 100)
# pp_check(fit.1.1, nsamples = 500)

fit.2 <-
  brm(data = melodic_df, family = negbinomial(),
      change_frequency ~ 1 + society,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15, adapt_delta = 0.99), 
      sample_prior = TRUE, file = "results/society_negbinom")

fit.2.1 <-
  brm(data = melodic_df, family = binomial,
      change_frequency | trials(total_notes) ~ 1 + society,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/society_bin")

fixef(fit.2.1) %>%
  exp()

# pp_check(fit.2.1, nsamples = 500)

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

fit.3.1 <-
  brm(data = melodic_df, family = binomial,
      change_frequency | trials(total_notes) ~ 1 + semitonal_distance + society,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/semitonaldistance_bin")

summary(fit.3.1)

fixef(fit.3.1) %>%
  exp()

pp_check(fit.3.1, nsamples = 500)

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


#### Mutability ####

fit.4.6 <-
  brm(data = melodic_df, family = binomial,
      change_frequency | trials(total_notes) ~ 1 + mutability_1 * mutability_2 + society,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/frequency_bin")

summary(fit.4.6)
pp_check(fit.4.6, nsamples = 100)
conditional_effects(fit.4.6, effects = "mutability_1:mutability_2",
                    int_conditions = list(mutability_2 = seq(0, 1, by = 0.2)))

conditional_effects(fit.4.6, conditions = 
                      data.frame(total_notes = unique(melodic_df$total_notes)))


fit.4.7 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ mutability_1 * mutability_2 + (1|society),
              c(prior(cauchy(0, 4), class = Intercept),
                prior(cauchy(0, 2), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99),
              sample_prior = TRUE, file = "results/mutabilityint_negbin")

fit.4.6 = add_criterion(fit.4.6, "loo", moment_match = TRUE)

conditional_effects(fit.4.6, effects = "mutability_1:mutability_2",
                    int_conditions = list(mutability_2 = seq(0, 1, by = 0.2)))


pp_check(fit.4.6, nsamples = 100)
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

fit.5.2 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ std_semitonaldistance + std_frequency1 * std_frequency2  + (1|society),
              c(prior(cauchy(0, 2), class = Intercept),
                prior(cauchy(0, 2), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 20000, warmup = 15000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/fullstd_negbin")

pp_check(fit.5.2, nsamples = 500)

fit.5.3 = brm(data = melodic_df, family = negbinomial(),
              change_frequency ~ semitonal_distance + mutability_1 * mutability_2  + (1|society),
              c(prior(cauchy(0, 2), class = Intercept),
                prior(cauchy(0, 2), class = b),
                prior(cauchy(0, 0.5), class = sd)),
              seed = 10, iter = 20000, warmup = 15000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/fullmutability_negbin")

fit.5.3 = add_criterion(fit.5.3, "loo", moment_match = TRUE)


pp_check(fit.5.3, nsamples = 100)

loo_compare(fit.5.3, fit.4.6)


fit.6 <-
  brm(data = melodic_df, family = binomial,
      change_frequency | trials(total_notes) ~ 1 + semitonal_distance + mutability_1 * mutability_2,
      prior(normal(0, 10), class = Intercept),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/full_bin")

summary(fit.6)

pp_check(fit.6, nsamples = 100)

melodic_df$label = paste(melodic_df$note1, "-", melodic_df$note2)

conditional_effects(fit.6, conditions = 
                      data.frame(total_notes = unique(melodic_df$total_notes)))

ggplot(melodic_df, aes(y = change_frequency / total_changes, 
                       x = mutability_1 + mutability_2 + mutability_1 * mutability_2, 
                       col = society,
                       size = semitonal_distance)) +
  geom_text(label = melodic_df$label)

ggplot(melodic_df, aes(y = frequency_1, x = mutability_1, col = society)) +
  geom_jitter()
