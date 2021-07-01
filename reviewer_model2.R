library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(dplyr)

melodic_df = read.csv('results/reviewer_modeldata.csv')
seed = 88878

# standardize predictor variables
# melodic_df$std_semitonaldistance = melodic_df$semitonal_distance - mean(melodic_df$semitonal_distance) / sd(melodic_df$semitonal_distance)
# melodic_df$std_count1 = melodic_df$count_1 - mean(melodic_df$count_1) / sd(melodic_df$count_1)
# melodic_df$std_count2 = melodic_df$count_2 - mean(melodic_df$count_2) / sd(melodic_df$count_2)

melodic_df$std_semitonaldistance = melodic_df$semitonal_distance / max(melodic_df$semitonal_distance)
melodic_df$std_count1 = melodic_df$count_1 / max(melodic_df$count_1)
melodic_df$std_count2 = melodic_df$count_2 / max(melodic_df$count_2)

# save standardized variables
write.csv(melodic_df, "results/reviewerstd_modeldata.csv")

## null model
fit.1 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ std_count1:std_count2,
      prior(normal(0, 1), coef = "std_count1:std_count2"),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/null_bin")

pred.1 = predict(fit.1)
residuals.1 = residuals(fit.1)
plot(y = melodic_df$mutation_count, x = pred.1[,1])
plot(melodic_df$mutation_count - pred.1[,1])
abline(h = -0)

fit.1.3 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2:society,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/societyrslopes_bin")

pp_check(fit.1.1)

fit.1.2 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2 + society,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = 'only', control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/nullrint_bin")

pp_check(fit.1.2, nsamples = 50) + xlim(c(0, 100))

fit.1 = add_criterion(fit.1, "loo", moment_match = TRUE)
fit.1.1 = add_criterion(fit.1.1, "loo", moment_match = TRUE)
fit.1.2 = add_criterion(fit.1.2, "loo", moment_match = TRUE)

loo_compare(fit.1, fit.1.1)

# semitones

fit.2 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2:society + std_semitonaldistance,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/semitones_bin")

pp_check(fit.2)

summary(fit.2)

# functional
fit.3 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2:society + functional_change,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/functional_bin")

pp_check(fit.3)
summary(fit.3)

# full model
fit.4 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2:society + std_semitonaldistance + functional_change,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/full_bin")

fit.4.1 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
         std_count1:std_count2:society:functional_change + std_semitonaldistance,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/fullrslopes_bin")

fit.4.2 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(total_notes) ~ 
        std_count1:std_count2:society + functional_change + std_semitonaldistance,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      file = "results/fullrslopes_bin2")

predict(fit.4.2, newdata = data.frame(std_count1 = 100000, std_count2 = 0, 
                                      total_notes = 11717, 
                                      functional_change = "NF-NF",
                                      std_semitonaldistance = 0.25,
                                      society = "English"))

pp_check(fit.4, nsamples = 100)

resid.4 = residuals(fit.4.2)
pred.4 = predict(fit.4.2)
plot(y = melodic_df$mutation_count, pred.4[,1])
points(y = melodic_df$mutation_count, x = melodic_df$mutation_count, col = "red")

# add model criterion
fit.1 = add_criterion(fit.1, "loo", moment_match = TRUE)
fit.1.3 = add_criterion(fit.1.3, "loo", moment_match = TRUE)
fit.2 = add_criterion(fit.2, "loo", moment_match = TRUE)
fit.3 = add_criterion(fit.3, "loo", moment_match = TRUE)
fit.4 = add_criterion(fit.4, "loo", moment_match = TRUE)

loo_compare(fit.1, fit.2, fit.3, fit.4)


post.4 = posterior_samples(fit.4)
