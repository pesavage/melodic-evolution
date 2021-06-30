library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(dplyr)
library(assertthat)
library(MASS)

seed = 88878

#### Data ####

melodic_df = read.csv('results/reviewer_modeldata.csv')

# Transformed variables
melodic_df$log_frequencyinteraction = log(melodic_df$count_1) * log(melodic_df$count_2)

# standardize variables
# melodic_df$std_semitonaldistance = melodic_df$semitonal_distance / max(melodic_df$semitonal_distance)
# melodic_df$std_count1 =  melodic_df$count_1 / max(melodic_df$count_1)
# melodic_df$std_count2 =  melodic_df$count_2 / max(melodic_df$count_2)
# melodic_df$std_freqinteraction = melodic_df$log_frequencyinteraction / max(melodic_df$log_frequencyinteraction)

melodic_df$std_semitonaldistance = (melodic_df$semitonal_distance - mean(melodic_df$semitonal_distance)) / sd(melodic_df$semitonal_distance)
melodic_df$std_count1 =  (log(melodic_df$count_1) - mean(log(melodic_df$count_1))) / sd(log(melodic_df$count_1))
melodic_df$std_count2 =  (log(melodic_df$count_2) - mean(log(melodic_df$count_2))) / sd(log(melodic_df$count_2))
melodic_df$std_freqinteraction = melodic_df$log_frequencyinteraction / max(melodic_df$log_frequencyinteraction)


# data tests
assert_that(nrow(melodic_df) == 137)
assert_that(sum(melodic_df$society == "Japanese") == 51)
assert_that(sum(melodic_df$society == "English") == 86)

#### Models ####
##
## Intercept Model
fit.1 = brm(data = melodic_df, family = binomial,
            mutation_count | trials(total_notes) ~ 1,
            prior(normal(0, 10), class = Intercept),
            seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
            control = list(max_treedepth = 15), sample_prior = TRUE,
            file = "results/intercept_bin")


# Societal difference

fit.2 = brm(data = melodic_df, family = binomial,
            mutation_count | trials(total_notes) ~ 1 + society,
            c(prior(normal(0, 10), class = Intercept), 
            prior(normal(0, 5), class = b)),
            seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
            control = list(max_treedepth = 15), sample_prior = TRUE,
            file = "results/society_bin")

# Semitonal distance model

fit.3 = brm(data = melodic_df, family = binomial,
            mutation_count | trials(total_notes) ~ 1 + semitonal_distance,
          c(prior(normal(0, 10), class = Intercept), 
            prior(normal(0, 5), class = b)),
          seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
          control = list(max_treedepth = 15), sample_prior = TRUE,
          file = "results/semitonaldistance_bin")

# Frequency

# Complete interaction
fit.4.1 = brm(data = melodic_df, family = binomial,
    mutation_count | trials(total_notes) ~ 1 + log(count_1) * log(count_2),
    c(prior(normal(0, 10), class = Intercept),
      prior(normal(0, 5), class = b)),
    seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
    control = list(max_treedepth = 15, adapt_delta = 0.99), 
    sample_prior = TRUE, file = "results/countfullint_bin",
    init_r = 0, save_pars = save_pars(all = TRUE))

# Interaction only
fit.4.2 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + log_frequencyinteraction,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/countintonly_bin",
              init_r = 0, save_pars = save_pars(all = TRUE))

# Minimum frequency
fit.4.2 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + minimum_count,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/mincount_bin",
              init_r = 0, save_pars = save_pars(all = TRUE))

fit.4.3 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + count_1 : count_2,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/countfullintnolog_bin",
              init_r = 0, save_pars = save_pars(all = TRUE))

#### Full Models

fit.5.1 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + log(count_1) * log(count_2) + semitonal_distance + society,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/fullcompleteint_bin",
              init_r = 0)

fit.5.2 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + log_frequencyinteraction + semitonal_distance + society,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/fullintonly_bin",
              init_r = 0)

fit.5.3 = brm(data = melodic_df, family = binomial,
              mutation_count | trials(total_notes) ~ 1 + log(count_1) * log(count_2) * semitonal_distance + society,
              c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 5), class = b)),
              seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
              control = list(max_treedepth = 15, adapt_delta = 0.99), 
              sample_prior = TRUE, file = "results/threewayint_bin",
              init_r = 0)

fit.5.1.1 = brm(data = melodic_df, family = binomial,
                mutation_count | trials(total_notes) ~ 1 + std_count1 * std_count2 + std_semitonaldistance + society,
                c(prior(normal(0, 10), class = Intercept),
                  prior(normal(0, 5), class = b)),
                seed = 10, iter = 10000, warmup = 5000, chains = 2, cores = 2,
                control = list(max_treedepth = 15, adapt_delta = 0.99), 
                sample_prior = TRUE, file = "results/fullcompleteint_stdbin",
                init_r = 0)

# Compare to a basic linear model
library(lmerTest)

fit.5.1.1 = glm(mutation_count / total_notes ~ semitonal_distance + count_1 * count_2 + society, data = melodic_df)

pred.5.1.1 = predict(fit.5.1.1)
post.5.1 = posterior_predict(fit.5.1)
postmean.5.1 = colMeans(post.5.1)

plot(y = melodic_df$mutation_count, x = pred.5.1.1 * melodic_df$total_notes,
     main = "Binomial (red) vs Gaussian (black) model", 
     ylab = "True response", xlab = "Response Estimate")
points(y = melodic_df$mutation_count, x = postmean.5.1, col = "red")
lines(y = melodic_df$mutation_count, x = melodic_df$mutation_count, col = "blue")
abline(v = 0, lty = "dashed")

mse.5.1 = 1 / length(postmean.5.1) * sum((melodic_df$mutation_count - postmean.5.1)^2)
mse.5.1.1 = 1 / length(pred.5.1.1) * sum((melodic_df$mutation_count - pred.5.1.1)^2)

bayes_R2(fit.5.1)
summary.5.1.1 = summary(fit.5.1.1)
1 - (summary.5.1.1$deviance/summary.5.1.1$null.deviance)


## 
melodic_english = melodic_df[melodic_df$society == "English",]
melodic_japanese = melodic_df[melodic_df$society == "Japanese",]

melodic_df$countint = melodic_df$count_1 * melodic_df$count_2

fit.7 = brm(data = melodic_df, family = binomial,
            mutation_count | trials(total_notes) ~ 0 + std_count1:std_count2 * std_semitonaldistance * society,
            # c(prior(normal(1, 1), class = b, coef = Intercept),
              c(prior(normal(0, 1), class = b, coef = std_count1:std_count2),
              prior(normal(0, 1), class = b, coef = std_semitonaldistance),
              prior(normal(0, 1), class = b, coef = societyJapanese)),
            seed = 10, iter = 2000, warmup = 1000, chains = 2, cores = 2,
            control = list(max_treedepth = 15), file = "results/std_int",
            sample_prior = TRUE,
            init_r = 0)

pp_check(fit.7, nsamples = 100)

fit.7.1 = glm(mutation_count ~ std_count1:std_count2 + std_semitonaldistance, family = poisson,
              data = melodic_df)

fit.7.2 = glm(mutation_count ~ count_1 : count_2 + semitonal_distance, family = poisson,
              data = melodic_japanese)

fit.7.3 = glm(mutation_count ~ count_1 : count_2 + society + semitonal_distance , 
              data = melodic_df, family = poisson)


pred.7.3 = predict(fit.7.3, type = "response")

post.7 = posterior_predict(fit.7)
postmean.7 = colMeans(post.7)


plot(y = melodic_df$mutation_count, x = postmean.7,
     main = "Binomial (red) vs Poisson (black) model", 
     ylab = "True response", xlab = "Response Estimate")


plot(y = melodic_df$mutation_count, x = pred.7.3,
     main = "Binomial (red) vs Poisson (black) model", 
     ylab = "True response", xlab = "Response Estimate")
points(y = melodic_df$mutation_count, x = postmean.5.1, col = "red")
lines(y = melodic_df$mutation_count, x = melodic_df$mutation_count, col = "blue")
abline(v = 0, lty = "dashed")

melodic_df$label = paste0(melodic_df$note1, "-", melodic_df$note2)

melodic_df$pred.7.3 = pred.7.3

ggplot(melodic_df, aes(y = mutation_count, x = count_1 * count_2,
                       label = label, size = semitonal_distance)) +
  geom_text() +
  geom_smooth(method='lm', se = FALSE) + 
  theme(legend.position = "none") +
  ylab("Number of Subsitutions") + xlab("Note 1 count x Note 2 count") +
  facet_wrap(~society, scales = "free") + 
  geom_text(aes(y = pred.7.3, x = count_1 * count_2, label = label), col = "red") + 
  ggtitle("Substitutions against count interaction in English and Japanese songs")


bayes_R2(fit.5.1)
1 - (fit.7.3$deviance/fit.7.3$null.deviance)

mse.5.1 = 1 / length(postmean.5.1) * sum((melodic_df$mutation_count - postmean.5.1)^2)
mse.7.3 = 1 / length(pred.7.3) * sum((melodic_df$mutation_count - pred.7.3)^2)
