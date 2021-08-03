library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(bayesplot)
library(dplyr)
library(projpred)

# melodic_df = read.csv('results/reviewer_modeldata.csv')
melodic_df = read.csv("results/reviewer_wonf_f_modeldata.csv")
seed = 88878

melodic_df$std_semitonaldistance = melodic_df$semitonal_distance / max(melodic_df$semitonal_distance)
melodic_df$std_count1 = melodic_df$count_1 / max(melodic_df$count_1)
melodic_df$std_count2 = melodic_df$count_2 / max(melodic_df$count_2)
melodic_df$society = as.factor(melodic_df$society)

# change levels
melodic_df$functional_change = factor(melodic_df$functional_change, levels = c("NF-NF", "F-NF", "F-F"))

# Null model
fit.1 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ std_count1:std_count2,
      prior(normal(0, 1), coef = "std_count1:std_count2"),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/null_bin_nonff")

#### Bi-Variate #### 
# "Y~society"            
fit.2.1 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2:society,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/societyrslopes_bin_nonff")

# "Y~semitonal_distance" 
fit.2.2 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2 + std_semitonaldistance,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = TRUE, control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/semitonal_bin_nonff")

# "Y~functional"  
fit.2.3 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2 + functionalchange_bin,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = TRUE, control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/functional_bin_nonff")

#### Tr-Variate ####
# "Y~society+semitonal_distance"  
fit.3.1 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2:society + std_semitonaldistance,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = TRUE, control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/societysemitonal_bin_nonff")


# "Y~society+functional"       
fit.3.2 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2:society + functionalchange_bin,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = TRUE, control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/societyfunctional_bin_nonff")

# "Y~semitonal_distance+functional"
fit.3.3 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2 + std_semitonaldistance + functionalchange_bin,
      prior(normal(0, 1), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      sample_prior = TRUE, control = list(max_treedepth = 15),
      save_pars = save_pars(all = TRUE),
      file = "results/semitonalfunctional_bin_nonff")

#### Full model ####
# "Y~society+semitonal_distance+functional"
fit.4 <-
  brm(data = melodic_df, family = binomial,
      mutation_count | trials(functional_total) ~ 
        std_count1:std_count2:society + functionalchange_bin + std_semitonaldistance,
      prior(normal(0, 5), class = b),
      seed = 10, iter = 6000, warmup = 3000, chains = 2, cores = 2,
      control = list(max_treedepth = 15), sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      file = "results/fullrslopes_bin_nonff")

summary(fit.4)


#### Model comparison ####
loo_comparison  = loo(fit.1, 
                        fit.2.1, fit.2.2, fit.2.3, 
                        fit.3.1, fit.3.2, fit.3.3,
                        fit.4, moment_match = TRUE)

plot_loo = data.frame(loo_comparison$diffs)
plot_loo$models = c("Society + Semitonal distance + Function",
                    "Semitonal distance + Function", 
                    "Society + Semitonal distance",
                    "Semitonal distance", 
                    "Society + Function",
                    "Function",
                    "Society",
                    "Null")

plot_loo$models = factor(plot_loo$models, 
                         levels = plot_loo$models[order(plot_loo$elpd_diff)])

plot_1 = ggplot(plot_loo, aes(x = elpd_diff, y = models)) + 
  geom_point() + 
  geom_errorbar(aes(xmin=elpd_diff-se_diff, xmax=elpd_diff+se_diff), 
                width=.2, position=position_dodge(.9)) + 
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("") + xlab("ELPD Difference") + 
  theme_light() + 
  theme(text = element_text(size=16, family="serif"))

ggsave(plot_1, filename = "figures/loo_compare_nonff.png")

mcmc_areas(as.matrix(fit.4)[,2:6])

bayes_R2(fit.4)
bayes_R2(fit.4) - bayes_R2(fit.3.1) # society + semitonal distance
bayes_R2(fit.4) - bayes_R2(fit.3.2) # society + Function

## 

post.4 = posterior_samples(fit.4)

quantile(post.4$b_std_semitonaldistance, c(.5, .025, .75))
quantile(inv_logit_scaled(post.4$b_std_semitonaldistance), c(.5, .025, .75))

newdata = data.frame(std_count1 = median(melodic_df$std_count1),
                     std_count2 = median(melodic_df$std_count2),
                     functionalchange_bin = 0,
                     std_semitonaldistance = c(0.0625, 0.1250, 0.1875, 0.2500, 0.3125, 0.3750, 0.4375),
                     society = "English", 
                     functional_total = 16240 + 16240)

pred.4 = predict(fit.4, newdata = newdata) %>% data.frame
pred.4$semitones = 1:7

plot_2 = ggplot(pred.4, aes(x = semitones, y = Estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                width=.2, position=position_dodge(.9))+
  xlab("Semitonal distance") + 
  ggtitle("Predicted substitutions with increases in semitonal distance") + 
  theme_light() + 
  theme(text = element_text(size=20, family="serif"))

ggsave(plot_2, filename = "figures/semitonal_distance_nonff.png")

# soceity 
post.4.1 = posterior_samples(fit.4)

diff_a =  post.4.1$`b_std_count1:std_count2:societyJapanese` -  post.4.1$`b_std_count1:std_count2:societyEnglish`

diff_b = inv_logit_scaled(post.4.1$`b_std_count1:std_count2:societyJapanese`) - inv_logit_scaled(post.4.1$`b_std_count1:std_count2:societyEnglish`)

quantile(diff_a, prob = c(.5, .025, .75))
quantile(diff_b, prob = c(.5, .025, .75))

## Function
newdata = data.frame(std_count1 = median(melodic_df$std_count1),
                     std_count2 = median(melodic_df$std_count2),
                     functionalchange_bin = c("NF-NF", "F-F", "F-NF"),
                     std_semitonaldistance = c(0.1250),
                     society = "English", 
                     functional_total = 16240 + 16240)

pred.4 = predict(fit.4, newdata = newdata) %>% data.frame
