library(ggplot2)
library(dplyr)
library(brms)
library(patchwork)

model_data = read.csv('results/model_data.csv')

model_data$std_semitonaldistance = model_data$semitonal_distance / max(model_data$semitonal_distance)

fit.4 = brm(data = model_data, family = binomial,
            substitution_count | trials(functional_total) ~ 
              frequency1:frequency2:society + 
              functional_change + std_semitonaldistance,
            file = "results/fullrslopes_bin2")


#Replacing with figure showing something like x axis substitution distance and y axis substitutions/(notecount 1x notecount2), with points coloured by function and fitted  lines for function and non-function notes. Perhaps another panel with a breakdown between functional types...or the note final result?

english_idx = model_data$society == "English"
english_median = median(c(model_data$frequency1[english_idx],
                        model_data$frequency2[english_idx]))
japanese_median = median(c(model_data$frequency1[!english_idx],
                          model_data$frequency2[!english_idx]))

english_total = floor(max(model_data$functional_total[english_idx]))
japanese_total = floor(max(model_data$functional_total[!english_idx]))


newdata = data.frame(functional_change = rep(c("s", "w"), each = 11),
                     std_semitonaldistance = 
                       rep(c(0.09090909, 0.18181818, 0.27272727, 
                             0.36363636, 0.45454545, 0.54545455,
                             0.63636364, 0.72727273, 0.81818182,
                             0.90909091, 1.00), 2))

english_data = newdata
english_data$society = "English"
english_data$frequency1 = english_median
english_data$frequency2 = english_median
english_data$functional_total = english_total

japanese_data = newdata
japanese_data$society = "Japanese"
japanese_data$frequency1 = japanese_median
japanese_data$frequency2 = japanese_median
japanese_data$functional_total = japanese_total

newdata = rbind(english_data, japanese_data)

pred.4 = predict(fit.4, newdata = newdata) %>% data.frame
pred.4 = cbind(pred.4, newdata)
pred.4$semitones = 1:11

#### Plot ####
example_plot = ggplot(data = model_data, aes(
  y = substitution_count,
  x = semitonal_distance,
  col = functional_change)) + 
  geom_jitter(width = 0.1, alpha = 0.7) + 
  geom_line(data = pred.4, aes(y = Estimate, 
                               x = semitones)) + 
  geom_ribbon(data = pred.4, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
              x = semitones,
              fill = functional_change), 
              linetype = 0,
              alpha = 0.2) + 
  facet_wrap(~society, scales = 'free') +
  ylab("Count of substitutions") + 
  xlab("Substitution distance (semitones)") + 
  theme_light() + 
  theme(legend.position = "bottom") + 
  scale_color_discrete(name="",
                       breaks=c("s", "w"),
                       labels=c("Strong Function", 
                                "Weak Function")) + 
  scale_fill_discrete(guide = "none")

example_plot
ggsave(filename = "figures/example_mainfigure.jpeg", plot = example_plot)

# Semitonal distance
semitonal_data = newdata 

pred.4 = predict(fit.4, newdata = semitonal_data) %>% data.frame
pred.4 = cbind(pred.4, semitonal_data)
pred.4$semitones = 1:7

pred.4 = pred.4 %>% 
  group_by(semitones, society) %>% 
  summarise(Estimate = mean(Estimate), 
            Q2.5 = mean(Q2.5),
            Q97.5 = mean(Q97.5))

plot_substitutions = ggplot(pred.4, aes(x = semitones, y = Estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                width=.2, position=position_dodge(.9)) +
  geom_point(data = model_data, aes(y = substitution_count, x = semitonal_distance), col = "grey") +
  xlab("Semitonal distance") + 
  ggtitle("Predicted substitutions with increases in semitonal distance") + 
  theme_light() + 
  ylab("Predicted substitutions") +
  theme(text = element_text(size=20, family="serif")) + 
  facet_wrap(~society)

# Strong - Weak function
function_data = newdata 

pred.function = predict(fit.4, newdata = function_data) %>% data.frame
pred.function = cbind(pred.function, function_data)
pred.function$note_function = ifelse(pred.function$functional_change == "s",
                                     "Strong function", "Weak function")

plot_function = ggplot(pred.function, aes(x = note_function, y = Estimate)) +
  geom_boxplot() + 
  xlab("Note function") + 
  ggtitle("Predicted substitutions for strong and weak functional notes") + 
  theme_light() + 
  ylab("Predicted substitutions") + 
  theme(text = element_text(size=20, family="serif"),
        axis.title.y = element_blank()) + 
  facet_wrap(~society)

main_plot = plot_substitutions / plot_function

ggsave(plot = main_plot, filename = "figures/main_effects.jpeg")
