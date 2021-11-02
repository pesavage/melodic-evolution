suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(brms)
  library(patchwork)
})

## Read in Data and models
model_data = read.csv('data/model_data.csv')

fit.5.english <-
  brm(data = english_df, family = poisson,
      substitution_count ~ functional_change + 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      file = "results/functionplusdistance_english")

fit.5.japanese <-
  brm(data = japanese_df, family = poisson,
      substitution_count ~ functional_change + 
        semitonal_distance + 
        frequency1:frequency2:functional_total + (1|notepair),
      file = "results/functionplusdistance_japanese")


## Function to prepare data for plots
data_prep = function(d, society, model){
  
  # Get data from one society
  society_idx = model_data$society == society
  d = d[society_idx,]
  
  # Figure out note frequency medians from those societies
  d_median = median(c(d$frequency1, d$frequency2))
  d_total = floor(max(d$functional_total))
  
  # Create a new dataset to predict from
  newdata = data.frame(
    functional_change = rep(c("s", "w"), each = 11),
    semitonal_distance = rep(1:11, 2),
    notepair = "CD")
  
  newdata$society = society
  newdata$frequency1 = d_median
  newdata$frequency2 = d_median
  newdata$functional_total = d_total
  
  
  pred = predict(model, 
                   newdata = newdata) %>% data.frame
  pred = cbind(newdata, pred)
  
  pred
}

# Get plot data
plot_english = data_prep(model_data, "English", fit.5.english)
plot_japanese = data_prep(model_data, "Japanese", fit.5.japanese)


#### Plot ####
english_plot = 
  ggplot(data = 
  model_data[model_data$society == "English",], 
    aes(
      y = substitution_count,
      x = semitonal_distance,
      col = functional_change)) + 
  geom_jitter(width = 0.1, alpha = 0.7) + 
  geom_line(data = plot_english, aes(y = Estimate, 
                               x = semitonal_distance)) + 
  geom_ribbon(data = plot_english, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
              x = semitonal_distance,
              fill = functional_change), 
              linetype = 0,
              alpha = 0.2) + 
  ylab("Count of substitutions") + 
  xlab("Substitution distance (semitones)") + 
  theme_light() + 
  theme(legend.position = "bottom") + 
  scale_color_discrete(name="",
                       breaks=c("s", "w"),
                       labels=c("Strong Function", 
                                "Weak Function")) + 
  scale_fill_discrete(guide = "none")


japanese_plot = 
  ggplot(data = 
           model_data[model_data$society == "Japanese",], 
         aes(
           y = substitution_count,
           x = semitonal_distance,
           col = functional_change)) + 
  geom_jitter(width = 0.1, alpha = 0.7) + 
  geom_line(data = plot_japanese, aes(y = Estimate, 
                                     x = semitonal_distance)) + 
  geom_ribbon(data = plot_japanese, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
                  x = semitonal_distance,
                  fill = functional_change), 
              linetype = 0,
              alpha = 0.2) + 
  ylab("") + 
  xlab("Substitution distance (semitones)") + 
  theme_light() + 
  theme(legend.position = "bottom") + 
  scale_color_discrete(name="",
                       breaks=c("s", "w"),
                       labels=c("Strong Function", 
                                "Weak Function")) + 
  scale_fill_discrete(guide = "none")

p_out = english_plot + japanese_plot +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave(filename = "figures/example_mainfigure.jpeg", plot = p_out)