library(dplyr)
library(ggplot2)
library(raincloudplots)


#### Data ####

melodic_df = read.csv('results/reviewer_modeldata.csv')

english_mutability = read.csv('results/english_mutability.csv')
colnames(english_mutability) = c("note", "mutability")
english_mutability$society = "English"

japanese_mutability = read.csv('results/japanese_mutability.csv')
colnames(japanese_mutability) = c("note", "mutability")
japanese_mutability$society = "Japanese"

mutability = rbind(japanese_mutability, english_mutability)

#melodic_df = left_join(melodic_df, mutability, )

#### Substitution distance and frequency ####

melodic_summary = melodic_df %>% 
  group_by(semitonal_distance, society) %>% 
  summarise(mean(mutation_count))

m_df = melodic_df[melodic_df$semitonal_distance <= 7,]
ggplot(m_df, aes(x = semitonal_distance, y = mutation_count, 
                       group = semitonal_distance, fill = factor(semitonal_distance))) + 
  ## add half-violin from {ggdist} package
  geom_boxplot(outlier.shape = NA, width=0.25) + 
  geom_jitter(aes(fill = factor(semitonal_distance)), 
              pch = 21, width = 0.2, height = 0) +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) + 
  theme(legend.position = "none")


##

melodic_df$label = paste0(melodic_df$note1, "-", melodic_df$note2)
melodic_df$functional_col = ifelse(melodic_df$functional_notes == 1, "red", "grey")

mde = melodic_df[melodic_df$society == "English",]

ggplot(melodic_df, aes(y = mutation_count, 
                       x = count_1 * count_2,
                       label = label, 
                       size = semitonal_distance)) +
  geom_smooth(method='lm', se = FALSE) + 
  geom_point(aes(col = functional_col), data = melodic_df) +
  theme(legend.position = "none") +
  ylab("Number of Subsitutions") + xlab("Note 1 count x Note 2 count") +
  facet_wrap(~society, scales = "free") + 
  ggtitle("Substitutions against count interaction in English and Japanese songs")


### 
library(readxl)

functional_distance = read_xlsx('data/SubstitutionMatrices.xlsx', 
                                sheet = "Functional intervals")

english_count = read.csv('results/english_notecounts.csv')
colnames(english_count) = c("note", "count")

d = left_join(functional_distance, english_count, by = c("Note1" = "note")) 
d = left_join(d, english_count, by = c("Note2" = "note"))

d$functional_col = ifelse(d$functional_notes == 1, "red", "grey")

ggplot(d, aes(y = eng_substitutions, 
                       x = count.x * count.y,
                       size = semitone_distance)) +
  geom_smooth(method='lm', se = FALSE) + 
  geom_point(aes(col = functional_col), data = d) +
  theme(legend.position = "none") +
  ylab("Number of Subsitutions") + xlab("Note 1 count x Note 2 count") + 
  ggtitle("Substitutions against count interaction")
