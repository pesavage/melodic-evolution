library(readxl)
library(dplyr)
library(ggplot2)

source('helper.R')

## Data needs these columns:
# Total = Frequency of note change (bidirectional)
# DistanceSemitone = Semitone distance
# NoteTotal1 = Frequency of the first note in a pair
# NoteTotal2 = Frequency  of the second note in a pair
# MinTotal = Minimum of NoteTotal1 and NoteTotal2

english_substitutionmatrix = read.csv('results/english_SubstitutionMatrix.csv', 
                               row.names = 1, nrows = 13) # last row is mutability 
english_output = get_substitutions(english_substitutionmatrix)

japanese_substitutionmatrix = read.csv('results/japanese_SubstitutionMatrix.csv', 
                                       row.names = 1, nrows = 13) # last row is mutability 
japanese_output = get_substitutions(japanese_substitutionmatrix)

## Combine English & Japanese
substitutions_long = rbind(english_output$substitutions, japanese_output$substitutions)
substitutions_long$society = rep(c("English", "Japanese"), each = nrow(english_output$substitutions))

## Total note counts
english_count = read.csv('results/english_notecounts.csv')
colnames(english_count) = c("note", "count")
english_count$society = "English"

japanese_count = read.csv('results/japanese_notecounts.csv')
colnames(japanese_count) = c("note", "count")
japanese_count$society = "Japanese"

counts = rbind(english_count, japanese_count)

# Load and Format semi-tonal data
semitonal_distance = read_xlsx("SubstitutionSize_distancematrices.xlsx", 
                               sheet = "Semitone distance matrix") %>%
  as.matrix()

semitonal_distance = semitonal_distance[,2:ncol(semitonal_distance)]
semitonal_distance = apply(semitonal_distance, 2, as.numeric)
rownames(semitonal_distance) = colnames(semitonal_distance)

semitonal_long = data.frame(col=colnames(semitonal_distance)[col(semitonal_distance)], 
                            row=rownames(semitonal_distance)[row(semitonal_distance)], 
                            semitonal_distance=c(semitonal_distance))

# Join datasets 
model_df = left_join(substitutions_long, semitonal_long, by = c("row", "col"))
model_df = model_df %>% filter(col != "-") %>% filter(row != "-")

# Add frequency
model_df = left_join(model_df, counts, by = c("col" = "note", "society" = "society")) 
model_df = left_join(model_df, counts, by = c("row" = "note", "society" = "society"))

colnames(model_df) = c("note1", "note2", "mutation_count", "society",
                       "semitonal_distance", "count_1", "count_2")

model_df$minimum_count = apply(model_df[,c("count_1", "count_2")], 1, min)

model_df = model_df[model_df$note1 != model_df$note2,]

# remove notes that don't occur
idx = model_df$count_1 == 0 | model_df$count_2 == 0
model_df = model_df[!idx,]

# Make semi-tonal adjustments based on manual calculations
new_semitones = read_xlsx('data/SubstitutionMatrices.xlsx', 
                          sheet = "Large intervals")

for(i in 1:nrow(new_semitones)){
  row = new_semitones[i,]
  existingrow_idx = which(row$Note1 == model_df$note1 & row$Note2 == model_df$note2)
  existing_row = model_df[existingrow_idx,]
  
  # reduce mutation_count
  existing_row$mutation_count[existing_row$society == "Japanese"] = 
    existing_row$mutation_count[existing_row$society == "Japanese"] - row$ja_substitutions
  existing_row$mutation_count[existing_row$society == "English"] = 
    existing_row$mutation_count[existing_row$society == "English"] - row$eng_substitutions
  
  # create new row 
  new_row = existing_row
  new_row$mutation_count[new_row$society == "Japanese"] = row$ja_substitutions
  new_row$mutation_count[new_row$society == "English"] = row$eng_substitutions
  new_row$semitonal_distance[new_row$society == "Japanese"] = row$semitone_distance
  new_row$semitonal_distance[new_row$society == "English"] = row$semitone_distance
  
  # put data back
  model_df[existingrow_idx,] = existing_row
  model_df = rbind(model_df, new_row)
}

functional_substitutions = read_xlsx('data/SubstitutionMatrices.xlsx', 
                                sheet = "Functional intervals")

model_df$functional_change = "NF-NF"

for(i in 1:nrow(functional_substitutions)){
  row = functional_substitutions[i,]
  
  existingrow_idx = which(row$Note1 == model_df$note1 & 
                            row$Note2 == model_df$note2 & 
                            row$semitone_distance == model_df$semitonal_distance &
                            model_df$functional_change == "NF-NF")
  existing_row = model_df[existingrow_idx,]
  
  # Change old substitutions
  # reduce mutation_count
  row_jmatch = existing_row$society == "Japanese" 
  row_ematch = existing_row$society == "English"
  existing_row$mutation_count[row_jmatch] = 
    existing_row$mutation_count[row_jmatch] - row$ja_substitutions
  existing_row$mutation_count[row_ematch] = 
    existing_row$mutation_count[row_ematch] - row$eng_substitutions
  
  
  
  # new row for F-F changes
  new_row = existing_row
  new_row$mutation_count[new_row$society == "Japanese"] = row$ja_substitutions
  new_row$mutation_count[new_row$society == "English"] = row$eng_substitutions
  new_row$semitonal_distance[new_row$society == "Japanese"] = row$semitone_distance
  new_row$semitonal_distance[new_row$society == "English"] = row$semitone_distance
  
  new_row$functional_change = ifelse(row$functional_notes == 1, "F-F", "F-NF")
  
  # Put rows back
  model_df[existingrow_idx,] = existing_row
  model_df = rbind(model_df, new_row)
  
  if(any(existing_row$mutation_count < 0)) stop("STOP")
}


assertthat::assert_that(nrow(model_df[model_df$note1 == "C" &
                                        model_df$note2 == "D" &
                                        model_df$semitonal_distance == 2 & 
                                        model_df$society == "English",]) == 3)

assertthat::assert_that(all(model_df$mutation_count[model_df$note1 == "C" &
                                                  model_df$note2 == "D" &
                                                  model_df$semitonal_distance == 2 & 
                                                  model_df$society == "English"] == c(55, 2, 3)))

## Total notes in each society 
total_notes = data.frame(society = c("English", "Japanese"),
                         total_notes = c(sum(english_count$count), sum(japanese_count$count)),
                         total_substitutions = tapply(model_df$mutation_count, model_df$society, sum))

model_df = left_join(model_df, total_notes, "society")


# total notes by function
function_total = data.frame(
  society = c(rep("English", 3), rep("Japanese", 3)),
  functional_change = rep(c("NF-NF", "F-F", "F-NF"), 2),
  functional_total = c(
    16240 + 16240, 5039 + 5039, 5039 + 16240,
    7424 + 7424, 2406 + 2406, 2406 + 7424
  )
)

model_df = left_join(model_df, function_total, 
                     by = c("society", "functional_change"))

model_df = model_df[model_df$functional_change != "F-NF",]
model_df$functionalchange_bin = ifelse(model_df$functional_change == "NF-NF", 0, 1)
model_df$functional_total = model_df$functional_total / 2

write.csv(model_df, "results/reviewer_wonf_f_modeldata.csv",
          row.names = FALSE, fileEncoding = 'utf-8')

# Make main figure
model_df$functional_change = ifelse(model_df$functional_change == "NF-NF", "Weaker function",
                                    "Strong function")
model_df$functionalchange_bin = ifelse(model_df$functional_change == "NF-NF", 0, 1)

plot_1 = ggplot(model_df, aes(y = mutation_count, 
                       x = count_1 * count_2, size = semitonal_distance)) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 0 + x) + 
  geom_point(shape = 21, aes(fill = functional_change), alpha = 0.8) +
  ylab("Number of Subsitutions") + 
  xlab("Note 1 count x Note 2 count") +
  facet_wrap(~society, scales = "free") + 
  scale_x_continuous(labels = scales::comma) +
  theme(legend.position = c(0.12, 0.94), legend.title = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.0)),
        text = element_text(size=14)) +
  scale_size_continuous(guide = "none") +
  ggtitle("Substitutions against count interaction", "English and Japanese songs")

ggsave(filename = 'figures/model_plot.png', plot = plot_1)


# With model response


plot_2 = ggplot(model_df, aes(y = mutation_count / functional_total, 
                              x = count_1 * count_2, size = semitonal_distance)) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 0 + x, col = "#6DA0FD") + 
  geom_point(shape = 21, aes(fill = functional_change), alpha = 0.8) +
  ylab("Number of Subsitutions") + 
  xlab("Note 1 count x Note 2 count") +
  facet_wrap(~society, scales = "free_x") + 
  scale_x_continuous(labels = scales::comma) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.0)),
        text = element_text(size=14)) +
  guides(size = guide_legend(override.aes = list(linetype = 0)))  +
  ggtitle("Substitutions against count interaction", "English and Japanese songs")

# annotate
ann_text1 <- data.frame(mutation_count = 55,
                       functional_total = 16240,
                       count_1 = 2669 + 490,
                       count_2 = 1501,
                       semitonal_distance = 2,
                       society = "English",
                       lab = "C-D \n 2 Semitones")

ann_text2 <- data.frame(mutation_count = 23,
                        functional_total = 7424,
                        count_1 = 1144 + 175,
                        count_2 = 876,
                        semitonal_distance = 2,
                        society = "Japanese",
                        lab = "C-D \n 2 Semitones")

ann_text3 <- data.frame(mutation_count = 1,
                        functional_total = 7274,
                        count_1 = 1144 + 175,
                        count_2 = 876,
                        semitonal_distance = 10,
                        society = "Japanese",
                        lab = "D-C \n 10 Semitones")

ann_text4 <- data.frame(mutation_count = 0,
                        functional_total = 16240,
                        count_1 = 2669 + 490,
                        count_2 = 1501,
                        semitonal_distance = 10,
                        society = "English",
                        lab = "D-C \n 10 Semitones")

plot_2 = plot_2 + 
  geom_label(data = ann_text1, label = ann_text1$lab, show.legend = FALSE, size = 2, col = "#F13C2C") + 
  geom_label(data = ann_text2, label = ann_text2$lab, show.legend = FALSE, size = 2, col = "#F13C2C") +
  geom_label(data = ann_text3, label = ann_text3$lab, show.legend = FALSE, size = 2, col = "#F13C2C") +
  geom_label(data = ann_text4, label = ann_text3$lab, show.legend = FALSE, size = 2, col = "#F13C2C")  

plot_2
ggsave(filename = 'figures/dividedresponse_plot.png', plot = plot_2)

