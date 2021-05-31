library(readxl)
library(dplyr)

## Data needs these columns:
# Total = Frequency of note change (bidirectional)
# DistanceSemitone = Semitone distance
# NoteTotal1 = Frequency of the first note in a pair
# NoteTotal2 = Frequency  of the second note in a pair
# MinTotal = Minimum of NoteTotal1 and NoteTotal2

# Loading the substitution matrix containing the frequencies 
# and substitutions 
substitution_matrix = read.csv('results/english_SubstitutionMatrix.csv', 
                               row.names = 1, nrows = 13) # last row is mutability 
substitution_matrix = as.matrix(substitution_matrix)
substitution_matrix = substitution_matrix[,1:ncol(substitution_matrix)]
colnames(substitution_matrix)[1] = "-"

substitutions_long = data.frame(col=colnames(substitution_matrix)[col(substitution_matrix)], 
                                row=rownames(substitution_matrix)[row(substitution_matrix)], 
                  change_frequency=c(substitution_matrix))

substitutions_long = substitutions_long[!is.na(substitutions_long$change_frequency),]
substitutions_long$change_frequency[substitutions_long$col == substitutions_long$row] = 0

unchanged_sites = data.frame(notes = colnames(substitution_matrix), 
                             total_frequency = diag(substitution_matrix))

# Load and Format semi-tonal data
semitonal_distance = read_xlsx("SubstitutionSize_distancematrices.xlsx", 
                               sheet = "Semitone distance matrix") %>%
  as.matrix()
rownames(semitonal_distance) = semitonal_distance[,1]
semitonal_distance = semitonal_distance[,2:ncol(semitonal_distance)]

semitonal_long = data.frame(col=colnames(semitonal_distance)[col(semitonal_distance)], 
                            row=rownames(semitonal_distance)[row(semitonal_distance)], 
                            semitonal_distance=c(semitonal_distance))

# Join datasets 

model_df = left_join(substitutions_long, semitonal_long, by = c("row", "col"))
model_df = model_df %>% filter(col != "-") %>% filter(row != "-")

model_df = left_join(model_df, unchanged_sites, by = c("row" = "notes")) %>%
  left_join(., unchanged_sites, by = c("col" = "notes"))

colnames(model_df) = c("note1", "note2", "change_frequency", 
                       "semitonal_distance", "frequency_1", "frequency_2")

model_df[,3:6] = apply(model_df[,3:6], 2, as.numeric)

model_df$minimum_frequency = min(model_df$frequency_1, model_df$frequency_2)

write.csv(model_df, "results/reviewer_modeldata.csv")
