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
get_substitutions = function(substitution_matrix){
  substitution_matrix = as.matrix(substitution_matrix)
  substitution_matrix = substitution_matrix[,1:ncol(substitution_matrix)]
  colnames(substitution_matrix)[1] = "-"
  
  substitutions_long = data.frame(col=colnames(substitution_matrix)[col(substitution_matrix)], 
                                  row=rownames(substitution_matrix)[row(substitution_matrix)], 
                                  change_frequency=c(substitution_matrix))
  
  substitutions_long = substitutions_long[!is.na(substitutions_long$change_frequency),]
  substitutions_long$change_frequency[substitutions_long$col == substitutions_long$row] = 0
  
  unchanged_sites = data.frame(notes = colnames(substitution_matrix), 
                               unchanged_sites = diag(substitution_matrix))

  list(substitutions = substitutions_long, unchanged_sites = unchanged_sites)  
}

english_substitutionmatrix = read.csv('results/english_SubstitutionMatrix.csv', 
                               row.names = 1, nrows = 13) # last row is mutability 
english_output = get_substitutions(english_substitutionmatrix)

japanese_substitutionmatrix = read.csv('results/japanese_SubstitutionMatrix.csv', 
                                       row.names = 1, nrows = 13) # last row is mutability 
japanese_output = get_substitutions(japanese_substitutionmatrix)

## Combine English & Japanese

substitutions_long = rbind(english_output$substitutions, japanese_output$substitutions)
substitutions_long$society = rep(c("English", "Japanese"), each = nrow(english_output$substitutions))

unchanged_sites = rbind(english_output$unchanged_sites, japanese_output$unchanged_sites)
unchanged_sites$society = rep(c("English", "Japanese"), each = nrow(english_output$unchanged_sites))

# total note counts
unchanged_sites$total_count = NA
for(i in 1:nrow(unchanged_sites)){
  note = unchanged_sites$notes[i]
  society = unchanged_sites$society[i]
  
  changed_sites = 
  
}

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

model_df = left_join(model_df, unchanged_sites, by = c("col" = "notes", "society" = "society")) 
model_df = left_join(model_df, unchanged_sites, by = c("row" = "notes", "society" = "society"))

colnames(model_df) = c("note1", "note2", "change_frequency", "society",
                       "semitonal_distance", "frequency_1", "frequency_2")

model_df$minimum_frequency = apply(model_df[,c("frequency_1", "frequency_2")], 1, min)

model_df = model_df[model_df$note1 != model_df$note2,]

total_changes = tapply(model_df$change_frequency, model_df$society, sum)

model_df$total_changes = rep(total_changes, each = 66)

# remove notes that don't occur
idx = model_df$frequency_1 == 0 | model_df$frequency_2 == 0
model_df = model_df[!idx,]

write.csv(model_df, "results/reviewer_modeldata.csv", 
          row.names = FALSE, fileEncoding = 'utf-8')
