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
                                  mutation_count=c(substitution_matrix))
  
  substitutions_long = substitutions_long[!is.na(substitutions_long$mutation_count),]
  substitutions_long$mutation_count[substitutions_long$col == substitutions_long$row] = 0

  list(substitutions = substitutions_long)  
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

## Total note counts
english_count = read.csv('results/english_notecounts.csv')
colnames(english_count) = c("note", "count")
english_count$society = "English"

japanese_count = read.csv('results/japanese_notecounts.csv')
colnames(japanese_count) = c("note", "count")
japanese_count$society = "Japanese"

counts = rbind(english_count, japanese_count)

# # mutability 
# english_mutability = read.csv('results/english_mutability.csv')
# colnames(english_mutability) = c("note", "mutability")
# english_mutability$society = "English"
# 
# japanese_mutability = read.csv('results/japanese_mutability.csv')
# colnames(japanese_mutability) = c("note", "mutability")
# japanese_mutability$society = "Japanese"
# 
# mutability = rbind(japanese_mutability, english_mutability)


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

# Add Mutability
# model_df = left_join(model_df, mutability, by = c("col" = "note", "society" = "society"))
# model_df = left_join(model_df, mutability, by = c("row" = "note", "society" = "society"))

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

## Calculate mutability
# english_unchanged = diag(as.matrix(english_substitutionmatrix[2:13,2:13]))
# english_changed = colSums(english_substitutionmatrix[,2:13], na.rm = TRUE) + english_substitutionmatrix[2:13,1]
# english_mutability = english_changed / english_changed + english_unchanged
# unchanged<-c(mat[2,2],mat[3,3],mat[4,4],mat[5,5],mat[6,6],mat[7,7],mat[8,8],mat[9,9],mat[10,10],mat[11,11],mat[12,12],mat[13,13])
# changed<-colSums(full.mat)[2:13] 
# total<-changed+unchangeda
# (mutability<-changed/total)

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

model_df = left_join

write.csv(model_df, "results/reviewer_modeldata.csv",
          row.names = FALSE, fileEncoding = 'utf-8')
