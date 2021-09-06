# This script creates the data for the model in substitution_models.R from
# the raw data file used to align songs

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(dplyr)
  library(assertthat)
})

source('helper.R')

get_data = function(df){
  # Coding columns
  # Columns 111 to 245 contain information on:
  # Substitution changes between notes; 
  # the semitonal distance of that change; 
  # whether that change was a strong functional change or weak functional 
  # However, it does not give information on near semitonal distance
  # and weak functional change. We calculate this from existing data. 
  coding_idx = 111:245
  
  notes = c("C", "d",	"D",	"e",	"E",	"F",	"g",
            "G",	"a",	"A",	"b",	"B")
  
  single_song = df[!duplicated(df$PairNo),]
  
  #### Substitutions ####
  # Extract information on notes, function, and distance from titles
  parsed_titles = sapply(colnames(df)[coding_idx], function(c) 
    unlist(
      str_match_all(c, "(^[A-Za-z]{1})([A-Za-z]{1})([0-9]{1,2})([w|s]$)"
      )))
  
  # Occurances 
  occurances = colSums(single_song[,coding_idx], na.rm = TRUE)
  
  # Make into a single data frame
  model_matrix = cbind(t(parsed_titles), occurances)
  model_matrix = as.data.frame(model_matrix)
  colnames(model_matrix) = c("ID", "note1", "note2", "semitonal_distance", "functional_change", "substitution_count")
  
  model_matrix$substitution_count = as.numeric(model_matrix$substitution_count)
  model_matrix$semitonal_distance = as.numeric(model_matrix$semitonal_distance)
  
  assert_that(all(table(model_matrix$ID) == 1))
  
  # Not all far strong pairs are in the dataframe. We add them here.
  farweak_names = model_matrix$ID[model_matrix$semitonal_distance >= 6 & model_matrix$functional_change == "w"]
  farstrong_names = str_replace(farweak_names, "w", "s")
  
  assert_that(length(farstrong_names) == 66)
  
  farstrong_missing = farstrong_names[!farstrong_names %in% model_matrix$ID]
  farstrong_settings = sapply(farstrong_missing, function(c) 
    unlist(
      str_match_all(c, "(^[A-Za-z]{1})([A-Za-z]{1})([0-9]{1,2})([w|s]$)"
      )))
  
  farstrong = cbind(t(farstrong_settings), 0)
  farstrong = as.data.frame(farstrong)
  colnames(farstrong) = colnames(model_matrix)
  
  
  farstrong$substitution_count = as.numeric(farstrong$substitution_count)
  farstrong$semitonal_distance = as.numeric(farstrong$semitonal_distance)
  
  ## add to main data
  model_matrix = rbind(model_matrix, farstrong)
  
  assert_that(all(table(model_matrix$ID) == 1))
  
  # Near weak counts are: 
  # total_counts - far & weak counts - near & strong - far & strong
  # But we need to calculate them from existing other data
  # nearweak_names = 
  #   str_extract(model_matrix$ID, "^[A-Za-z]{2}[0-6]s") %>% 
  #   unique() %>% 
  #   na.omit %>% 
  #   str_replace_all(., "s", "w")
  
  # Columns titled with two note letter codes only, indicate total substitution counts
  # i.e. not subset by function type of semitonal distance. 
  totalsubstitution_idx = str_detect(colnames(single_song), "^[A-Za-z]{2}$")
  total_counts = colSums(single_song[,totalsubstitution_idx], na.rm = TRUE)
  
  # We calculate the subtotal for all other types fro model_matrix to create
  # a subtotal
  st = model_matrix %>% 
    group_by(note1, note2) %>% 
    summarise(subtotal = sum(substitution_count))
  
  subtotal = st$subtotal
  names(subtotal) = paste0(st$note1, st$note2)
  
  # Reorder notes to be the same as total_counts
  subtotal = subtotal[names(total_counts)]
  
  near_weak = total_counts - subtotal
  
  assert_that(all(near_weak >= 0))
  
  nearstrong_idx = str_detect(colnames(df), "[A-Za-z]{2}[0-6]s")
  nearstrong_names = colnames(df)[nearstrong_idx]
  nearweak_names = str_replace(nearstrong_names, "s", "w")
  
  nearweak_df = do.call("rbind",
                        str_match_all(nearweak_names, "(^[A-Za-z]{1})([A-Za-z]{1})([0-9]{1,2})([w|s]$)"
                        ))
  
  nearweak_df = as.data.frame(nearweak_df)
  colnames(nearweak_df) = c("ID", "note1", "note2", "semitonal_distance", "functional_change")
  nearweak_df$substitution_count = near_weak
  
  # Some notes have the same near and far distance (6 semitones up & down).
  # We only count these once. 
  nearweak_df = nearweak_df[!nearweak_df$ID %in% 
                              c("Cg6w", "Da6w", "eA6w", 
                                "Eb6w", "FB6w", "dG6w"),]
  
  model_matrix = rbind(model_matrix, nearweak_df)
  
  # All note ids should occur once
  assert_that(all(table(model_matrix$ID) == 1))
  
  #### Total counts #### 
  # The baseline model is the product of note occurrences. 
  # We need a count of note occurrences. 
  # Note counts are determined by counting the number of times notes occur
  # in one of the song pairs + insertions between the songs
  
  song_counts = 
    sapply(notes, function(n) {
      total = sum(str_count(df$Full.note.sequence..unaligned., n),
                  na.rm = TRUE)
      ornamental_mutations = sum(str_count(df$Ornamental.mutations, n),
                                 na.rm = TRUE)
      final_mutations = sum(str_count(df$Final.mutations, n),
                            na.rm = TRUE)
      stress_mutations = sum(str_count(df$Stress.mutations, n),
                             na.rm = TRUE)
      unstressed_mutations = sum(str_count(df$Unstressed.mutations, n),
                                 na.rm = TRUE)
      # Total - all mutations gives the note count
      total - (ornamental_mutations + final_mutations + stress_mutations + unstressed_mutations)
    })
  
  insertion_cols = str_detect(colnames(df), "^[A-Za-z]{1}\\.$")
  insertion_counts = colSums(df[,insertion_cols], na.rm = TRUE)
  
  note_totals = song_counts + insertion_counts
  
  notecounts_df = data.frame(note = names(note_totals), 
                                  value = note_totals)
  
  ## Add to the dataframe
  model_matrix = left_join(model_matrix, 
                           notecounts_df, by = c("note1" = "note")) %>% 
    left_join(., notecounts_df, by = c("note2" = "note"))
  
  
  model_matrix = 
    rename(model_matrix, frequency1 = value.x, frequency2 = value.y)
  
  # Functional totals
  # We divide the count of a substitution by the frequency of a particular type
  # To ensure that we are not only finding effects due to a difference in 
  # base-rate of types. 
  
  total_notes = sum(
    sapply(df$Full.note.sequence..unaligned., nchar), 
    na.rm = TRUE)
  
  ornamental_notes = sum(
    sapply(df$Ornamental.notes, nchar), 
    na.rm = TRUE)
  
  final_notes = sum(
    sapply(df$Final.note, nchar), 
    na.rm = TRUE)
  
  stressed_notes = sum(
    sapply(df$Stressed.notes, nchar), 
    na.rm = TRUE)
  
  unstressed_notes = total_notes - ornamental_notes - final_notes - stressed_notes
  
  strong_function = final_notes + stressed_notes
  weak_function   = ornamental_notes + unstressed_notes
  
  typetotal_df = data.frame(functional_change = c("w", "s"), 
                            functional_total = c(weak_function, strong_function))
  
  # # dividing by two accounts for the fact that one substitution 
  # # involves two melodies
  # typetotal_df$functional_total = typetotal_df$functional_total / 2
  
  model_matrix = left_join(model_matrix, typetotal_df, by = "functional_change")
  
  
  ## Standardization
  model_matrix$frequency1 = 
    model_matrix$frequency1 / 
    max(model_matrix[,c("frequency1", "frequency2")])
  model_matrix$frequency2 = 
    model_matrix$frequency2 / 
    max(model_matrix[,c("frequency1", "frequency2")])
  
  model_matrix$functionaltotal_normal = model_matrix$functional_total
  model_matrix$functional_total = 
    model_matrix$functional_total / max(model_matrix$functional_total)
  
  model_matrix
}

## Run function
raw_data = read_xlsx("MelodicEvoSeq.xlsx", .name_repair = "universal")

english_raw  = raw_data[raw_data$Language == "English",]
japanese_raw = raw_data[raw_data$Language == "Japanese",]

english_modeldata  = get_data(english_raw)
japanese_modeldata = get_data(japanese_raw)

model_data = rbind(english_modeldata, japanese_modeldata)
model_data$society = rep(c("English", "Japanese"), 
                         each = nrow(english_modeldata))

model_data$notepair = paste0(model_data$note1, model_data$note2)

# We should not have any semitonal distances below 16
assert_that(!any(model_data$semitonal_distance == 16))

# Each society should occur 264 times
assert_that(all(table(model_data$society) == 252))

# There are 66 note pairs for the 8 different categories:
# 1) English & Japanese societies, 
# 2) Close and far semitone distances, 
# 3) Strong and weak functional notes
# 6 note pairs only have one semitonal distance (6 semitones), because 
# the note pair has the same distance up and down the scale. Which occur 
# twice in each society, four times total
assert_that(nrow(model_data) == (66 * 8) - 4 * 6)

write.csv(model_data, 'data/model_data.csv')
