library(readxl)
library(stringr)

source('helper.R')

raw_data = read_xlsx("MelodicEvoSeq.xlsx", .name_repair = "universal")

# individual songs
raw_songs = raw_data[!duplicated(raw_data$PairNo),]

# by society & single songs
raw_english = raw_songs[raw_songs$Language == "English",]
raw_japanese = raw_songs[raw_songs$Language == "Japanese",]


end_dataframe = data.frame(
  note1              = character(),
  note2              = character(),
  society            = character(),
  functional_change  = factor(),
  mutation_count     = integer(), 
  count_1            = integer(),
  count_2            = integer(),
  functional_total   = integer()
)[1:132,]

#### Add note pairs ####
notes = c("C", "d", "D", "e", "E", "F", "g", "G", "a", "A", "b", "B")
all_pairs = t(combn(notes, 2))

end_dataframe$note1 = c(all_pairs[,1], all_pairs[,1])
end_dataframe$note2 = c(all_pairs[,2], all_pairs[,2])

#### Add Society ####
end_dataframe$society = rep(c("English", "Japanese"), each = 66)

#### Add default semi-tonal distances ####
# Load and Format semi-tonal data
semitonal_distance = read_xlsx("SubstitutionSize_distancematrices.xlsx", 
                               sheet = "Semitone distance matrix") %>%
   as.matrix()

semitonal_distance = semitonal_distance[,2:ncol(semitonal_distance)]
semitonal_distance = apply(semitonal_distance, 2, as.numeric)
rownames(semitonal_distance) = colnames(semitonal_distance)

semitonal_distance[upper.tri(semitonal_distance)] = 
  t(semitonal_distance)[upper.tri(semitonal_distance)]

semitonal_long = data.frame(col=colnames(semitonal_distance)[col(semitonal_distance)], 
                            row=rownames(semitonal_distance)[row(semitonal_distance)], 
                            semitonal_distance=c(semitonal_distance))

# Join datasets 
end_dataframe = left_join(end_dataframe, semitonal_long, by = c("note1" = "row", 
                                                                "note2" = "col"))

# Get counts of note pair substitutions
english_notecounts = get_notecounts(raw_english, notes)
japanese_notecounts = get_notecounts(raw_japanese, notes)

# Get counts of note occurrences 
unique_notes = function(x){
  z = unlist(strsplit(x = x, split = ""))
  paste(z, collapse = "")
}

english_notes = raw_english %>% 
  group_by(PairNo) %>% 
  summarize(notes = unique_notes(Full.note.sequence..aligned.)
              )
english_splitnotes = lapply(english_notes, function(x) 
  unlist(strsplit(x = x, split = "")))
english_total = unlist(english_splitnotes)
table(english_total)

japanese_notes = raw_japanese$Full.note.sequence..aligned.
japanese_splitnotes = lapply(japanese_notes, function(x) 
  unlist(strsplit(x = x, split = "")))
japanese_total = unlist(japanese_splitnotes)
table(japanese_total)


# Number of Notes in each song
n_full = str_length(raw_data$Full.note.sequence..unaligned.)
# number of ornamental notes in each song
n_ornamental = str_length(raw_data$Ornamental.notes)
# Number of final notes
n_final = str_length(raw_data$Final.note)
# Number of stressed notes
n_stressed = str_length(raw_data$Stressed.notes)
# number of unstressed notes
n_unstressed = n_full - (n_ornamental + n_final + n_stressed)

strong_function = str_length(raw_data$Final.mutations) + str_length(raw_data$Stress.mutations)

weak_function = 

mut<-s[,c(1,13:20)]
mut$nFull<-str_length(mut[,2])
mut$nOrn<-str_length(mut[,3])
mut$nFin<-str_length(mut[,4])
mut$nStr<-str_length(mut[,5])
mut$nUnStr<-mut$nFull-rowSums(mut[,11:13])
mut$nOrnMut<-str_length(mut[,6])
mut$nFinMut<-str_length(mut[,7])
mut$nStrMut<-str_length(mut[,8])
mut$nUnStrMut<-str_length(mut[,9])

mut$FinMutRate<-mut$nFinMut/mut$nFin
mut$StrMutRate<-mut$nStrMut/mut$nStr
mut$UnStrRate<-mut$nUnStrMut/mut$nUnStr
mut$OrnMutRate<-mut$nOrnMut/mut$nOrn
mut$StrongFunctionRate<-(mut$nFinMut+mut$nStrMut)/(mut$nFin+mut$nStr)
mut$WeakFunctionRate<-(mut$nUnStrMut+mut$nOrnMut)/(mut$nUnStr+mut$nOrn)