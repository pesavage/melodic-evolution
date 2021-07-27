## Custom functions for "Sequence alignment of folk song melodies reveals cross-cultural mechanisms of musical evolution"


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

get_notecounts = function(d, notes){
  n_notes = vector(length = length(notes))
  names(n_notes) = notes
  for(note in notes){
    note_occurance = 
      sum(
        str_count(d$Full.note.sequence..aligned., pattern = note),
        na.rm = TRUE)
    
    ornamental_mutations = 
      sum(
        str_count(d$Ornamental.mutations, pattern = note),
        na.rm = TRUE)
    
    final_mutations = 
      sum(
        str_count(d$Final.mutations, pattern = note),
        na.rm = TRUE)
    
    stress_mutations = 
      sum(
        str_count(d$Stress.mutations, pattern = note),
        na.rm = TRUE)
    
    unstressed_mutations =
      sum(
        str_count(d$Unstressed.mutations, pattern = note),
        na.rm = TRUE)
    
    n_notes[note] = note_occurance - 
      ornamental_mutations - 
      final_mutations - 
      stress_mutations -
      unstressed_mutations
  }
  
  n_notes
}
