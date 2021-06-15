# Notes
# A B C

set.seed(543)

notes = LETTERS[1:3]
probability_change = matrix(c(1.00, 0.25, 0.25,
                              0.25, 1.00, 0.50,
                              0.25, 0.50, 1.00 ), 
                            nrow = 3, 
                            dimnames = list(notes, notes))

frequency = c(0.01, 0.25, 0.74)
names(frequency) = notes

n_notes = 50 # number of notes in a song
fake_songs = matrix(NA, nrow = n_notes, ncol = n_notes)
for(i in 1:n_notes){
  fake_songs[i,] = sample(notes, n_notes, replace = TRUE, prob = frequency)
}
  
of = apply(fake_songs, 1, table)
original_frequencies = matrix(0, nrow = n_notes, ncol = length(notes), dimnames = list(c(), notes))
for(i in seq_along(of)) original_frequencies[i,names(of[[i]])] = of[[i]]

all(rowSums(original_frequencies) == 50)
# Original frequency of notes
colSums(original_frequencies)

iterations = 100
changes = matrix(NA, nrow = iterations, ncol = length(notes))
mutated_songs = fake_songs
note_frequencies = list()
for(i in 1:iterations){ # number of generations
  for(j in 1:nrow(mutated_songs)){ # number of songs
    song = mutated_songs[j,]    
    mutated_song = sapply(song, function(s){
    sample(notes, 1, prob = probability_change[,s])
    })
    mutated_songs[j,] = mutated_song
  }
  nf = t(apply(mutated_songs, 1, table))
  note_frequencies[[i]] = colSums(nf)
}


# Original frequency of notes
colSums(original_frequencies)
# New Frequencies
note_frequencies[[length(note_frequencies)]]

lastgeneration_frequency = note_frequencies[[length(note_frequencies)]]

outer(lastgeneration_frequency, lastgeneration_frequency, "/")

head(note_frequencies)



