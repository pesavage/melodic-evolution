#subscript of MelodicEvo.R

full<-read.csv("MelodicEvoSeq.csv",header=TRUE) #Import all 10,000+ sequences
#d<-read.csv("MelodicEvolutionFullSequenceAndMetadata.csv",headers=TRUE)

#load file

#load all ~10,000 English and Japanese melodies  (without phrase gaps)
melodies<-dataframe2fas(full[,c(1,14)], file="melodies.fasta")
melodies<-seqinr::read.fasta("melodies.fasta", seqtype = "AA")
#melodies<-seqinr::read.fasta("10064MelodySequences.fasta", seqtype = "AA")

#restrict to only those melodies previously identified to be highly related
rel<-melodies[1:656]

#Calculate similarities between all pairs using gap opening penalty (GOP) = .8 and gap extension penalty (GEP) = .2, including mode (best for within-family alignment in Savage & Atkinson 2015)

all_pairs <- combn(1:length(rel), 2)

align_from_index <- function(rel, index.a, index.b){
  seq1 <- c2s(rel[[index.a]])
  seq2 <- c2s(rel[[index.b]])
  return(pid(pairwiseAlignment(seq1, seq2,type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2) ,type="PID4"))
}

res <- apply(all_pairs, 2, function(indices) align_from_index(rel, indices[1], indices[2]) )

attributes(res) <- attributes(dist(1:length(rel))) 
res <- as.matrix(res)

rel.dist<-as.dist(1-(res/100))
write.table(as.matrix(rel.dist),"HighlyRelatedDist.txt")

#restrict to only English melodies (This uses old fasta file - update)
eng<-melodies[1:4125]

#Calculate similarities between all 8.5 million pairs using gap opening penalty (GOP) = .8 and gap extension penalty (GEP) = .2, including mode (best for within-family alignment in Savage & Atkinson 2015 - takes over a week!)

all_pairs <- combn(1:length(eng), 2)

align_from_index <- function(eng, index.a, index.b){
  seq1 <- c2s(eng[[index.a]])
  seq2 <- c2s(eng[[index.b]])
  return(pid(pairwiseAlignment(seq1, seq2,type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2) ,type="PID4"))
}

res <- apply(all_pairs, 2, function(indices) align_from_index(eng, indices[1], indices[2]) )

attributes(res) <- attributes(dist(1:length(eng))) 
res <- as.matrix(res)

eng.dist<-as.dist(1-(res/100))
write.table(as.matrix(eng.dist),"EnglishDist.txt")


#Now the same for only Japanese folk songs

ja<-melodies[4126:10062]

#Calculate similarities between all 8.5 million pairs using gap opening penalty (GOP) = .8 and gap extension penalty (GEP) = .2, including mode (best for within-family alignment in Savage & Atkinson 2015 - takes several weeks!)

all_pairs <- combn(1:length(ja), 2)

align_from_index <- function(ja, index.a, index.b){
  seq1 <- c2s(ja[[index.a]])
  seq2 <- c2s(ja[[index.b]])
  return(pid(pairwiseAlignment(seq1, seq2,type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2) ,type="PID4"))
}

res <- apply(all_pairs, 2, function(indices) align_from_index(ja, indices[1], indices[2]) )

attributes(res) <- attributes(dist(1:length(ja))) 
res <- as.matrix(res)

ja.dist<-as.dist(1-(res/100))
write.table(as.matrix(ja.dist),"JapanDist.txt")







