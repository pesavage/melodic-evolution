setwd("/Users/pesavage/Documents/Research/Papers/Unpublished/English and Japanese folk song evolution automated analysis")

#open packages
library(plotrix)
library(seqinr)
library(Biostrings)
library(seriation)
library(tidyr)
library(ggplot2)
library(dplyr)

#evolutionary rates of highly-related NHK ballad pairs
all.mut<-read.csv("AllMutationRatesNHK.csv") #use this for Japanese

#functional role
mut <- subset(all.mut, Recording=="Older")


#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#Make violin plot
data_wide <- mut[ , 1:4]
data_wide %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",size=.2) + geom_jitter(binaxis='y', stackdir='center', size=0.4,position=position_jitter(0.15)) + ylim(0,1) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t tests
t.test(mut[,4],mut[,3],alternative="greater",paired=TRUE)
t.test(mut[,3],mut[,2],alternative="greater",paired=TRUE)
t.test(mut[,2],mut[,1],alternative="greater",paired=TRUE)
colSums(!is.na(mut[1:6]))

#mutational distance
interval<-colSums(all.mut[,7:12],na.rm=TRUE)
cor.test(interval,c(2:7),method="spearman",alternative="less")
x <- barplot(interval,names.arg=c("2nd","3rd","4th","5th","6th","7th"),ylab=expression(paste("Number of substitutions")))

#mutational distance (grouped by # of semitones)
semitone<-colSums(all.mut[,27:37],na.rm=TRUE)
cor.test(semitone,c(1:11),method="spearman",alternative="less")

barplot(semitone,names.arg=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"),ylab=expression(paste("Number of substitutions")))


#突然変異立の計算（図３．４）evolutionary rates of highly-related Child ballad pairs
library(plotrix)
all.mut<-read.csv("AllMutationRates.csv")

#functional role
mut <- subset(all.mut, Recording=="Older")

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#Make violin plot
data_wide <- mut[ , 1:4]
data_wide %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",size=.2) + geom_jitter(binaxis='y', stackdir='center', size=0.4,position=position_jitter(0.15)) + ylim(0,2) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

t.test(mut[,4],mut[,3],alternative="greater",paired=TRUE)
t.test(mut[,3],mut[,2],alternative="greater",paired=TRUE)
t.test(mut[,2],mut[,1],alternative="greater",paired=TRUE)
colSums(!is.na(mut[1:6]))

#mutational distance
interval<-colSums(all.mut[,7:12],na.rm=TRUE)
cor.test(interval,c(2:7),method="spearman",alternative="less")
x <- barplot(interval,names.arg=c("2nd","3rd","4th","5th","6th","7th"),ylab=expression(paste("Number of substitutions")))

#mutational distance (grouped by # of semitones)
semitone<-colSums(all.mut[,87:97],na.rm=TRUE) #Note these columns are numbered differently from the Japanese columns
cor.test(semitone,c(1:11),method="spearman",alternative="less")

barplot(semitone,names.arg=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"),ylab=expression(paste("Number of substitutions")))

#########
##confirmatory analyses (of functional role and mutational distance): repeat above, changing definition of "mut" sample as follows:
old <- subset(all.mut, Recording=="Older")

mut <- subset(old, Period=="post-1907")
mut <- subset(old, Period=="pre-1907")
mut <- subset(old, Region=="Britain")
mut <- subset(old, Region=="Americas")
mut <- subset(all.mut, Recording=="Younger") #(this reverses the ancestry assumption)

#transmission fidelity
mut <- subset(all.mut, Recording=="Younger") #(NB: this uses younger variant to classify as oral vs. written. This assumption is reversed later [see Methods].)

#create subsets
oral <- subset(mut, Transmission=="Oral")
written <- subset(mut, Transmission=="Written")
my.values<-c(mean(oral[,6]), mean(written[,6]))
write.csv(my.values,"averages.csv")

err1<- c(std.error(oral[,6]), std.error(written[,6]))
x <- barplot(my.values, ylim=c(0,0.1),cex.names=1.2,cex.axis=1.2,cex.lab=1.2,names.arg=c("Oral","Written"),ylab=expression(paste("Mutation rate (",mu,")")))
arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)

t.test(oral[,6],written[,6],alternative="greater")

c(length(oral[,6]), length(written[,6]))

#confirmatory analyses (of functional role and mutational distance): repeat above, changing definition of "mut" sample as follows [note: this is the same as for functional position and mutational distance, except that "younger" and "older" are reversed for the ancestry assumption]
old <- subset(all.mut, Recording=="Younger")

mut <- subset(old, Period=="post-1907")
mut <- subset(old, Period=="pre-1907")
mut <- subset(old, Region=="Britain")
mut <- subset(old, Region=="Americas")
mut <- subset(all.mut, Recording=="Older")

#checking using only 34 pairs containing grace notes: same as above, but replace "all.mut" sample as follows:
all.mut<-read.csv("34GraceNotePairMutationRates.csv")

#compare indels vs. substitutions
chisq.test(c(564,368))
#X-squared = 41.219, df = 1, p-value = 1.361e-10




#（旋律配列と類似性計算)（図３．３))　melodic sequence alignment:

#Convert spreadsheet to .fasta by saving as WINDOWS-FORMATTED .txt file and then converting using this site: http://sequenceconversion.bugaco.com/converter/biology/sequences/tab_to_fasta.php



#load file

#load all 10,064 English and Japanese melodies  (without phrase gaps)
melodies<-read.fasta("10064MelodySequences.fasta", seqtype = "AA")

#restrict to only English melodies
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

ja<-melodies[4126:10064]

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







#to read back in distance matrix and labels:

#English folk songs
eng.d<-as.dist(as.matrix(read.table("EnglishDist.txt"))) 

#Reorder matrix to put tune families together (took ~2 minutes with 4125x4125 matrix using OLO method, only ~15sec with HC method) (code modified from Hahsler et al., 2008, J. Statistical Software)

##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
eng.res <- dissplot(eng.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(eng.res$x_reordered),"Englishreordered.txt")

#plot histogram
hist(eng.d)
abline(v=0.15, lty=2)


#Japanese folk songs
ja.d<-as.dist(as.matrix(read.table("JapanDist.txt"))) 


##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
ja.res <- dissplot(ja.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(ja.res$x_reordered),"Japanesereordered.txt")

#plot histogram
hist(ja.d)
abline(v=0.15, lty=2)






#Supplementary analyses

#統計学的有意の計算　significance (one pair only)
#load file
#full (no phrase gaps)
bronson<-read.fasta("Bronson4184NoPhraseGaps.fasta", seqtype = "AA")

#save sequences as strings
s3 <- c2s(bronson[[1]])
s4<- c2s(bronson[[2]])

#make sure all CAPS
s3<- toupper(s3)
s4<-toupper(s4)

(globalAligns3s4<-pairwiseAlignment(s3, s4, type="global", substitutionMatrix = NULL, gapOpening = -12,gapExtension = -6))

# Print out the optimal global alignment and its score

pid(globalAligns3s4,type="PID4")
#find percent identity


globalAligns3s4 <- pairwiseAlignment(s3, s4, type="global", substitutionMatrix = NULL, gapOpening = -12,gapExtension = -6, scoreOnly = TRUE)

generateSeqsWithMultinomialModel <- function(inputsequence, X)   
{     
 # Change the input sequence into a vector of letters      
require("seqinr") # This function requires the SeqinR package.      
inputsequencevector <- s2c(inputsequence)      # Find the frequencies of the letters in the input sequence "inputsequencevector":      
mylength <- length(inputsequencevector)      
mytable <- table(inputsequencevector)      # Find the names of the letters in the sequence      
letters <- rownames(mytable)      
numletters <- length(letters)      
probabilities <- numeric() # Make a vector to store the probabilities of letters      
for (i in 1:numletters)      
{         
letter <- letters[i]         
count <- mytable[[i]]        
 probabilities[i] <- count/mylength      
}      
# Make X random sequences using the multinomial model with probabilities "probabilities"      
seqs <- numeric(X)      
for (j in 1:X)      
{         
seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement         
seq <- c2s(seq)         
seqs[j] <- seq      
}      
# Return the vector of random sequences      
return(seqs)   
}


randomseqs <- generateSeqsWithMultinomialModel(s4,100) 
randomscores <- double(100) # Create a numeric vector with 100 elements
for (i in 1:100)
   {
      score <- pairwiseAlignment(s3, randomseqs[i], type="global", substitutionMatrix = NULL,
        gapOpening = -12, gapExtension = -6, scoreOnly = TRUE)
      randomscores[i] <- score
   }

hist(randomscores, col="red") # Draw a red histogram

sum(randomscores >= globalAligns3s4)
# [1] 0 
#i.e., 0/100 random sequences were greater than observed match score
#i.e., p<0.01


