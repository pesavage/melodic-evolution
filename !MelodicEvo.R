setwd("/Users/pesavage/Documents/Research/Papers/Unpublished/English and Japanese folk song evolution automated analysis")

#open packages (install first as required - #If Biostrings is not yet installed, enter the following commented out code: 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
library(plotrix)
library(seqinr)
library(Biostrings) 

BiocManager::install("Biostrings")If not yet installed, follow installation instructions at 
library(seriation)
library(tidyr)
library(ggplot2)
library(dplyr)
library(varhandle)

#source("PID.R") #uncomment this to recalculate distance matrices of PID among the 10,064 melodic variants, but be aware that this will take up to a month on a standard computer!

#read back in distance matrices calculated by "PID.R" subscript and identify highly related (>=85%) melodic variants:

#English folk songs
eng.d<-as.dist(as.matrix(read.table("EnglishDist.txt"))) 

#plot histogram
hist(eng.d)
abline(v=0.15, lty=2)

#sort by similarity
eng.m<-as.matrix(eng.d)
xy <- t(combn(colnames(eng.m), 2))
d<-data.frame(xy, dist=eng.m[xy])
d<-unfactor(d) #fixes ID numbers being coded as factors rather than numbers
r<-subset(d,dist<=.15)

attach(r)
r.ord <- r[order(dist),]
detach(r)
length(r.ord$dist) #774 highly related pairs
length(unique(c(r.ord[,1],r.ord[,2])))#651 unique melodic variants

#filter to remove identical pairs
r.ord<-subset(r.ord,dist>0)
length(r.ord$dist) #reduces to 760 non-identical pairs)

#filter to only use more related pairs in cases of melodies highly related to multiple other melodies
r.ord<-r.ord %>% distinct(X1, .keep_all = TRUE)
length(r.ord$dist) #reduces to 383 highly related pairs not with same X1 melody

r.ord<-r.ord %>% distinct(X2, .keep_all = TRUE)
length(r.ord$dist) #reduces to 316 highly related pairs also not with same X2 melody

#Remove remaining duplicates across two columns
dup<-duplicated(c(r.ord[,1],r.ord[,2]))[(length(r.ord$dist)+1):(length(r.ord$dist)*2)]
r.ord<-r.ord[!dup, ]
length(r.ord$dist) #reduces to 242 highly related unique [non-shared] pairs)


write.csv(r.ord,file="EngHighlyRelatedPairs.csv")



##########

#Japanese folk songs
ja.d<-as.dist(as.matrix(read.table("JapanDist.txt"))) 

#plot histogram
hist(ja.d)
abline(v=0.15, lty=2)

#sort by similarity
ja.m<-as.matrix(ja.d)
xy <- t(combn(colnames(ja.m), 2))
d<-data.frame(xy, dist=ja.m[xy])
d<-unfactor(d) #fixes ID numbers being coded as factors rather than numbers
r<-subset(d,dist<=.15)

attach(r)
r.ord <- r[order(dist),]
detach(r)
length(r.ord$dist) #143 highly related pairs
length(unique(c(r.ord[,1],r.ord[,2])))#220 unique melodic variants

#filter to remove identical pairs
r.ord<-subset(r.ord,dist>0)
length(r.ord$dist) #reduces to 131 non-identical pairs)

#filter to only use more related pairs in cases of melodies highly related to multiple other melodies
r.ord<-r.ord %>% distinct(X1, .keep_all = TRUE)
length(r.ord$dist) #reduces to 110 highly related pairs not with same X1 melody

r.ord<-r.ord %>% distinct(X2, .keep_all = TRUE)
length(r.ord$dist) #reduces to 101 highly related pairs not with same X1 melody

#Remove remaining duplicates across two columns
dup<-duplicated(c(r.ord[,1],r.ord[,2]))[(length(r.ord$dist)+1):(length(r.ord$dist)*2)]
r.ord<-r.ord[!dup, ]
length(r.ord$dist) #reduces to 91 highly related unique [non-shared] pairs)


write.csv(r.ord,file="JaHighlyRelatedPairs.csv")

#Merge highly related pair info with full sequence/metadata info:
full<-read.csv("All_Tune_Family_Metadata_And_Aligned_Sequences.csv",header=TRUE)
j<-read.csv("JaHighlyRelatedPairs.csv",header=TRUE)
e<-read.csv("EngHighlyRelatedPairs.csv",header=TRUE)
e<-rbind(e,j)
e$X<-row.names(e)
full<-full[,1:8]
full<-merge(full, e, by.x = "Overall.ID.No.", by.y = "X1", all.x = TRUE)
full<-merge(full, e, by.x = "Overall.ID.No.", by.y = "X2", all.x = TRUE)
full$PairNo<-ifelse(is.na(full$X.x), full$X.y, full$X.x)
full$PairID<-ifelse(is.na(full$X2), full$X1, full$X2)
full$PID<-ifelse(is.na(full$dist.x), full$dist.y, full$dist.x)

write.csv(full,file="FullSequencesWRelatedDataMerged.csv")


#######
#Reorder matrices to put tune families together (took ~2 minutes with 4125x4125 matrix using OLO method, only ~15sec with HC method) (code modified from Hahsler et al., 2008, J. Statistical Software)

##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
eng.res <- dissplot(eng.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(eng.res$x_reordered),"Englishreordered.txt")

##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
ja.res <- dissplot(ja.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(ja.res$x_reordered),"Japanesereordered.txt")


##########Calculate evolutionary rates of highly-related melodic variant pairs
###Japanese (NHK)
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


##########
######To automatically align a given pair of melodies:
#set working drive

#(if not yet installed/loaded, you will need to install/load them first) (copied for ease of training)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("seqinr")

library(seqinr)
library(Biostrings)

#load all 10,064 English and Japanese melodies
melodies<-read.fasta("10064MelodySequences.fasta", seqtype = "AA")

#Set ID number of top and bottom melodies (this example uses ID numbers 9464 and 9465, but these should be replaced with the ID numbers of interest)
top.melody<-9464
bottom.melody<-9465

#save sequences as strings
s3 <- c2s(melodies[[top.melody]])
s4<- c2s(melodies[[bottom.melody]])

(globalAligns3s4<-pairwiseAlignment(s3, s4, type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2))

###########
#Supplementary analyses

#統計学的有意の計算　significance (one pair only)
#load file
#full (no phrase gaps)
bronson<-read.fasta("Bronson4184NoPhraseGaps.fasta", seqtype = "AA")

#save sequences as strings
s3 <- c2s(bronson[[1]])
s4<- c2s(bronson[[2]])

(globalAligns3s4<-pairwiseAlignment(s3, s4, type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2))

# Print out the optimal global alignment and its score

pid(globalAligns3s4,type="PID4")
#find percent identity


globalAligns3s4 <- pairwiseAlignment(s3, s4, type="global", substitutionMatrix = NULL, gapOpening = -.8,gapExtension = -.2, scoreOnly = TRUE)

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
        gapOpening = -.8, gapExtension = -.2, scoreOnly = TRUE)
      randomscores[i] <- score
   }

hist(randomscores, col="red") # Draw a red histogram

sum(randomscores >= globalAligns3s4)
# [1] 0 
#i.e., 0/100 random sequences were greater than observed match score
#i.e., p<0.01


