setwd("/Users/pesavage/Documents/Research/Papers/Unpublished/Savage et al. (in prep) English and Japanese folk song evolution automated analysis/melodic-evolution-master")

#open packages (install first as required - #If Biostrings is not yet installed, enter the following commented out code: 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("Biostrings") #If not yet installed, follow installation instructions 

#open packages
library(plotrix)
library(seqinr)
library(Biostrings)
library(seriation)
library(tidyr)
library(ggplot2)
library(dplyr)
library(varhandle)
library(stringr)
library("seqRFLP")
library(pwr)
library(lsr)
library(sp)
library(RColorBrewer)
library(phangorn)

#source("PID.R") #uncomment this to recalculate distance matrices of PID among the 10,064 melodic variants, but be aware that this will take up to a month on a standard computer!

#source("Dist.R") # uncomment this to re-dentify highly related pairs of melodies

##########Calculate evolutionary rates of highly-related melodic variant pairs

full<-read.csv("MelodicEvoSeq.csv",header=TRUE,row.names=1) #Import all 10,000+ sequences
d<-s <- subset(full, PairNo>0)  #Restrict to only highly related pairs
s<-d

#Calculate mutation rates for different functional types

#Full English subset
s <- subset(d, Language=="English")
source("MelodicEvoAnalysis.R")

####FOLLOWING CODE IS COPIED HERE TEMPORARILY BECAUSE IT DOESN'T PRINT WHEN RUN VIA THE PREVOUS "SOURCE" COMMAND:
#Calculate mutation rates for different functional types

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
mut[is.na(mut)] <- 0
mut$FinMutRate<-mut$nFinMut/mut$nFin
mut$StrMutRate<-mut$nStrMut/mut$nStr
mut$UnStrRate<-mut$nUnStrMut/mut$nUnStr
mut$OrnMutRate<-mut$nOrnMut/mut$nOrn
mut$StrongFunctionRate<-(mut$nFinMut+mut$nStrMut)/(mut$nFin+mut$nStr)
mut$WeakFunctionRate<-(mut$nUnStrMut+mut$nOrnMut)/(mut$nUnStr+mut$nOrn)

sum(mut[,15:18]) #3462 mutations total
length(mut$FinMutRate)
sum(!is.na(mut$OrnMutRate)) #149 melodies with ornamental notes
sum(is.na(mut$OrnMutRate)) #507 melodies without ornamental notes

####Compare/plot

#mutational distance 
m <- s[!duplicated(mut$PairNo),] #only using one value per pair, since substitution numbers are identical between variant
semitone<-colSums(m[,21:31],na.rm=TRUE) 

#grouped by 2-7 interval size
interval<-c(sum(semitone[1:2]),sum(semitone[3:4]),sum(semitone[5:6]),sum(semitone[7]),sum(semitone[8:9]),sum(semitone[10:11]))
print(cor.test(interval,c(2:7),method="spearman",alternative="less"))
x <- plot(c(2:7),log10(interval),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
axis(1, at=2:7, labels=c("2nd","3rd","4th","5th","6th","7th"))

#grouped by # of semitones
print(cor.test(semitone,c(1:11),method="spearman",alternative="less"))
x <- plot(c(1:11),log10(semitone),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
axis(1, at=1:11, labels=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"))


#Function
#Strong vs. weak
data_wide <- mut[ , c(23:24,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=3)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$StrongFunctionRate) #328 pairs


#Make violin plot

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#plot
data_wide[,1:2] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,0.4) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t tests
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE)) #paired t-test
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=FALSE)) #unpaired t-test


#For all four functional types
data_wide <- mut[ , c(19:22,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$FinMutRate) #328 pairs
sum(!is.na(data_wide$OrnMutRate)) #104 pairs with ornamental notes
sum(is.na(data_wide$OrnMutRate)) #224 pairs without ornamental notes


#Make violin plot

data_wide[,1:4] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,1) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t-tests
print(t.test(data_wide[,4],data_wide[,3],alternative="greater",paired=TRUE))
print(t.test(data_wide[,3],data_wide[,2],alternative="greater",paired=TRUE))
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE))
colMeans(data_wide,na.rm=TRUE)
length(subset(data_wide, FinMutRate>0)$FinMutRate) #number of pairs with final mutations

#########

#Full Japanese subset
s <- subset(d, Language=="Japanese")
source("MelodicEvoAnalysis.R")

#Calculate mutation rates for different functional types

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
mut[is.na(mut)] <- 0
mut$FinMutRate<-mut$nFinMut/mut$nFin
mut$StrMutRate<-mut$nStrMut/mut$nStr
mut$UnStrRate<-mut$nUnStrMut/mut$nUnStr
mut$OrnMutRate<-mut$nOrnMut/mut$nOrn
mut$StrongFunctionRate<-(mut$nFinMut+mut$nStrMut)/(mut$nFin+mut$nStr)
mut$WeakFunctionRate<-(mut$nUnStrMut+mut$nOrnMut)/(mut$nUnStr+mut$nOrn)

sum(mut[,15:18]) #3462 mutations total
length(mut$FinMutRate)
sum(!is.na(mut$OrnMutRate)) #149 melodies with ornamental notes
sum(is.na(mut$OrnMutRate)) #507 melodies without ornamental notes

####Compare/plot

#mutational distance 
m <- s[!duplicated(mut$PairNo),] #only using one value per pair, since substitution numbers are identical between variant
semitone<-colSums(m[,21:31],na.rm=TRUE) 

#grouped by 2-7 interval size
interval<-c(sum(semitone[1:2]),sum(semitone[3:4]),sum(semitone[5:6]),sum(semitone[7]),sum(semitone[8:9]),sum(semitone[10:11]))
print(cor.test(interval,c(2:7),method="spearman",alternative="less"))
x <- plot(c(2:7),log10(interval),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
axis(1, at=2:7, labels=c("2nd","3rd","4th","5th","6th","7th"))

#grouped by # of semitones
print(cor.test(semitone,c(1:11),method="spearman",alternative="less"))
x <- plot(c(1:11),log10(semitone),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
axis(1, at=1:11, labels=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"))


#Function
#Strong vs. weak
data_wide <- mut[ , c(23:24,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=3)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$StrongFunctionRate) #328 pairs


#Make violin plot

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#plot
data_wide[,1:2] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,0.4) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t tests
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE)) #paired t-test
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=FALSE)) #unpaired t-test


#For all four functional types
data_wide <- mut[ , c(19:22,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$FinMutRate) #328 pairs
sum(!is.na(data_wide$OrnMutRate)) #104 pairs with ornamental notes
sum(is.na(data_wide$OrnMutRate)) #224 pairs without ornamental notes


#Make violin plot

data_wide[,1:4] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,1) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t-tests
print(t.test(data_wide[,4],data_wide[,3],alternative="greater",paired=TRUE))
print(t.test(data_wide[,3],data_wide[,2],alternative="greater",paired=TRUE))
print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE))
colMeans(data_wide,na.rm=TRUE)
length(subset(data_wide, FinMutRate>0)$FinMutRate) #number of pairs with final mutations


#########

#Sensitivity analyses 
#For English sample
e <- subset(d, Language=="English")

#Time
s <- subset(e, Year<median(e$Year,na.rm=TRUE)) #older sample
source("MelodicEvoAnalysis.R")

s <- subset(e, Year>=median(e$Year,na.rm=TRUE)) #newer sample
source("MelodicEvoAnalysis.R")

#Singer
s <- subset(e, Same.singer=="Y") #21 pairs 
source("MelodicEvoAnalysis.R")

s <- subset(e, Same.singer=="N") #160 pairs
source("MelodicEvoAnalysis.R")

#Coder
s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")  #122 pairs
source("MelodicEvoAnalysis.R")

s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")  #120 pairs
source("MelodicEvoAnalysis.R")

#For Japanese sample
j <- subset(d, Language=="Japanese")

#Time
s <- subset(j, Year<median(j$Year,na.rm=TRUE)) 
source("MelodicEvoAnalysis.R")

s <- subset(j, Year>=median(j$Year,na.rm=TRUE)) 
source("MelodicEvoAnalysis.R")

#Singer
s <- subset(j, Same.singer=="Y") #23 pairs
source("MelodicEvoAnalysis.R")

s <- subset(j, Same.singer=="N") #21 pairs
source("MelodicEvoAnalysis.R")

#Coder
s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")  #46 pairs
source("MelodicEvoAnalysis.R")

s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")  #40 pairs
source("MelodicEvoAnalysis.R")

#Descriptive stats for sensitivity analyses
sens<-read.csv("sensitivity.csv")
colMeans(sens,na.rm=TRUE) #Means: Substitution -0.9343750 ; Frequency  -0.7442857 ;    Function; 7.6227273 
std.error(sens[,1],na.rm=TRUE) #Substitution r SE = 0.01655719
std.error(sens[,2],na.rm=TRUE) #Frequency r SE = 0.05159771
std.error(sens[,3],na.rm=TRUE) #Function t SE =  0.9414673


##Map samples
#source("MelodyMap.R") #(Uncomment this once I've finished tweaking map scripts)