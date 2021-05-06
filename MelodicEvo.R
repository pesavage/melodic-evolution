source("install.R")

#If Biostrings is already installed you can comment out the following three lines of code: 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 BiocManager::install("Biostrings") #If not yet installed, follow installation instructions 

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
library(seqRFLP)
library(pwr)
library(lsr)
library(sp)
library(RColorBrewer)
library(phangorn)
library(GADMTools)

# To calculate distance matrices  
#source(PID.R) #uncomment this to recalculate distance matrices of PID among the 10,064 melodic variants, but be aware that this will take up to a month on a standard computer!
#source(Dist.R) # uncomment this to re-dentify highly related pairs of melodies

##### Calculate evolutionary rates of highly-related melodic variant pairs ####
 # Import all 10,000+ sequences. 
 # The highly related pairs were automatically identified using the scripts in "PID.R" and "Dist.R", but a lot of manual work was required to align related pairs, code functional positions, and count the numbers and sizes of all mutation types to create the "MelodicEvoSeq.csv" file used for subsequent analyses
full<-read.csv("MelodicEvoSeq.csv", header=TRUE, row.names=1)
d <- subset(full, PairNo>0)  #Restrict to only highly related pairs

## Source analysis function
source("MelodicEvoAnalysis.R")

##### Calculate mutation rates for different functional types ####
# Full English subset
english <- subset(d, Language=="English")
MelodicEvoAnalysis(english)

#########

# Full Japanese subset
japanese <- subset(d, Language=="Japanese")
MelodicEvoAnalysis(japanese)

#########

#### Sensitivity analyses ####
# For English sample
e <- subset(d, Language=="English")

#Time
english_old <- subset(e, Year<median(e$Year,na.rm=TRUE)) #older sample
MelodicEvoAnalysis(english_old)

english_new <- subset(e, Year>=median(e$Year,na.rm=TRUE)) #newer sample
MelodicEvoAnalysis(english_new)

#Singer
singer <- subset(e, Same.singer=="Y") 
MelodicEvoAnalysis(singer)

not_singer <- subset(e, Same.singer=="N") 
MelodicEvoAnalysis(not_singer)

#Coder
coder_pes <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")
MelodicEvoAnalysis(coder_pes)

coder_gc <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")
MelodicEvoAnalysis(coder_gc)

#For Japanese sample
j <- subset(d, Language=="Japanese")

#Time
japanese_old <- subset(j, Year<median(j$Year,na.rm=TRUE)) 
MelodicEvoAnalysis(japanese_old)

japanese_new <- subset(j, Year>=median(j$Year,na.rm=TRUE)) 
MelodicEvoAnalysis(japanese_new)

#Singer
singer_japanese <- subset(j, Same.singer=="Y")
MelodicEvoAnalysis(singer_japanese)

notsinger_japanese <- subset(j, Same.singer=="N")
MelodicEvoAnalysis(notsinger_japanese)

#Coder
japancoder_pes <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")
MelodicEvoAnalysis(japancoder_pes)

japancoder_gc <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")
MelodicEvoAnalysis(japancoder_gc)

#Descriptive stats for sensitivity analyses
sens<-read.csv("sensitivity.csv")
colMeans(sens,na.rm=TRUE) #Means: Substitution -0.9312500 ; Frequency  -0.7071429 ;    Function; 7.6863636 
std.error(sens[,1],na.rm=TRUE) #Substitution r SE = 0.0163012
std.error(sens[,2],na.rm=TRUE) #Frequency r SE = 0.07389054
std.error(sens[,3],na.rm=TRUE) #Function t SE =  0.9381823


##Map samples
#source(MelodyMap.R) #(Uncomment this once I've finished tweaking map scripts)

# Write R version and packages
# If you are reproducing these results, then please check you have the same 
# Package versions as in the sessionInfo.txt file
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")