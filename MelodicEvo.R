install.packages(c('plotrix', 'seqinr', 'Biostrings', 'seriation', 'tidyr', 'ggplot2', 'dplyr', 'varhandle', 'stringr', 'seqRFLP', 'pwr', 'lsr', 'sp', 'RColorBrewer', 'phangorn', 'GADMTools'))

#open packages (install first as required - #If Biostrings is not yet installed, enter the following commented out code: 
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

#Files
MelodicEvoAnalysis.R<-"https://raw.githubusercontent.com/pesavage/melodic-evolution/master/MelodicEvoAnalysis.R"
MelodyMap.R<-"https://raw.githubusercontent.com/pesavage/melodic-evolution/master/MelodyMap.R"
PID.R<-"https://raw.githubusercontent.com/pesavage/melodic-evolution/master/PID.R"
Dist.R<-"https://raw.githubusercontent.com/pesavage/melodic-evolution/master/Dist.R"
sensitivity.csv<-"https://raw.githubusercontent.com/pesavage/melodic-evolution/master/sensitivity.csv"

source(PID.R) #uncomment this to recalculate distance matrices of PID among the 10,064 melodic variants, but be aware that this will take up to a month on a standard computer!

#source(Dist.R) # uncomment this to re-dentify highly related pairs of melodies

##########Calculate evolutionary rates of highly-related melodic variant pairs

full<-read.csv("https://raw.githubusercontent.com/pesavage/melodic-evolution/master/MelodicEvoSeq.csv",header=TRUE,row.names=1) #Imports all 10,000+ sequences. The highly related pairs were automatically identified using the scripts in "PID.R" and "Dist.R", but a lot of manual work was required to align related pairs, code functional positions, and count the numbers and sizes of all mutation types to create the "MelodicEvoSeq.csv" file used for subsequent analyses   
s <- d <- subset(full, PairNo>0)  #Restrict to only highly related pairs


#Calculate mutation rates for different functional types

#Full English subset
s <- subset(d, Language=="English")
source(MelodicEvoAnalysis.R)

#########

#Full Japanese subset
s <- subset(d, Language=="Japanese")
source(MelodicEvoAnalysis.R)


#########

#Sensitivity analyses 
#For English sample
e <- subset(d, Language=="English")

#Time
s <- subset(e, Year<median(e$Year,na.rm=TRUE)) #older sample
source(MelodicEvoAnalysis.R)

s <- subset(e, Year>=median(e$Year,na.rm=TRUE)) #newer sample
source(MelodicEvoAnalysis.R)

#Singer
s <- subset(e, Same.singer=="Y") 
source(MelodicEvoAnalysis.R)

s <- subset(e, Same.singer=="N") 
source(MelodicEvoAnalysis.R)

#Coder
s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")
source(MelodicEvoAnalysis.R)

s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")
source(MelodicEvoAnalysis.R)

#For Japanese sample
j <- subset(d, Language=="Japanese")

#Time
s <- subset(j, Year<median(j$Year,na.rm=TRUE)) 
source(MelodicEvoAnalysis.R)

s <- subset(j, Year>=median(j$Year,na.rm=TRUE)) 
source(MelodicEvoAnalysis.R)

#Singer
s <- subset(j, Same.singer=="Y")
source(MelodicEvoAnalysis.R)

s <- subset(j, Same.singer=="N")
source(MelodicEvoAnalysis.R)

#Coder
s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")
source(MelodicEvoAnalysis.R)

s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")
source(MelodicEvoAnalysis.R)

#Descriptive stats for sensitivity analyses
sens<-read.csv(sensitivity.csv)
colMeans(sens,na.rm=TRUE) #Means: Substitution -0.9312500 ; Frequency  -0.7071429 ;    Function; 7.6863636 
std.error(sens[,1],na.rm=TRUE) #Substitution r SE = 0.0163012
std.error(sens[,2],na.rm=TRUE) #Frequency r SE = 0.07389054
std.error(sens[,3],na.rm=TRUE) #Function t SE =  0.9381823


##Map samples
#source(MelodyMap.R) #(Uncomment this once I've finished tweaking map scripts)

# Print R version and packages
sessionInfo()
Sys.time()
