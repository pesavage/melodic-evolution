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

#########

#Full Japanese subset
s <- subset(d, Language=="Japanese")
source("MelodicEvoAnalysis.R")


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