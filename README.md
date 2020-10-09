# melodic-evolution
Code and source data for Savage, Chiba, Currie, Suzuki, &amp; Atkinson "Sequence alignment of folk song melodies reveals cross-cultural mechanisms of musical evolution"

To run the analyses, run the master file "MelodicEvo.R" in R or R Studio. Several of the sub-scripts are commented out because they have been pre-calculated. In particular, the distance matrices and identification of highly related melodies have been pre-calculated, so you don't need to run the "PID.R" or "Dist.R" subscripts. To recalculate distance matrices and re-identify highly related pairs, uncomment the following line of code from "MelodicEvo.R" after downloading the full English and Japanese distance matrices from https://osf.io/nhvzw/ (but be aware that this will take up to a month on a standard computer):

source("PID.R")
source("Dist.R")
