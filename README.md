# melodic-evolution
Code and source data for Savage, Currie, &amp; Atkinson "Automatic sequence alignment reveals cross-cultural mechanisms of musical evolution"

To run the analyses, download the two large distance matrices from https://osf.io/nhvzw/, then open "MelodicEvo.R" and copy and paste the relevant sections into R or R studio (modifying working directory as necessary). The distance matrices have been pre-calculated using the "PID.R" subscript - to recalculate distance matrices, uncomment the following line of code from "!MelodicEvo.R" (but be aware that this will take up to a month on a standard computer!):

#source("PID.R")
