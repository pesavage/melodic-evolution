list.of.packages <- c("plotrix", 
                      "seqinr",
                      "Biostrings",
                      "seriation",
                      "tidyr",
                      "ggplot2",
                      "dplyr",
                      "varhandle",
                      "stringr",
                      "seqRFLP",
                      "pwr",
                      "lsr",
                      "sp",
                      "RColorBrewer",
                      "phangorn",
                      "GADMTools",
                      "BiocManager")

new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

