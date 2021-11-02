# script to install necessary R packages
# These will depend on the R version. Please check sessionInfo.txt

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
                      "BiocManager",
                      "ggthemes",
                      "boot",
                      "dplyr",
                      "ggthemes",
                      "readxl",
                      "brms",
                      "tidybayes",
                      "bayesplot",
                      "projpred",
                      "assertthat",
                      "patchwork")

new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

