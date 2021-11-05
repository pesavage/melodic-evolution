# script to install necessary R packages
# These will depend on the R version. Please check sessionInfo.txt


if (!requireNamespace(pacman, quietly = TRUE)){
  install.packages(pacman)
}
library(pacman)

pacman::p_load(plotrix, 
          seqinr,
          Biostrings,
          seriation,
          tidyr,
          ggplot2,
          dplyr,
          varhandle,
          stringr,
          seqRFLP,
          pwr,
          lsr,
          sp,
          RColorBrewer,
          phangorn,
          GADMTools,
          BiocManager,
          ggthemes,
          boot,
          dplyr,
          ggthemes,
          readxl,
          brms,
          tidybayes,
          bayesplot,
          projpred,
          assertthat,
          # for maps
          patchwork,
          sf,
          rnaturalearth)
