suppressPackageStartupMessages({
  ## Tests
  library(testthat) 
  # For analysis
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
  library(readxl)
})

source("analysis/MelodicEvoFunction.R")

full = suppressMessages(
  read_xlsx("data/MelodicEvoSeq.xlsx", .name_repair = "universal")
)
d    = subset(full, PairNo>0)  # Restrict to only highly related pairs

english = subset(d, Language=="English")

e = MelodicEvoAnalysis(english, "english")

# Useful objects
notes = c("C", "d",	"D",	"e",	"E",	"F",	"g",
          "G",	"a",	"A",	"b",	"B")

# Check English counts match what is expected 
testthat::test_that("Mutation counts & Rates check", {
  expect_equal(sum(e$mut[,c("n_ornamentalmutations", "n_finalmutations", 
                            "n_stressedmutations", "n_unstressedmutations")]), 
               2385)
  expect_equal(sum(!is.na(e$mut$ornamentalmutation_rate)), 55)
  expect_equal(sum(is.na(e$mut$ornamentalmutation_rate)), 429)
})

# Intervals and Semitones check
testthat::test_that("Interval & Semitones check", {
  expect_equal(unname(e$interval), c(326, 202, 73, 22, 13, 2))
  expect_equal(unname(e$semitone), c(70, 256, 135, 67, 73, 0, 22, 6 ,7 ,2, 0))
})

# Note counts check
test_that("Note Counts check", {
  expect_equal(names(e$song_counts), notes)
  expect_equal(unname(e$song_counts), c(2231, 0, 1235, 325, 1077, 1000, 1, 
                                        2157, 4, 794, 430, 192))
})

# Note totals check
test_that("Note Counts check", {
  expect_equal(names(e$total), notes)
  expect_equal(unname(e$total), c(2670, 1, 1501, 410, 1368, 1285, 2, 
                                        2610, 14, 1014, 586, 262))
})
