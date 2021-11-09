#### subscript of "MelodicEvo.R"
source('./helper.R')

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  single_song <- mean(x)
  ymin <- single_song-std.error(x)*1.96
  ymax <- single_song+std.error(x)*1.96
  return(c(y=single_song,ymin=ymin,ymax=ymax))
}

# Calculate mutation rates for different functional types
MelodicEvoAnalysis = function(s, name){
  require(stringr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggthemes)
  require(boot)
  
  s$n_notes      = str_length(s$Full.note.sequence..unaligned.)
  s$n_ornamental = str_length(s$Ornamental.notes)
  s$n_final      = str_length(s$Final.note)
  s$n_stressed   = str_length(s$Stressed.notes)
  s$n_unstressed = s$n_notes - 
    rowSums(s[,c("n_ornamental", "n_final", "n_stressed")], na.rm = TRUE) 
  
  s$n_ornamentalmutations = str_length(s$Ornamental.mutations)
  s$n_finalmutations      = str_length(s$Final.mutations)
  s$n_stressedmutations   = str_length(s$Stress.mutations)
  s$n_unstressedmutations = str_length(s$Unstressed.mutations)
  
  # If no mutations, change to 0 instead of NA
  s$n_ornamentalmutations = ifelse(is.na(s$n_ornamentalmutations), 0,
                                   s$n_ornamentalmutations)
  s$n_finalmutations      = ifelse(is.na(s$n_finalmutations), 0,
                                   s$n_finalmutations)
  s$n_stressedmutations   = ifelse(is.na(s$n_stressedmutations), 0,
                                   s$n_stressedmutations)
  s$n_unstressedmutations = ifelse(is.na(s$n_unstressedmutations), 0,
                                   s$n_unstressedmutations)
  
  s$finalmutation_rate      = s$n_finalmutations / s$n_final
  s$stressedmutation_rate   = s$n_stressedmutations / s$n_stressed
  s$unstressedmutation_rate = s$n_unstressedmutations / s$n_unstressed
  s$ornamentalmutation_rate = s$n_ornamentalmutations / s$n_ornamental
  
  s$strongfunction_rate = (s$n_finalmutations + s$n_stressedmutations) / (s$n_final + s$n_stressed)
  s$weakfunction_rate   = (s$n_unstressedmutations + s$n_ornamentalmutations) /
    rowSums(s[,c("n_unstressed", "n_ornamental")], na.rm = TRUE)
  
  single_song = s[!duplicated(s$PairNo),]
  
  # Get columns of semitonal distances
  semitonal_columns = str_detect(colnames(single_song), 
                                 "\\.semitone\\.substitutions")
  # Calculate semitonal distance frequency
  semitone = colSums(single_song[,semitonal_columns],na.rm=TRUE) 
  
  # Calculate semitonal bootstrapped 95% CI intervals
  sumFunc = function(x,i){sum(x[i], na.rm = TRUE)}
  bootSum = apply(single_song[,semitonal_columns], 2, 
                    function(x) {
                    boot(x, sumFunc, R = 1000)
                    })
  semitonalci_df = sapply(bootSum, function(b) boot.ci(b, type = "norm")$normal)
  semitonalci_df = data.frame(t(semitonalci_df))
  colnames(semitonalci_df) = c("interval", "low", "high")
  
  semitonalci_df$semitone = 1:11
  semitonalci_df$high = semitonalci_df$high + 1
  
  # grouped by 2-7 interval size
  int_df = single_song[,semitonal_columns]
  interval_df = data.frame(second = rowSums(int_df[,1:2], na.rm = TRUE),
                           third = rowSums(int_df[,3:4], na.rm = TRUE),
                           fourth = rowSums(int_df[,5:6], na.rm = TRUE),
                           fifth = int_df[,7],
                           sixth = rowSums(int_df[,8:9], na.rm = TRUE),
                           seventh = rowSums(int_df[,10:11], na.rm = TRUE))
                           
  interval_substitutions = colSums(interval_df, na.rm = TRUE)
  
  
  ## Interval bootstrapped 95% CI
  bootSum = apply(interval_df, 2, 
                  function(x) boot(x, sumFunc, R = 1000))
  intervalci_df = sapply(bootSum, function(b) boot.ci(b, type = "norm")$normal)
  intervalci_df = data.frame(t(intervalci_df))
  colnames(intervalci_df) = c("ci_interval", "low", "high")
  intervalci_df$interval = 2:7
  
  print(cor.test(interval_substitutions, 
                 c(2:7), 
                 method="spearman",
                 alternative="less", 
                 exact = FALSE))
  
  print(cor.test(interval_substitutions, 
                 c(2:7), 
                 method="pearson",
                 alternative="less", 
                 exact = FALSE))
  
  # Graph of Intervals by Substitution count
  jpeg(paste0("figures/NumberSubstitutions_byintervaldistance_", name, ".jpeg"))
  plot(c(2:7),
       interval_substitutions,
       pch = 16, xaxt = "n", 
       ylab = "Number of substitutions",
       xlab = "Substitution distance (intervals)",
       ylim = c(0, max(intervalci_df$high)))
  arrows(x0 = intervalci_df$interval, y0 = intervalci_df$low, 
         x1 = intervalci_df$interval, y1 = intervalci_df$high, 
         length = 0.05, 
         angle = 90, 
         code = 3)
  axis(1, at = 2:7, labels = c("2nd", "3rd", "4th", "5th", "6th", "7th"))
  dev.off()
  
  #grouped by # of semitones
  print(cor.test(semitone,
                 c(1:11),
                 method="spearman",
                 alternative="less", 
                 exact = FALSE))
  
  print(cor.test(semitone,
                 c(1:11),
                 method="pearson",
                 alternative="less", 
                 exact = FALSE))
  
  jpeg(paste0("figures/NumberSubstitutions_bysemitonedistance_", name, ".jpeg"))
  plot(c(1:11), 
       semitone,
       pch = 16, xaxt = "n",
       ylab = "Number of substitutions",
       xlab = "Substitution distance (semitones)",
       ylim = c(0, max(semitonalci_df$high)))
  arrows(x0 = semitonalci_df$semitone, y0 = semitonalci_df$low, 
         x1 = semitonalci_df$semitone, y1 = semitonalci_df$high, 
         length = 0.05, 
         angle = 90, 
         code = 3)
  axis(1, at = 1:11, labels = c("1", "2", "3", "4",
                                "5", "6", "7", "8",
                                "9", "10", "11"))
  dev.off()

  #### Function ####
  #Strong vs. weak
  strong_weak = aggregate(s[, c("PairNo", "strongfunction_rate", 
                              "weakfunction_rate","PID")], 
                   by = list(PairNo = s$PairNo),
                   FUN = mean)

  #Check sample sizes
  length(strong_weak$strongfunction_rate) #no. of pairs
  
  #Make violin plot
  strong_weak<-strong_weak[,c("PairNo","strongfunction_rate", "weakfunction_rate","PID")]
  #plot
  violin_df = strong_weak %>% 
    dplyr::select(strongfunction_rate, weakfunction_rate) %>% 
    gather(key="MeasureType", value="Val") %>% 
    mutate(Val = (1 - Val) * 100)
  
  
  p1 = ggplot(violin_df, 
         aes(x = reorder(MeasureType, Val), 
             y = Val, 
             fill = MeasureType)) +
    geom_violin()  + 
    geom_jitter(size = 1,
                position = position_jitter(0.3)) + 
    stat_summary(fun.data = data_summary, 
                 geom = "pointrange", 
                 color = "red", 
                 size = .6) +
    #geom_hline(yintercept = 85, colour = "red", linetype = "dashed") + 
    theme_base() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none") +
    xlab("") + 
    ylab("Stability (% identity)") + 
    scale_x_discrete(labels = c("strongfunction_rate" = "Strong Function", 
                                "weakfunction_rate" = "Weak Function")) + 
    scale_y_continuous(labels=function(x) paste0(x,"%")) 
  
  
  ggsave(filename = 
           paste0("figures/Function_byevolutionaryrate", name, "_1.jpeg"), plot = p1)
  
  #### t-tests ####
  ## paired t-test
  print(t.test(strong_weak$weakfunction_rate,
               strong_weak$strongfunction_rate,
               alternative = "greater",
               paired = TRUE,
               na.action = "na.omit"))
  print(cohensD( x = strong_weak$weakfunction_rate, y = strong_weak$strongfunction_rate, method = "paired"))
  
  # unpaired t-test
  print(t.test(strong_weak$weakfunction_rate,
               strong_weak$strongfunction_rate,
               alternative = "greater",
               paired=FALSE)) 
  
 #Plot difference between strong vs. weak function as function of PID
  jpeg(paste0("figures/FunctionalDifferenceVsPID_", name, ".jpeg"))
  plot(strong_weak$PID,
       ((1 - strong_weak$strongfunction_rate) * 100)- ((1 - strong_weak$weakfunction_rate) * 100),
       pch = 16, 
       ylab = "Strong function % identity - weaker function % identity",
       xlab = "Overall % identity")
  
  dev.off()
  
  print(cor.test(strong_weak$PID,
                 ((1 - strong_weak$strongfunction_rate) * 100)- ((1 - strong_weak$weakfunction_rate) * 100)),method="spearman")
  
  print(cor.test(strong_weak$PID,
                 ((1 - strong_weak$strongfunction_rate) * 100)- ((1 - strong_weak$weakfunction_rate) * 100)),method="pearson")
  
   # For all four functional types
  functional_types = aggregate(s[, c("finalmutation_rate", 
                              "stressedmutation_rate",
                              "unstressedmutation_rate", 
                              "ornamentalmutation_rate",
                              "PairNo")], 
                        by = list(PairNo = s$PairNo),
                        FUN = mean, na.rm = TRUE)
  
  #Check sample sizes
  #no. of pairs
  length(functional_types$finalmutation_rate)
  #no. of pairs with ornamental notes
  sum(!is.na(functional_types$ornamentalmutation_rate))  
  #no. of pairs without ornamental notes
  sum(is.na(functional_types$ornamentalmutation_rate)) 
  functional_types<-functional_types[,c("finalmutation_rate",                             
                                        "stressedmutation_rate",
                                        "unstressedmutation_rate", 
                                        "ornamentalmutation_rate",
                                        "PairNo")]
  #Make violin plot
  violin_df = functional_types %>% 
    select(finalmutation_rate, stressedmutation_rate,
           unstressedmutation_rate, ornamentalmutation_rate) %>% 
    gather(key="MeasureType", value="Val") %>% 
    mutate(Val = (1 - Val) * 100)
  
  violin_df$MeasureType = factor(violin_df$MeasureType, 
                                 levels = c("ornamentalmutation_rate",
                                            "unstressedmutation_rate",
                                            "stressedmutation_rate",
                                            "finalmutation_rate"))
  
  p2 = ggplot(violin_df, 
              aes(x = MeasureType, 
                  y = Val, 
                  fill = MeasureType)) +
    geom_violin()  + 
    geom_jitter(size = 1,
                position = position_jitter(0.3)) + 
    stat_summary(fun.data = data_summary, 
                 geom = "pointrange", 
                 color = "red",
                 size = .6) +
    #geom_hline(yintercept = 85, colour = "red", linetype = "dashed") + 
    theme_base() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none") +
    xlab("") + 
    ylab("Stability (% melodic identity)") + 
    scale_x_discrete(labels = c("unstressedmutation_rate" = "Unstressed",
                                "stressedmutation_rate" = "Stressed", 
                                "finalmutation_rate" = "Final",
                                "ornamentalmutation_rate" = "Ornamental")) +
    scale_y_continuous(labels=function(x) paste0(x,"%"))
  
  ggsave(filename = 
           paste0("figures/Function_byevolutionaryrate", name, "_2.jpeg"), plot = p2)
  
  ### Functional t-tests ####
  print(t.test(functional_types$ornamentalmutation_rate,
               functional_types$unstressedmutation_rate,
               alternative = "greater", 
               paired = TRUE,
               na.action = "na.omit"))
  
  print(t.test(functional_types$unstressedmutation_rate,
               functional_types$stressedmutation_rate,
               alternative = "greater",
               paired = TRUE,
               na.action = "na.omit"))
  
  print(t.test(functional_types$stressedmutation_rate,
               functional_types$finalmutation_rate,
               alternative = "greater", 
               paired = TRUE,
               na.action = "na.omit"))
  
  colMeans(functional_types, na.rm = TRUE)
  #number of pairs with final mutations
  length(
    subset(functional_types, finalmutation_rate > 0)$finalmutation_rate
    ) 
  
  #### Substitution matrix ####
  # Get all columns representing a deletion in the form of a letter 
  # followed by a .
  indel_columns = str_detect(colnames(single_song), "^[A-Za-z]{1}\\.$")
  indel = colSums(single_song[,indel_columns],na.rm=TRUE) 
  sum(indel) # total no. of indels
  
  # Get all substitution columns, in the form of two letters
  sub_columns = str_detect(colnames(single_song), "^[A-Za-z]{2}$")
  sub = colSums(single_song[,sub_columns],na.rm=TRUE) 
  sum(sub) # total no. of substitutions
  
  ###The following section creates a substitution matrix in which the diagonal represents the number of times a given note appears unchanged in both melodies from a highly related pair, while the other cells represent the sum of all substitutions involving each possible pair of notes. The diagonal is calculated by subtracting the sum of all substitutions and indels from the total number of times a given note appears. This information is used to manually create Fig. 3C.
  notes = c("C", "d",	"D",	"e",	"E",	"F",	"g",
            "G",	"a",	"A",	"b",	"B")
  
  song_counts = 
    sapply(notes, function(n) {
      total = 
        sum(str_count(single_song$Full.note.sequence..unaligned., n), 
                  na.rm = TRUE)
      ornamental_mutations = 
        sum(str_count(single_song$Ornamental.mutations, n), 
                                 na.rm = TRUE)
      final_mutations = 
        sum(str_count(single_song$Final.mutations, n), 
                            na.rm = TRUE)
      stress_mutations = 
        sum(str_count(single_song$Stress.mutations, n), 
                             na.rm = TRUE)
      unstressed_mutations = sum(str_count(single_song$Unstressed.mutations, 
                                           n), na.rm = TRUE)
      # Total - all mutations gives the note count
      total - (ornamental_mutations + final_mutations + 
                 stress_mutations + unstressed_mutations)
    })
  
  mat = matrix(NA, 
               ncol = length(notes), 
               nrow = length(notes),
               dimnames = list(notes, notes))
  
  mat[lower.tri(mat)] = sub
  diag(mat) = song_counts
  # Add indels
  mat_indel = rbind(rep(NA, ncol(mat)), mat)
  mat_indel = cbind(c(0, indel), mat_indel)
  mat_indel[upper.tri(mat_indel)] = t(mat_indel)[upper.tri(mat_indel)]
  
  #### Calculate Mutability #### [these analyses were not included in the final version but are included in earlier PsyArXiv preprints]
  mat[upper.tri(mat)] = t(mat)[upper.tri(mat)]

  print(colSums(mat_indel[,2:ncol(mat_indel)], na.rm = TRUE))
  
  changed = colSums(mat_indel[,2:ncol(mat_indel)], na.rm = TRUE) - song_counts
  
  print(colSums(mat_indel[,2:ncol(mat_indel)]))
  
  total = changed + song_counts
  print("CHANGED")
  print(changed)
  print("MUTABILITY")
  (mutability = changed / total)
  
  write.csv(mutability, paste0("results/", name, "_mutability.csv"))
  write.csv(total, paste0("results/", name, "_notecounts.csv"))
  
  mat = rbind(mat, c(mutability))
  write.csv(mat,
            paste0("results/", name, "_SubstitutionMatrix.csv")
  )
  
  #Test correlation between note frequency and mutability
  jpeg(paste0("figures/Frequency_byMutability", name, ".jpeg"))
  plot(log10(total),
       log10(mutability),
       pch = 16,
       ylim = c(log10(.1), log10(1)),
       xlim = c(log10(1), log10(10000)), 
       xaxt = "n", yaxt = "n", 
       ylab = "Mutability",
       xlab = "Note frequency")
  text(log10(total), 
       log10(mutability), 
       names(total), 
       cex = 1.5, 
       pos = 2, 
       col = "red")
  axis(2, 
       at = c(log10(1), log10(.5), log10(.2), log10(.1)), 
       labels = c(1, .5, .2, .1))
  axis(1, 
       at = c(log10(1), log10(10), log10(100), log10(1000), log10(10000)), 
       labels = c(1, 10, 100, 1000, 10000))
  dev.off()
  
  print(cor.test(total,
                 mutability,
                 method="spearman",
                 alternative="less", 
                 exact = FALSE))
  
  #print descriptive stats for melody length
  print(mean(s$n_notes))
  print(median(s$n_notes))
  print(min(s$n_notes))
  print(max(s$n_notes))
  
  #####Calculate scale frequencies:
  # extract unique notes/scales
  # note occurrence in each song
  note_occurences = sapply(notes, 
         function(note) 
           str_count(single_song$Full.note.sequence..unaligned., note))
  colnames(note_occurences) = notes
  
  # single_song$sC<-ifelse(single_song$C>0,"C","")
  
  # for(i in 1:length(single_song$C)){
  #   single_song[i,123]<-ifelse(single_song[i,111]==0,single_song[i,122],paste0(c(single_song[i,122],"d"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,124]<-ifelse(single_song[i,112]==0,single_song[i,123],paste0(c(single_song[i,123],"D"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,125]<-ifelse(single_song[i,113]==0,single_song[i,124],paste0(c(single_song[i,124],"e"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,126]<-ifelse(single_song[i,114]==0,single_song[i,125],paste0(c(single_song[i,125],"E"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,127]<-ifelse(single_song[i,115]==0,single_song[i,126],paste0(c(single_song[i,126],"F"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,128]<-ifelse(single_song[i,116]==0,single_song[i,127],paste0(c(single_song[i,127],"g"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,129]<-ifelse(single_song[i,117]==0,single_song[i,128],paste0(c(single_song[i,128],"G"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,130]<-ifelse(single_song[i,118]==0,single_song[i,129],paste0(c(single_song[i,129],"a"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,131]<-ifelse(single_song[i,119]==0,single_song[i,130],paste0(c(single_song[i,130],"A"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,132]<-ifelse(single_song[i,120]==0,single_song[i,131],paste0(c(single_song[i,131],"b"),collapse = ""))
  # }
  # for(i in 1:length(single_song$C)){
  #   single_song[i,133]<-ifelse(single_song[i,121]==0,single_song[i,132],paste0(c(single_song[i,132],"B"),collapse = ""))
  # }
  # names(single_song)[133] <- "scale"
  # single_song$scaleNum<-str_length(single_song$scale)
  
  # jpeg(paste0("figures/Frequency_ofScale", name, ".jpeg"))
  # barplot(table(single_song$scaleNum))
  # dev.off()
  # 
  # # barplot of scales ordered by frequency
  # scale = as.data.frame(
  #   sort(
  #     table(single_song$scale), decreasing = TRUE
  #   )
  # )
  # jpeg(paste0("figures/Frequency_ofScale_sorted", name, ".jpeg"))
  # barplot(scale$Freq,las=2,names.arg=scale$Var1, cex.names=.7, main = name)
  # dev.off()
  
  # Outputs for testing
  invisible(
    list(mut = s, interval = interval_substitutions, semitone = semitone, 
       strongweak = strong_weak,
       functional_types = functional_types, song_counts = song_counts,
       total = total, note_occurences = note_occurences)
    )
}
