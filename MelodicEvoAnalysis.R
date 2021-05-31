#### subscript of "MelodicEvo.R"

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Calculate mutation rates for different functional types
MelodicEvoAnalysis = function(s, name){
  mut<-s[,c(1,13:20)]
  mut$nFull<-str_length(mut[,2])
  mut$nOrn<-str_length(mut[,3])
  mut$nFin<-str_length(mut[,4])
  mut$nStr<-str_length(mut[,5])
  mut$nUnStr<-mut$nFull-rowSums(mut[,11:13])
  mut$nOrnMut<-str_length(mut[,6])
  mut$nFinMut<-str_length(mut[,7])
  mut$nStrMut<-str_length(mut[,8])
  mut$nUnStrMut<-str_length(mut[,9])
  
  mut$FinMutRate<-mut$nFinMut/mut$nFin
  mut$StrMutRate<-mut$nStrMut/mut$nStr
  mut$UnStrRate<-mut$nUnStrMut/mut$nUnStr
  mut$OrnMutRate<-mut$nOrnMut/mut$nOrn
  mut$StrongFunctionRate<-(mut$nFinMut+mut$nStrMut)/(mut$nFin+mut$nStr)
  mut$WeakFunctionRate<-(mut$nUnStrMut+mut$nOrnMut)/(mut$nUnStr+mut$nOrn)
  
  sum(mut[,15:18]) #no. of total mutations 
  length(mut$FinMutRate)
  sum(!is.na(mut$OrnMutRate)) #no. of melodies with ornamental notes
  sum(is.na(mut$OrnMutRate)) #no. of melodies without ornamental notes
  
  ####Compare/plot
  
  #mutational distance 
  m <- s[!duplicated(mut$PairNo),] #only using one value per pair, since substitution numbers are identical between variant
  semitone<-colSums(m[,21:31],na.rm=TRUE) 
  
  semitones_df = m[,21:31]
  save(semitones_df, file = "semitones.RData")
  
  #grouped by 2-7 interval size
  interval<-c(sum(semitone[1:2]),sum(semitone[3:4]),sum(semitone[5:6]),sum(semitone[7]),sum(semitone[8:9]),sum(semitone[10:11]))
  print(cor.test(interval,c(2:7),method="spearman",alternative="less"))
  
  jpeg(paste0("figures/NumberSubstitutions_byintervaldistance_", name, ".jpeg"))
  plot(c(2:7),log10(interval),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
  axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
  axis(1, at=2:7, labels=c("2nd","3rd","4th","5th","6th","7th"))
  dev.off()
  
  #grouped by # of semitones
  print(cor.test(semitone,c(1:11),method="spearman",alternative="less"))
  
  jpeg(paste0("figures/NumberSubstitutions_bysemitonedistance_", name, ".jpeg"))
  plot(c(1:11),log10(semitone),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
  axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
  axis(1, at=1:11, labels=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"))
  dev.off()
  
  #Function
  #Strong vs. weak
  data_wide <- mut[ , c(23:24,1)]
  
  #Average rates for each pair
  out <- matrix(NA, nrow=0, ncol=3)
  for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
    rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
    out <- rbind(out,rates)
  }
  data_wide<-as.data.frame(out)
  
  #Check sample sizes
  length(data_wide$StrongFunctionRate) #no. of pairs

  #Make violin plot
  
  #plot
  violin_df = data_wide[,1:2] %>% 
          gather(key="MeasureType", value="Val")
  
  jpeg(paste0("figures/Function_byevolutionaryrate", name, "_1.jpeg"))
  ggplot(violin_df, aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
          geom_violin() +
          stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + 
          geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + 
          ylim(0,0.4) + 
          theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))
  dev.off()
  
  #t tests
  print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE)) #paired t-test
  print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=FALSE)) #unpaired t-test
  
  
  #For all four functional types
  data_wide <- mut[ , c(19:22,1)]
  
  #Average rates for each pair
  out <- matrix(NA, nrow=0, ncol=5)
  for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
    rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
    out <- rbind(out,rates)
  }
  data_wide<-as.data.frame(out)
  
  #Check sample sizes
  length(data_wide$FinMutRate) #no. of pairs
  sum(!is.na(data_wide$OrnMutRate)) #no. of pairs with ornamental notes
  sum(is.na(data_wide$OrnMutRate)) #no. of pairs without ornamental notes
  
  #Make violin plot
  violin_df = data_wide[,1:4] %>% 
          gather(key="MeasureType", value="Val") 
  
  jpeg(paste0("figures/Function_byevolutionaryrate", name, "_2.jpeg"))
  ggplot(violin_df, aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
          geom_violin() + 
          stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + 
          geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + 
          ylim(0,1) + 
          theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))
  dev.off()
  
  #t-tests
  print(t.test(data_wide[,4],data_wide[,3],alternative="greater",paired=TRUE))
  print(t.test(data_wide[,3],data_wide[,2],alternative="greater",paired=TRUE))
  print(t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE))
  colMeans(data_wide,na.rm=TRUE)
  length(subset(data_wide, FinMutRate>0)$FinMutRate) #number of pairs with final mutations
  
  ####Substitution matrix:
  indel<-colSums(m[,32:43],na.rm=TRUE) 
  sum(indel) #total no. of indels
  sub<-colSums(m[,44:109],na.rm=TRUE) 
  sum(sub) #total no. of substitutions
  
  ###The following section creates a substitution matrix in which the diagonal represents the number of times a given note appears unchanged in both melodies from a highly related pair, while the other cells represent the sum of all substitutions involving each possible pair of notes. The diagonal is calculated by subtracting the sum of all substitutions and indels from the total number of times a given note appears.
  n<-"C"
  Cn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"d"
  dn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"D"
  Dn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"e"
  en<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"E"
  En<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"F"
  Fn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"g"
  gn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"G"
  Gn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"a"
  an<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"A"
  An<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"b"
  bn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  n<-"B"
  Bn<-sum(str_count(m[,13], n))-(sum(c(str_count(m[,17], n),str_count(m[,18], n),str_count(m[,19], n),str_count(m[,20], n))))
  
  mat<-cbind(c(0,indel),c(NA,Cn,sub[1:11]),c(rep(NA,2),dn,sub[12:21]),c(rep(NA,3),Dn,sub[22:30]),c(rep(NA,4),en,sub[31:38]),c(rep(NA,5),En,sub[39:45]),c(rep(NA,6),Fn,sub[46:51]),c(rep(NA,7),gn,sub[52:56]),c(rep(NA,8),Gn,sub[57:60]),c(rep(NA,9),an,sub[61:63]),c(rep(NA,10),An,sub[64:65]),c(rep(NA,11),bn,sub[66]),c(rep(NA,12),Bn))
  rownames(mat)<-c("-","C","d","D","e","E","F","g","G","a","A","b","B")
  colnames(mat)<-c("-","C","d","D","e","E","F","g","G","a","A","b","B")
  
  #calculate mutability
  full.mat<-as.matrix(as.dist(as.matrix(mat)))
  sub.mat<-full.mat[2:13,2:13] #exclude substitutions from matrix calculations
  unchanged<-c(mat[2,2],mat[3,3],mat[4,4],mat[5,5],mat[6,6],mat[7,7],mat[8,8],mat[9,9],mat[10,10],mat[11,11],mat[12,12],mat[13,13])
  changed<-colSums(full.mat)[2:13] 
  total<-changed+unchanged
  (mutability<-changed/total)
  
  mat<-rbind(mat,c(NA,mutability))
  write.csv(mat,
            paste0("results/", name, "_SubstitutionMatrix.csv")
            )
  
  #Test correlation between note frequency and mutability
  jpeg(paste0("figures/Frequency_byMutability", name, ".jpeg"))
  plot(log10(total),log10(mutability),pch=16,ylim=c(log10(.1),log10(1)),xlim=c(log10(1),log10(10000)),xaxt="n",yaxt="n",ylab="Mutability",xlab="Note frequency")
  text(log10(total),log10(mutability), names(total), cex=1.5, pos=2, col="red")
  axis(2, at=c(log10(1),log10(.5),log10(.2),log10(.1)), labels=c(1,.5,.2,.1))
  axis(1, at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000)), labels=c(1,10,100,1000,10000))
  dev.off()
  
  print(cor.test(total,mutability,method="spearman",alternative="less"))
  
  
  #####Calculate scale frequencies:
  #extract unique notes/scales
  m$C<-str_count(m[,13], "C")
  m$d<-str_count(m[,13], "d")
  m$D<-str_count(m[,13], "D")
  m$e<-str_count(m[,13], "e")
  m$E<-str_count(m[,13], "E")
  m$F<-str_count(m[,13], "F")
  m$g<-str_count(m[,13], "g")
  m$G<-str_count(m[,13], "G")
  m$a<-str_count(m[,13], "a")
  m$A<-str_count(m[,13], "A")
  m$b<-str_count(m[,13], "b")
  m$B<-str_count(m[,13], "B")
  m$sC<-ifelse(m$C>0,"C","")
  
  for(i in 1:length(m$C)){
    m[i,123]<-ifelse(m[i,111]==0,m[i,122],paste0(c(m[i,122],"d"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,124]<-ifelse(m[i,112]==0,m[i,123],paste0(c(m[i,123],"D"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,125]<-ifelse(m[i,113]==0,m[i,124],paste0(c(m[i,124],"e"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,126]<-ifelse(m[i,114]==0,m[i,125],paste0(c(m[i,125],"E"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,127]<-ifelse(m[i,115]==0,m[i,126],paste0(c(m[i,126],"F"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,128]<-ifelse(m[i,116]==0,m[i,127],paste0(c(m[i,127],"g"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,129]<-ifelse(m[i,117]==0,m[i,128],paste0(c(m[i,128],"G"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,130]<-ifelse(m[i,118]==0,m[i,129],paste0(c(m[i,129],"a"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,131]<-ifelse(m[i,119]==0,m[i,130],paste0(c(m[i,130],"A"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,132]<-ifelse(m[i,120]==0,m[i,131],paste0(c(m[i,131],"b"),collapse = ""))
  }
  for(i in 1:length(m$C)){
    m[i,133]<-ifelse(m[i,121]==0,m[i,132],paste0(c(m[i,132],"B"),collapse = ""))
  }
  names(m)[133] <- "scale"
  m$scaleNum<-str_length(m$scale)
  
  jpeg(paste0("figures/Frequency_ofScale", name, ".jpeg"))
  barplot(table(m$scaleNum))
  dev.off()
  
  #barplot of scales ordered by frequency
  scale<-as.data.frame(
    sort(
      table(m$scale), decreasing = TRUE
      )
    )
  jpeg(paste0("figures/Frequency_ofScale_sorted", name, ".jpeg"))
  barplot(scale$Freq,las=2,names.arg=scale$Var1, cex.names=.7, main = name)
  dev.off()
}
