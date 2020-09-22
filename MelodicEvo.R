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

#source("PID.R") #uncomment this to recalculate distance matrices of PID among the 10,064 melodic variants, but be aware that this will take up to a month on a standard computer!

#read back in distance matrices calculated by "PID.R" subscript and identify highly related (>=85%) melodic variants:

#English folk songs
eng.d<-as.dist(as.matrix(read.table("EnglishDist.txt"))) 

#plot histogram
hist(eng.d)
abline(v=0.15, lty=2)

#sort by similarity
eng.m<-as.matrix(eng.d)
xy <- t(combn(colnames(eng.m), 2))
d<-data.frame(xy, dist=eng.m[xy])
d<-unfactor(d) #fixes ID numbers being coded as factors rather than numbers
r<-subset(d,dist<=.15)

attach(r)
r.ord <- r[order(dist),]
detach(r)
length(r.ord$dist) #774 highly related pairs
length(unique(c(r.ord[,1],r.ord[,2])))#651 unique melodic variants

#filter to remove identical pairs
r.ord<-subset(r.ord,dist>0)
length(r.ord$dist) #reduces to 760 non-identical pairs)

#filter to only use more related pairs in cases of melodies highly related to multiple other melodies
r.ord<-r.ord %>% distinct(X1, .keep_all = TRUE)
length(r.ord$dist) #reduces to 383 highly related pairs not with same X1 melody

r.ord<-r.ord %>% distinct(X2, .keep_all = TRUE)
length(r.ord$dist) #reduces to 316 highly related pairs also not with same X2 melody

#Remove remaining duplicates across two columns
dup<-duplicated(c(r.ord[,1],r.ord[,2]))[(length(r.ord$dist)+1):(length(r.ord$dist)*2)]
r.ord<-r.ord[!dup, ]
length(r.ord$dist) #reduces to 242 highly related unique [non-shared] pairs)


write.csv(r.ord,file="EngHighlyRelatedPairs.csv")



##########

#Japanese folk songs
ja.d<-as.dist(as.matrix(read.table("JapanDist.txt"))) 

#plot histogram
hist(ja.d)
abline(v=0.15, lty=2)

#sort by similarity
ja.m<-as.matrix(ja.d)
xy <- t(combn(colnames(ja.m), 2))
d<-data.frame(xy, dist=ja.m[xy])
d<-unfactor(d) #fixes ID numbers being coded as factors rather than numbers
r<-subset(d,dist<=.15)

attach(r)
r.ord <- r[order(dist),]
detach(r)
length(r.ord$dist) #143 highly related pairs
length(unique(c(r.ord[,1],r.ord[,2])))#220 unique melodic variants

#filter to remove identical pairs
r.ord<-subset(r.ord,dist>0)
length(r.ord$dist) #reduces to 131 non-identical pairs)

#filter to only use more related pairs in cases of melodies highly related to multiple other melodies
r.ord<-r.ord %>% distinct(X1, .keep_all = TRUE)
length(r.ord$dist) #reduces to 110 highly related pairs not with same X1 melody

r.ord<-r.ord %>% distinct(X2, .keep_all = TRUE)
length(r.ord$dist) #reduces to 101 highly related pairs not with same X1 melody

#Remove remaining duplicates across two columns
dup<-duplicated(c(r.ord[,1],r.ord[,2]))[(length(r.ord$dist)+1):(length(r.ord$dist)*2)]
r.ord<-r.ord[!dup, ]
length(r.ord$dist) #reduces to 91 highly related unique [non-shared] pairs)


write.csv(r.ord,file="JaHighlyRelatedPairs.csv")

#Merge highly related pair info with full sequence/metadata info:
full<-read.csv("All_Tune_Family_Metadata_And_Aligned_Sequences.csv",header=TRUE)
j<-read.csv("JaHighlyRelatedPairs.csv",header=TRUE)
e<-read.csv("EngHighlyRelatedPairs.csv",header=TRUE)
e<-rbind(e,j)
e$X<-row.names(e)
full<-full[,1:8]
full<-merge(full, e, by.x = "Overall.ID.No.", by.y = "X1", all.x = TRUE)
full<-merge(full, e, by.x = "Overall.ID.No.", by.y = "X2", all.x = TRUE)
full$PairNo<-ifelse(is.na(full$X.x), full$X.y, full$X.x)
full$PairID<-ifelse(is.na(full$X2), full$X1, full$X2)
full$PID<-ifelse(is.na(full$dist.x), full$dist.y, full$dist.x)

write.csv(full,file="FullSequencesWRelatedDataMerged.csv")


#######
#Reorder matrices to put tune families together (took ~2 minutes with 4125x4125 matrix using OLO method, only ~15sec with HC method) (code modified from Hahsler et al., 2008, J. Statistical Software)

##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
eng.res <- dissplot(eng.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(eng.res$x_reordered),"Englishreordered.txt")

##plot reordered matrix (using hierarchical clustering with optimal leaf reordering)
ja.res <- dissplot(ja.d, method="OLO", options = list(main = "Dissimilarity plot with seriation", col= colorRampPalette(c("yellow", "red"))( 20 )))

#write re-ordered distance matrix
write.table(as.matrix(ja.res$x_reordered),"Japanesereordered.txt")


##########Calculate evolutionary rates of highly-related melodic variant pairs

#New automated code (initially for full sample)
full<-read.csv("MelodicEvoSeq.csv",header=TRUE,row.names=1) #Import all 10,000+ sequences
d<-s <- subset(full, PairNo>0)  #Restrict to only highly related pairs
s<-d

#Calculate mutation rates for different functional types
mut<-s[,1:10]
mut$nFull<-str_length(mut[,2])
mut$nOrn<-str_length(mut[,3])
mut$nFin<-str_length(mut[,4])
mut$nStr<-str_length(mut[,5])
mut$nUnStr<-mut$nFull-rowSums(mut[,12:14])
mut$nOrnMut<-str_length(mut[,6])
mut$nFinMut<-str_length(mut[,7])
mut$nStrMut<-str_length(mut[,8])
mut$nUnStrMut<-str_length(mut[,9])
mut[is.na(mut)] <- 0
mut$FinMutRate<-mut$nFinMut/mut$nFin
mut$StrMutRate<-mut$nStrMut/mut$nStr
mut$UnStrRate<-mut$nUnStrMut/mut$nUnStr
mut$OrnMutRate<-mut$nOrnMut/mut$nOrn
mut$StrongFunctionRate<-(mut$nFinMut+mut$nStrMut)/(mut$nFin+mut$nStr)
mut$WeakFunctionRate<-(mut$nUnStrMut+mut$nOrnMut)/(mut$nUnStr+mut$nOrn)

sum(mut[,16:19]) #3462 (total number of mutations)
length(mut$FinMutRate)
sum(!is.na(mut$OrnMutRate)) #149 melodies with ornamental notes
sum(is.na(mut$OrnMutRate)) #507 melodies without ornamental notes

#Power analysis
#Hypothesis 1:
pwr.t.test(n=86, sig.level=.0125, power = 0.8, type = "paired", alternative = "greater") #d = 0.3374611 (for 86 Japanese pairs)
pwr.t.test(n=242, sig.level=.0125, power = 0.8, type = "paired", alternative = "greater") #d = 0.199224 (for 242 English pairs)
#Hypothesis 2:
pwr.r.test(n=6, sig.level=.0125, power = 0.8, alternative = "greater") #r = 0.9374896 for 2nd-7th grouping of mutational distance
pwr.r.test(n=11, sig.level=.0125, power = 0.8, alternative = "greater") #r = 0.7869453 for 0-11 raw semitone analysis of mutational distance

#Estimate effect sizes based on previous analyses:
#English
est<-read.csv("AllMutationRates.csv")
est$StrongFunctionRate<-rowMeans(est[,1:2])
est$WeakFunctionRate<-rowMeans(est[,3:4],na.rm=TRUE)
cohensD(est$StrongFunctionRate,est$WeakFunctionRate) # d = 0.6932685
#Japanese:
est<-read.csv("AllMutationRatesNHK.csv")
est$StrongFunctionRate<-rowMeans(est[,1:2])
est$WeakFunctionRate<-rowMeans(est[,3:4],na.rm=TRUE)
cohensD(est$StrongFunctionRate,est$WeakFunctionRate) # d = 1.31745

####Compare/plot

#mutational distance (grouped by # of semitones)
m <- s[!duplicated(mut$PairNo),] #only using one value per pair, since substitution numbers are identical between variant
semitone<-colSums(m[,18:28],na.rm=TRUE) #Note these columns are numbered differently from the Japanese columns
cor.test(semitone,c(1:11),method="spearman",alternative="less")

barplot(semitone,names.arg=c("1(m2nd)","2(M2nd)","3(m3rd)","4(M3rd)","5(P5th)","6(A4/D5)", "7(P5th)","8(m6th)", "9(M6th)","10(m7th)", "11(M7th)"),ylab=expression(paste("Number of substitutions")))

#mutational distance (grouped by 2-7 interval size)
interval<-c(sum(semitone[1:2]),sum(semitone[3:4]),sum(semitone[5:6]),sum(semitone[7]),sum(semitone[8:9]),sum(semitone[10:11]))
cor.test(interval,c(2:7),method="spearman",alternative="less")
#x <- barplot(interval,names.arg=c("2nd","3rd","4th","5th","6th","7th"),ylab=expression(paste("Number of substitutions")))
x <- plot(c(2:7),log10(interval),ylim=c(0,3),pch=16,xaxt="n",yaxt="n",ylab="Number of substitutions (log scale)",xlab="Substitution distance")
axis(2, at=c(0,1,2,3), labels=c(1,10,100,1000))
axis(1, at=2:7, labels=c("2nd","3rd","4th","5th","6th","7th"))

#For strong vs. weak function only
data_wide <- mut[ , c(24:25,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=3)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$StrongFunctionRate) #328 pairs

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#Make violin plot

data_wide[,1:2] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,0.4) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t tests
t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE)
colSums(!is.na(data_wide[1:2]))




#For all four functional types
data_wide <- mut[ , c(20:23,1)]

#Average rates for each pair
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(mut[!duplicated(mut$PairNo),]$PairNo)){
  rates<-colMeans(data_wide[(i*2-1):(i*2),],na.rm=TRUE)
  out <- rbind(out,rates)
}
data_wide<-as.data.frame(out)

#Check sample sizes
length(data_wide$FinMutRate) #328 pairs
sum(!is.na(data_wide$OrnMutRate)) #104 pairs with ornamental notes
sum(is.na(data_wide$OrnMutRate)) #224 pairs without ornamental notes

#Define function for mean and 95% confidence interval
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-std.error(x)*1.96
  ymax <- m+std.error(x)*1.96
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#Make violin plot

data_wide[,1:4] %>% 
  gather(key="MeasureType", value="Val") %>%
  ggplot( aes(x=reorder(MeasureType, Val), y=Val, fill=MeasureType)) +
  geom_violin() +stat_summary(fun.data=data_summary, geom="pointrange",color="red",width=1,size=.6) + geom_jitter(binaxis='y', stackdir='center', size=1,position=position_jitter(0.3)) + ylim(0,1) + theme(axis.text=element_text(size=21),axis.title=element_text(size=23,face="bold"))

#t-tests
t.test(data_wide[,4],data_wide[,3],alternative="greater",paired=TRUE)
t.test(data_wide[,3],data_wide[,2],alternative="greater",paired=TRUE)
t.test(data_wide[,2],data_wide[,1],alternative="greater",paired=TRUE)
colMeans(data_wide,na.rm=TRUE)
length(subset(data_wide, FinMutRate>0)$FinMutRate) #number of pairs with final mutations


####Exploratory analysis of substitution matrix:
indel<-colSums(m[,29:40],na.rm=TRUE) 
sum(indel) #1798 indels total
sub<-colSums(m[,41:106],na.rm=TRUE) 
sum(sub) #769 substitutions total

###Transition matrix
n<-"C"
Cn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"d"
dn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"D"
Dn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"e"
en<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"E"
En<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"F"
Fn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"g"
gn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"G"
Gn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"a"
an<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"A"
An<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"b"
bn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))
n<-"B"
Bn<-sum(str_count(m[,2], n))-(sum(c(str_count(m[,6], n),str_count(m[,7], n),str_count(m[,8], n),str_count(m[,9], n))))

mat<-cbind(c(0,indel),c(NA,Cn,sub[1:11]),c(rep(NA,2),dn,sub[12:21]),c(rep(NA,3),Dn,sub[22:30]),c(rep(NA,4),en,sub[31:38]),c(rep(NA,5),En,sub[39:45]),c(rep(NA,6),Fn,sub[46:51]),c(rep(NA,7),gn,sub[52:56]),c(rep(NA,8),Gn,sub[57:60]),c(rep(NA,9),an,sub[61:63]),c(rep(NA,10),An,sub[64:65]),c(rep(NA,11),bn,sub[66]),c(rep(NA,12),Bn))
rownames(mat)<-c("-","C","d","D","e","E","F","g","G","a","A","b","B")
colnames(mat)<-c("-","C","d","D","e","E","F","g","G","a","A","b","B")
plot(1:12,c(mat[2,2],mat[3,3],mat[4,4],mat[5,5],mat[6,6],mat[7,7],mat[8,8],mat[9,9],mat[10,10],mat[11,11],mat[12,12],mat[13,13]),type="line")
write.csv(mat,"SubstitutionMatrix.csv") #Rename after running English and Japanese subsets


full.mat<-as.matrix(as.dist(as.matrix(mat)))
sub.mat<-full.mat[2:13,2:13] #exclude substitutions from matrix calculations
unchanged<-c(mat[2,2],mat[3,3],mat[4,4],mat[5,5],mat[6,6],mat[7,7],mat[8,8],mat[9,9],mat[10,10],mat[11,11],mat[12,12],mat[13,13])
changed<-colSums(sub.mat) 
total<-changed+unchanged
mutability<-changed/total
trans.mat<-sub.mat*mutability/total #transition matrxi
for(i in 1:12){
  trans.mat[i,i]<-1-mutability[i]
}
write.csv(signif(trans.mat,digits=3),"TransitionMatrix.csv") #Rename after running English and Japanese subsets

#Test correlation between note frequency and mutability
#plot(log(total),mutability,ylim=c(0,max(mutability,na.rm=TRUE)),xlim=c(0,max(log(total),na.rm=TRUE)))
plot(log10(total),log10(mutability),pch=16,ylim=c(log10(.05),log10(1)),xlim=c(log10(1),log10(10000)),xaxt="n",yaxt="n",ylab="Mutability",xlab="Note frequency")
text(log10(total),log10(mutability), names(total), cex=1.5, pos=2, col="red")
axis(2, at=c(log10(1),log10(.5),log10(.2),log10(.1),log10(.05)), labels=c(1,.5,.2,.1,.05))
axis(1, at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000)), labels=c(1,10,100,1000,10000))
cor.test(total,mutability,method="spearman",alternative="less")

#####Calculate scale frequencies:
#extract unique notes/scales
#unique(strsplit(as.character(m[1,2]), "")[[1]])
#m[1,143:148]<-unique(strsplit(as.character(m[1,2]), "")[[1]])
m$C<-str_count(m[,2], "C")
m$d<-str_count(m[,2], "d")
m$D<-str_count(m[,2], "D")
m$e<-str_count(m[,2], "e")
m$E<-str_count(m[,2], "E")
m$F<-str_count(m[,2], "F")
m$g<-str_count(m[,2], "g")
m$G<-str_count(m[,2], "G")
m$a<-str_count(m[,2], "a")
m$A<-str_count(m[,2], "A")
m$b<-str_count(m[,2], "b")
m$B<-str_count(m[,2], "B")
m$sC<-ifelse(m$C>0,"C","")

for(i in 1:length(m$C)){
  m[i,143]<-ifelse(m[i,131]==0,m[i,142],paste0(c(m[i,142],"d"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,144]<-ifelse(m[i,132]==0,m[i,143],paste0(c(m[i,143],"D"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,145]<-ifelse(m[i,133]==0,m[i,144],paste0(c(m[i,144],"e"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,146]<-ifelse(m[i,134]==0,m[i,145],paste0(c(m[i,145],"E"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,147]<-ifelse(m[i,135]==0,m[i,146],paste0(c(m[i,146],"F"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,148]<-ifelse(m[i,136]==0,m[i,147],paste0(c(m[i,147],"g"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,149]<-ifelse(m[i,137]==0,m[i,148],paste0(c(m[i,148],"G"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,150]<-ifelse(m[i,138]==0,m[i,149],paste0(c(m[i,149],"a"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,151]<-ifelse(m[i,139]==0,m[i,150],paste0(c(m[i,150],"A"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,152]<-ifelse(m[i,140]==0,m[i,151],paste0(c(m[i,151],"b"),collapse = ""))
}
for(i in 1:length(m$C)){
  m[i,153]<-ifelse(m[i,141]==0,m[i,152],paste0(c(m[i,152],"B"),collapse = ""))
}
names(m)[153] <- "scale"
m$scaleNum<-str_length(m$scale)
barplot(table(m$scaleNum)) #hist(m$scaleNum)

#barplot of scales ordered by frequency
map<-as.data.frame(table(m$scale))
attach(map)
map <- map[order(-Freq),]
detach(map)
barplot(map$Freq,las=2,names.arg=map$Var1, cex.names=.7)


table(m$scale,m$Language)
#double barplot (draft code, not working, not sure if needed)
#scale<-as.matrix(table(m$scale,m$Language))
#for(i in 1:length(scale[,1])){
#  scale[i,3]<-ifelse(scale[i,2]==0,"red","blue")
#}

###Subset analyses (repeat above from "mut<-s[,1:10]", changing definition of sample as follows):
#Full English subset
s <- subset(d, Language=="English")
#Full Japanese subset
s <- subset(d, Language=="Japanese")

#Sensitivity analyses (repeat above from "mut<-s[,1:10]", changing definition of sample as follows, ignoring inferential statistics):
#For English sample
e <- subset(d, Language=="English")

#Time
s <- subset(e, Year<median(e$Year,na.rm=TRUE)) 
s <- subset(e, Year>=median(e$Year,na.rm=TRUE)) 
#Singer
s <- subset(e, Same.singer...Y.N...=="Y") #21 pairs 
s <- subset(e, Same.singer...Y.N...=="N") #160 pairs
s <- subset(e, Same.singer...Y.N...=="?") #58 pairs
#Coder
s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")  #122 pairs
s <- subset(e, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")  #120 pairs


#For Japanese sample
j <- subset(d, Language=="Japanese")

#Time
s <- subset(j, Year<median(j$Year,na.rm=TRUE)) 
s <- subset(j, Year>=median(j$Year,na.rm=TRUE)) 
#Singer
s <- subset(j, Same.singer...Y.N...=="Y") #23 pairs
s <- subset(j, Same.singer...Y.N...=="N") #21 pairs
s <- subset(j, Same.singer...Y.N...=="?") #39 pairs
#Coder
s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="PES")  #46 pairs
s <- subset(j, Alignment.functional.coding.performed.by..PES...Patrick.E..Savage..GC...Gakuto.Chiba.=="GC")  #40 pairs





###Map sample

full<-read.csv("FullMelodies.csv",header=TRUE)

#barplot by state/prefecture
map<-as.data.frame(table(full$Country..County.State.))
attach(map)
map <- map[order(-Freq),]
detach(map)
barplot(map$Freq,las=2,names.arg=map$Var1, cex.names=.7)

#barplot by Child ballad no.
e.full <- subset(full, Language=="English")
map<-as.data.frame(table(e.full$Child.Ballad.no..NHK.Volume.no.))
attach(map)
map <- map[order(-Freq),]
detach(map)
barplot(map$Freq,las=2,names.arg=map$Var1, cex.names=.7)
