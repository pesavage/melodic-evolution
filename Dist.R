#subscript of MelodicEvo.R that identifies highly related pairs of melodies

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
length(r.ord$dist) #141 highly related pairs
length(unique(c(r.ord[,1],r.ord[,2])))#216 unique melodic variants

#filter to remove identical pairs
r.ord<-subset(r.ord,dist>0)
length(r.ord$dist) #reduces to 129 non-identical pairs)

#filter to only use more related pairs in cases of melodies highly related to multiple other melodies
r.ord<-r.ord %>% distinct(X1, .keep_all = TRUE)
length(r.ord$dist) #reduces to 108 highly related pairs not with same X1 melody

r.ord<-r.ord %>% distinct(X2, .keep_all = TRUE)
length(r.ord$dist) #reduces to 99 highly related pairs not with same X1 melody

#Remove remaining duplicates across two columns
dup<-duplicated(c(r.ord[,1],r.ord[,2]))[(length(r.ord$dist)+1):(length(r.ord$dist)*2)]
r.ord<-r.ord[!dup, ]
length(r.ord$dist) #reduces to 89 highly related unique [non-shared] pairs (later reduced to 86 after correcting errors in 3 pairs)


write.csv(r.ord,file="JaHighlyRelatedPairs.csv")

#Merge highly related pair info with full sequence/metadata info:
full<-read.csv("MelodicEvoSeq.csv",header=TRUE) #Import all sequences
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
full$PID<-(1-full$PID)*100

write.csv(full,file="FullSequencesWRelatedDataMerged.csv")

##Used FullSequencesWRelatedDataMerged.csv as basis for manually aligning sequences. After eliminating several highly related pairs that were included due to data entry errrs in original sequence entry, saved the full file as "MelodicEvoSeq.csv"


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