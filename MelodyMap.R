##Map samples (subscript of "MelodicEvo.R")

###Map samples

#barplot by state/prefecture
map<-as.data.frame(table(full$NAME_1))
map <- rename(map, NAME_1 = Var1)
attach(map)
map <- map[order(-Freq),]
detach(map)
barplot(map$Freq,las=2,names.arg=map$Var1, cex.names=.7)
write.csv(map,"CountySampleNos.csv")

map<-read.csv("CountySampleNos.csv",row.names=1)

#Japan
ja<- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_JPN_1_sp.rds")) 
#ja$NAME_1 <- as.factor(iconv(as.character(ja$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
ja@data<- merge(ja@data, map,by="NAME_1",all.x=TRUE)
ja@data$Freq[is.na(ja@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(ja, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(ja@data$NAME_1)), col = "#081D58", main = "Melody sample size (Japan)")

#USA
us <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_USA_1_sp.rds")) 
#us$NAME_1 <- as.factor(iconv(as.character(us$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
us@data<- merge(us@data, map,by="NAME_1",all.x=TRUE)
us@data$Freq[is.na(us@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(us, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(us@data$NAME_1)), col = "#081D58", main = "Melody sample size")

#Canada
can <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_CAN_1_sp.rds")) 
can$NAME_1 <- as.factor(iconv(as.character(can$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
can@data<- merge(can@data, map,by="NAME_1",all.x=TRUE)
can@data$Freq[is.na(can@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(can, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(us@data$NAME_1)), col = "#081D58", main = "Melody sample size")

#UK
uk <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_1_sp.rds")) 
#uk$NAME_1 <- as.factor(iconv(as.character(uk$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
uk@data<- merge(uk@data, map,by="NAME_1",all.x=TRUE)
uk@data$Freq[is.na(uk@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(uk, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(uk@data$NAME_1)), col = "#081D58", main = "Melody sample size")


###Need to combine UK, Ireland, & Isle of Man into one "British Isles" map

#Ireland
irl <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_IRL_0_sp.rds")) 
irl$NAME_1 <- irl$NAME_0
#irl$NAME_1 <- as.factor(iconv(as.character(irl$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
irl@data<- merge(irl@data, map,by="NAME_1",all.x=TRUE)
irl@data$Freq[is.na(irl@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(irl, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(uk@data$NAME_1)), col = "#081D58", main = "Melody sample size")

#Isle of Man
  imn <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_IMN_0_sp.rds")) 
imn$NAME_1 <- imn$NAME_0
imn$NAME_1 <- as.factor(iconv(as.character(imn$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
imn@data<- merge(imn@data, map,by="NAME_1",all.x=TRUE)
imn@data$Freq[is.na(imn@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(imn, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(uk@data$NAME_1)), col = "#081D58", main = "Melody sample size")



#####Following is some code I've attempted to join multiple countries into a single plot, so far without complete success:
## load a file from GADM (you just have to specify the countries "special part" of the file name, like "ARG" for Argentina. Optionally you can specify which level you want to have
Brit = gadm_sp_loadCountries(c("GBR","IRL"), level=1, basefile="./")
NAm = gadm_sp_loadCountries(c("USA","CAN"), level=1, basefile="./")

NAm$spdf@data<- merge(NAm$spdf@data, map,by.x="NAME_1",by.y="Var1",all.x=TRUE)
NAm$spdf@data$Freq[is.na(NAm$spdf@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(NAm, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(NAm$spdf@data$NAME_1)), col = "#081D58", main = "Melody sample size (N. America)")

DAT<- merge(NAm$spdf@data,map,by.x="NAME_1",by.y="Var1",all.x=TRUE)

