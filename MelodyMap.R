##Map samples (subscript of "MelodicEvo.R")
#source("MelodyMap.R")

#Japan
ja<- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_JPN_1_sp.rds")) 
#ja$NAME_1 <- as.factor(iconv(as.character(ja$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
ja@data<- merge(ja@data, map,by.x="NAME_1",by.y="Var1",all.x=TRUE)
ja@data$Freq[is.na(ja@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(ja, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(ja@data$NAME_1)), col = "#081D58", main = "Melody sample size (Japan)")


#UK
uk <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_1_sp.rds")) 
#uk$NAME_1 <- as.factor(iconv(as.character(uk$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
uk@data<- merge(uk@data, map,by.x="NAME_1",by.y="Var1",all.x=TRUE)
uk@data$Freq[is.na(uk@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(uk, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(uk@data$NAME_1)), col = "#081D58", main = "Melody sample size (Britain)")

#Ireland
irl <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_IRL_0_sp.rds")) 
irl$NAME_1 <- irl$NAME_0
#irl$NAME_1 <- as.factor(iconv(as.character(irl$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
gbr<-
  
  #Isle of Man (NB: here we only need NAME_0,not NAME_1)
  imn <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_IMN_1_sp.rds")) 
imn$NAME_1 <- imn$NAME_0
imn$NAME_1 <- as.factor(iconv(as.character(imn$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.



#USA
us <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_USA_1_sp.rds")) 
#us$NAME_1 <- as.factor(iconv(as.character(us$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
us@data<- merge(us@data, map,by.x="NAME_1",by.y="Var1",all.x=TRUE)
us@data$Freq[is.na(us@data$Freq)] <- 0 #(Add 0 for states without samples)
spplot(us, "Freq", col.regions = colorRampPalette(brewer.pal(9, "YlGnBu"))(length(us@data$NAME_1)), col = "#081D58", main = "Melody sample size (US)")

#Canada
can <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_CAN_1_sp.rds")) 
can$NAME_1 <- as.factor(iconv(as.character(can$NAME_1), "latin1", "UTF-8"))  # Convert the encoding to UTF-8 in order to avoid the problems with 'tildes' and 'eñes'.
