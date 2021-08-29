##### import library #####
library(dplyr)
library(GADMTools)


##### load dataset #####
full<-read.csv("https://raw.githubusercontent.com/pesavage/melodic-evolution/master/MelodicEvoSeq.csv",header=TRUE,row.names=1) #Import all 10,000+ sequences
s <- d <- subset(full, PairNo>0)  #Restrict to only highly related pairs

map <- as.data.frame(table(full$NAME_1))
map <- rename(map, NAME_1 = Var1)
attach(map)
map <- map[order(-Freq),]
detach(map)
barplot(map$Freq,las=2,names.arg=map$Var1, cex.names=.7)
write.csv(map,"MapSampleNos.csv")
map <- read.csv("MapSampleNos.csv",row.names=1)

sub.map <- as.data.frame(table(d$NAME_1))
sub.map <- rename(sub.map, NAME_1 = Var1)
attach(sub.map)
sub.map <- sub.map[order(-Freq),]
detach(sub.map)
barplot(sub.map$Freq,las=2,names.arg=map$Var1, cex.names=.7)
write.csv(sub.map,"SubMapSampleNos.csv")
sub.map <- read.csv("SubMapSampleNos.csv",row.names=1)


##### Figure of Japan (Single country) #####

BASEFILE <- "./GADM"
BASEURL <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/"
ja <- gadm_sp.loadCountries("JPN", level = 1, basefile=BASEFILE, baseurl=BASEURL, simplify=0.02)

jards<-full.jards<-readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_JPN_1_sp.rds")) 
jards@data <- merge(jards@data, map, by="NAME_1", all.x=TRUE)
jards@data$Freq[is.na(jards@data$Freq)] <- 0
jadat <- data.frame(NAME_1=jards@data$NAME_1, Freq=jards@data$Freq)

choropleth(ja, jadat, adm.join = "NAME_1",
               value = "Freq",
               breaks = "sd",
               palette="Oranges",
               legend = "Number of melodies",
               title="Number of melodies per region")


##### Figure of US & Canada (Multiple countries) #####
BASEFILE <- "./GADM"
BASEURL <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/"
usacan <- gadm_sp.loadCountries(c("USA", "CAN"), level = 1, basefile=BASEFILE, baseurl=BASEURL, simplify=0.02)

usrds <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_USA_1_sp.rds")) 
usrds@data <- merge(usrds@data, map, by="NAME_1", all.x=TRUE)
usrds@data$Freq[is.na(usrds@data$Freq)] <- 0
usdat <- data.frame(NAME_1=usrds@data$NAME_1, Freq=usrds@data$Freq)

cards <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_CAN_1_sp.rds")) 
cards@data <- merge(cards@data, map, by="NAME_1", all.x=TRUE)
cards@data$Freq[is.na(cards@data$Freq)] <- 0
cadat <- data.frame(NAME_1=cards@data$NAME_1, Freq=cards@data$Freq)

uscadat <- rbind(usdat, cadat)

choropleth(usacan, uscadat, adm.join = "NAME_1",
               value = "Freq",
               breaks = "sd",
               palette="Oranges",
               legend = "Number of melodies",
               title="Number of melodies per region")

##### Figure of UK (Single country) #####
BASEFILE <- "./GADM"
BASEURL <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/"
uk <- gadm_sp.loadCountries("GBR", level = 1, basefile=BASEFILE, baseurl=BASEURL, simplify=0.02)

ukrds<-readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_1_sp.rds")) 
ukrds@data <- merge(ukrds@data, map, by="NAME_1", all.x=TRUE)
ukrds@data$Freq[is.na(ukrds@data$Freq)] <- 0
ukdat <- data.frame(NAME_1=ukrds@data$NAME_1, Freq=ukrds@data$Freq)

choropleth(uk, ukdat, adm.join = "NAME_1",
           value = "Freq",
           breaks = "sd",
           palette="Oranges",
           legend = "Number of melodies",
           title="Number of melodies per region")

##Map highly related melodies only

#Japan
sub.jadat <- merge(jadat, sub.map, by="NAME_1", all.x=TRUE)
sub.jadat$Freq.y[is.na(sub.jadat$Freq.y)] <- 0

choropleth(ja, sub.jadat, adm.join = "NAME_1",
           value = "Freq.y",
           breaks = "sd",
           palette="Oranges",
           legend = "Number of melodies",
           title="Number of melodies per region")

#US + Canada
sub.uscadat <- merge(uscadat, sub.map, by="NAME_1", all.x=TRUE)
sub.uscadat$Freq.y[is.na(sub.uscadat$Freq.y)] <- 0

choropleth(usacan, sub.uscadat, adm.join = "NAME_1",
           value = "Freq.y",
           breaks = "sd",
           palette="Oranges",
           legend = "Number of melodies",
           title="Number of melodies per region")

#UK
sub.ukdat <- merge(ukdat, sub.map, by="NAME_1", all.x=TRUE)
sub.ukdat$Freq.y[is.na(sub.ukdat$Freq.y)] <- 0
choropleth(uk, sub.ukdat, adm.join = "NAME_1",
           value = "Freq.y",
           breaks = "sd",
           palette="Oranges",
           legend = "Number of melodies",
           title="Number of melodies per region")
