##############################################################################/
##############################################################################/
#R code for analyzing the bioassays of Fusicoccum amygdali populations
##############################################################################/
##############################################################################/

#loading the data
source("load_CydPomEvol_data.R")


##############################################################################/
#Regression analysis of percentage of sporulation for the populations####
##############################################################################/

library(plotrix)
library(Rcpp)
library(magrittr)
library(rlang)
library(tibble)
library(abind)
library(car)
library(gtools)
library(multcomp)
library(sandwich)
library(zoo)
library(mvtnorm)
library(TH.data)
library(scales)
library(drc)


DATA<-read.table("data/DATAX176.txt",header=T,sep="\t",dec=",")


#########Calcul des DL#############

BouTuSpi<-DATA[DATA$pop=="bournes"&DATA$protocole=="tube"&DATA$insecticide=="spinosad",]
BouPlaSpi<-DATA[DATA$pop=="bournes"&DATA$protocole=="plaque"&DATA$insecticide=="spinosad",]
BouTuDel<-DATA[DATA$pop=="bournes"&DATA$protocole=="tube"&DATA$insecticide=="delta",]
BouPlaDel<-DATA[DATA$pop=="bournes"&DATA$protocole=="plaque"&DATA$insecticide=="delta",]
BouTuRyx<-DATA[DATA$pop=="bournes"&DATA$protocole=="tube"&DATA$insecticide=="rynaxypyr",]
BouPlaRyx<-DATA[DATA$pop=="bournes"&DATA$protocole=="plaque"&DATA$insecticide=="rynaxypyr",]
svTuSpi<-DATA[DATA$pop=="sv"&DATA$protocole=="tube"&DATA$insecticide=="spinosad",]
svPlaSpi<-DATA[DATA$pop=="sv"&DATA$protocole=="plaque"&DATA$insecticide=="spinosad",]
svTuDel<-DATA[DATA$pop=="sv"&DATA$protocole=="tube"&DATA$insecticide=="delta",]
svPlaDel<-DATA[DATA$pop=="sv"&DATA$protocole=="plaque"&DATA$insecticide=="delta",]



drm.BouTuSpi<-drm(mort/total~dose,weights=total,data=BouTuSpi,fct=LN.3u(),type="binomial")
drm.BouPlaSpi<-drm(mort/total~dose,weights=total,data=BouPlaSpi,fct=LN.3u(),type="binomial")
drm.BouTuDel<-drm(mort/total~dose,weights=total,data=BouTuDel,fct=LN.3u(),type="binomial")
drm.BouPlaDel<-drm(mort/total~dose,weights=total,data=BouPlaDel,fct=LN.3u(),type="binomial")
drm.BouTuRyx<-drm(mort/total~dose,weights=total,data=BouTuRyx,fct=LN.3u(),type="binomial")
drm.BouPlaRyx<-drm(mort/total~dose,weights=total,data=BouPlaRyx,fct=LN.3u(),type="binomial")

drm.svTuSpi<-drm(mort/total~dose,weights=total,data=svTuSpi,fct=LN.3u(),type="binomial")
drm.svPlaSpi<-drm(mort/total~dose,weights=total,data=svPlaSpi,fct=LN.3u(),type="binomial")
drm.svTuDel<-drm(mort/total~dose,weights=total,data=svTuDel,fct=LN.3u(),type="binomial")
drm.svPlaDel<-drm(mort/total~dose,weights=total,data=svPlaDel,fct=LN.3u(),type="binomial")


LD60_BouTuSpi<-ED(drm.BouTuSpi,60,interval="delta",reference="control")
LD40_BouTuSpi<-ED(drm.BouTuSpi,40,interval="delta",reference="control")
LD50_BouTuSpi<-ED(drm.BouTuSpi,50,interval="delta",reference="control")


LD50_BouPlaSpi<-ED(drm.BouPlaSpi,50,interval="delta",reference="control")
LD50_BouTuSpi<-ED(drm.BouTuSpi,50,interval="delta",reference="control")
LD50_BouTuDel<-ED(drm.BouTuDel,50,interval="delta",reference="control")
LD50_BouPlaDel<-ED(drm.BouPlaDel,50,interval="delta",reference="control")
LD50_BouTuRyx<-ED(drm.BouTuRyx,50,interval="delta",reference="control")
LD50_BouPlaRyx<-ED(drm.BouPlaRyx,50,interval="delta",reference="control")


LD50_svTuSpi<-ED(drm.svTuSpi,50,interval="delta",reference="control")
LD50_svPlaSpi<-ED(drm.svPlaSpi,50,interval="delta",reference="control")
LD50_svTuDel<-ED(drm.svTuDel,50,interval="delta",reference="control")
LD50_svPlaDel<-ED(drm.svPlaDel,50,interval="delta",reference="control")

####Plots#####

plot(drm.BouTuSpi,type="confidence",main="Spinosad Bournes G6",ylim=c(0.02,1),xlim=c(0,100),col="blue", tck = 1)
plot(drm.BouTuSpi,type="obs",col="blue",add=TRUE)
plot(drm.BouPlaSpi,type="confidence",col="red",add=TRUE)
plot(drm.BouPlaSpi,type="obs",col="red",add=TRUE)

text(30,0.3,"DL50 Tube: 9.4",col="blue",add=TRUE)

plot(drm.BouTuDel,type="confidence",main="Deltamethrine Bournes G6",ylim=c(0.02,1),xlim=c(0,100),col="blue", tck = 1)
plot(drm.BouTuDel,type="obs",col="blue",add=TRUE)
plot(drm.BouPlaDel,type="confidence",col="red",add=TRUE)
plot(drm.BouPlaDel,type="obs",col="red",add=TRUE)

text(10,0.4,"DL40 Tube: 0.003",col="blue",add=TRUE)
text(10,0.3,"DL50 Tube: 0.005",col="blue",add=TRUE)
text(10,0.2,"DL50 Plaque: 0.21",col="red",add=TRUE)

plot(drm.BouTuRyx,type="confidence",main="Rynaxypyr Bournes G6",ylim=c(0.02,1),xlim=c(0,100),col="blue", tck = 1)
plot(drm.BouTuRyx,type="obs",col="blue",add=TRUE)
plot(drm.BouPlaRyx,type="confidence",col="red",add=TRUE)
plot(drm.BouPlaRyx,type="obs",col="red",add=TRUE)

text(10,0.4,"DL40 Tube: 0.003",col="blue",add=TRUE)
text(10,0.3,"DL50 Tube: 0.005",col="blue",add=TRUE)
text(10,0.2,"DL50 Plaque: 0.21",col="red",add=TRUE)


plot(drm.svTuSpi,type="confidence",main="Spinosad sv",ylim=c(0.02,1),xlim=c(0,100),col="blue", tck = 1)
plot(drm.svTuSpi,type="obs",col="blue",add=TRUE)
plot(drm.svPlaSpi,type="confidence",col="red",add=TRUE)
plot(drm.svPlaSpi,type="obs",col="red",add=TRUE)

text(10,0.3,"DL50 Tube: 4.6",col="blue",add=TRUE)
text(10,0.2,"DL50 Plaque: 4.4",col="red",add=TRUE)


plot(drm.svTuDel,type="confidence",main="Deltamethrine sv",ylim=c(0.02,1),xlim=c(0,100),col="blue", tck = 1)
plot(drm.svTuDel,type="obs",col="blue",add=TRUE)
plot(drm.svPlaDel,type="confidence",col="red",add=TRUE)
plot(drm.svPlaDel,type="obs",col="red",add=TRUE)

text(10,0.3,"DL50 Tube: 0.008",col="blue",add=TRUE)
text(10,0.2,"DL50 Plaque: 0.05",col="red",add=TRUE)

##############################################################################/
#END
##############################################################################/