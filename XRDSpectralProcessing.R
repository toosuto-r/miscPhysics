library("ggplot2")
library("Peaks")
library.dynam('Peaks', 'Peaks', lib.loc=NULL) 

fName<-"C:/Users/Ryan/Documents/Physics/Nanopillars/Nanopillars XRD/Ga As TRIAL 3.txt"
specData<-read.table(fName,header=FALSE,sep=" ")
names(specData)<-c("twoTheta","raw","proc")

foundPeaks<-SpectrumSearch(specData$proc,threshold=0.5,background = TRUE)

peakPoints<-data.frame(twoTheta=specData[foundPeaks$pos,1],proc=foundPeaks$y[foundPeaks$pos])

peakPlot<-ggplot(data=specData,aes(twoTheta,proc))+geom_line(aes(twoTheta,proc))
peakPlot<-peakPlot+geom_line(aes(specData$twoTheta,foundPeaks$y),color="red")
peakPlot<-peakPlot+geom_point(data=peakPoints,aes(twoTheta,proc))


lam<-1.5418e-10

peakPoints