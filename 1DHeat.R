library(ggplot2)
library(RColorBrewer)

ptm<-proc.time()


#define the width of the cylinder (for heat capacity) and the length
cylWidth<-0.1
cylLength<-2

# define a starting temperature for the first element (degrees C)
tempStart<-600

# define the x step size
stepSize<-0.1
timeStep<-0.0001

# specific heat capacity of pyrex
cp<-753

# density of pyrex
rho<-2.21

# thermal conductivity of pyrex
k<-1.005

alpha<-k/(cp*rho)

# discretize the cylinder - the values don't matter right now
cylMesh<-seq(0,cylLength,stepSize)
cylHeat<-rep(0,length(cylMesh))

# determine the amount of heat energy above the surrounding this has
uStart<-cp*rho*stepSize*pi*(cylWidth/2)^2*(tempStart-20)

#set the start point
cylHeat[1]<-uStart

thisHeat<-cylHeat

startTime<-0
timeLim<-10000

temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+20

tempFrame<-data.frame(cylMesh,temp)

if ((timeStep/stepSize^2)<=0.5){
  
  for (p in seq(startTime,timeLim,timeStep)){
    
    lastHeat<-c(thisHeat[1],thisHeat[1:length(thisHeat)-1])
    nextHeat<-c(thisHeat[2:length(thisHeat)],thisHeat[length(thisHeat)])
    
    nextHeat<-thisHeat+alpha*(nextHeat-2*thisHeat+lastHeat)*(timeStep/stepSize^2)
    #heat2[1]<-uStart
    
    thisHeat<-nextHeat
    thisHeat[1]<-uStart
    
    temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+20
    
    
    if (p %% 100==0){
      tempFrame<-cbind(tempFrame,temp)
    }
    
  }
} else {
  cat("The square of the spatial step should be no more than twice the time step (currently ",(timeStep/stepSize^2),")\n",sep="")
  cat("Recommended maximum spatial step for this time step: ",((timeStep/0.5)^0.5),"\n", sep="")
  cat("Recommended minimum time step for this spatial step: ",(0.5*stepSize^2),"\n", sep="")
}
temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+20

tempMelt<-melt(tempFrame,id="cylMesh")


# define a colour palette to be interpolated - use yellow/orange/red and reverse later
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))
pall <- cols(length(unique(tempMelt$variable)))
dev.new()

# plot the data frame up, using the previously mentioned colour palette
timePlot<-ggplot(tempMelt,aes(cylMesh,y=value,color=variable))+geom_line()
timePlot+scale_colour_manual(values=rev(pall))+labs(x="x (m)", y=expression(paste("Temperature (",degree,"C) ")))+theme(legend.position="none")

print(proc.time()-ptm)