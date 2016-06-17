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

# define time parameters (in s)
startTime<-0
timeLim<-10000

# specific heat capacity of pyrex
cp<-753

# density of pyrex
rho<-2.21

# thermal conductivity of pyrex
k<-1.005

roomTemp<-20

# record spatial temp profile every no. of seconds:
timeSnapShot<-100

# loss rate
lossRate<-0

alpha<-k/(cp*rho)

# discretize the cylinder - the values don't matter right now
cylMesh<-seq(0,cylLength,stepSize)
cylHeat<-rep(0,length(cylMesh))

# determine the amount of heat energy above the surrounding this has
uStart<-cp*rho*stepSize*pi*(cylWidth/2)^2*(tempStart-roomTemp)

#set the start point
cylHeat[1]<-uStart

thisHeat<-cylHeat

temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+20

tempFrame<-data.frame(cylMesh,temp)

# only allow it to run if it meets the convergence criteria
if ((timeStep/stepSize^2)<=0.5){
  
  # for each time use the explicit method (forward time differential, central second
  # order spatial differential) to numerically solve and iterate the 1D heat eqn
  for (p in seq(startTime,timeLim,timeStep)){
    
    # create vectors holding all of the spatial heat values (current, offset left and offset right by one)
    # vectors instead of loop sidesteps the possibility of non-simultaneous solving, i.e. double-updating points
    
    # use this term to establish the starting condition - left of start is always set Temp
    lastHeat<-c(thisHeat[1],thisHeat[1:length(thisHeat)-1])
    
    # use this term to establish no gradient at the end
    nextHeat<-c(thisHeat[2:length(thisHeat)],thisHeat[length(thisHeat)])
    
    futureHeat<-thisHeat+alpha*(nextHeat-2*thisHeat+lastHeat)*(timeStep/stepSize^2)-lossRate*thisHeat*(timeStep/stepSize^2)
    
    thisHeat<-futureHeat
    
    # effectively "re-heat" the first element back to the start
    thisHeat[1]<-uStart
    
    # use the energy deposited to calc rise above room temp (20)
    temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+roomTemp
    
    # every some number of runs, take an output snapshot
    if (p %% timeSnapShot==0){
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