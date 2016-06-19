# check if the required packages are there, and install if they aren't

if(!require("ggplot2",character.only = TRUE)){
  install.packages("ggplot2",dep=TRUE)
  if(!require("ggplot2",character.only = TRUE)) stop("Package not found")
}
if(!require("RColorBrewer",character.only = TRUE)){
  install.packages("RColorBrewer",dep=TRUE)
  if(!require("RColorBrewer",character.only = TRUE)) stop("Package not found")
}
if(!require("reshape2",character.only = TRUE)){
  install.packages("reshape2",dep=TRUE)
  if(!require("reshape2",character.only = TRUE)) stop("Package not found")
}



ptm<-proc.time()


#define the width of the cylinder (for heat capacity) and the length (in cm)
cylWidth<-20
cylLength<-200

# define a starting temperature for the first element (degrees C)
tempStart<-600

# define the x step size (m)
stepSize<-0.001
timeStep<-0.0000001

#define the x step size (cm)
stepSize<-5
timeStep<-10

# define time parameters (in s)
startTime<-0
timeLim<-60*60*24

# specific heat capacity of pyrex (K/kg/degC)
cp<-753
# and in J/g/degC
cp<-cp/1000

# density of pyrex (g/cm^3)
rho<-2.21

#density in g/mm^3
rho<-rho/(10^3)

# thermal conductivity of pyrex (W/m/K)
k<-1.14
# and in W/mm/K
k<-k/1000

roomTemp<-20

# record spatial temp profile every no. of seconds:
timeSnapShot<-60

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
  
  p<-startTime  
  # for each time use the explicit method (forward time differential, central second
  # order spatial differential) to numerically solve and iterate the 1D heat eqn
  while (p<timeLim){
    
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
    
    rm(lastHeat)
    rm(nextHeat)
    rm(futureHeat)
    
    # every some number of runs, take an output snapshot 
    # (sometimes returns rouding errors)
    #print(checkP)
    checkP<-(p/timeStep)%%(timeSnapShot/timeStep)
    #print(checkP)
    if (checkP<=0.5){
      print(p)
      tempFrame<-cbind(tempFrame,temp)
    }
    gc()
    
    #print(p)
    p<-p+timeStep
    
  }
  
} else {
  cat("The square of the spatial step should be no more than twice the time step (currently ",(timeStep/stepSize^2),")\n",sep="")
  cat("Recommended maximum spatial step for this time step: ",((timeStep/0.5)^0.5),"\n", sep="")
  cat("Recommended minimum time step for this spatial step: ",(0.5*stepSize^2),"\n", sep="")
}
temp<-thisHeat/(cp*rho*stepSize*pi*(cylWidth/2)^2)+20

names(tempFrame)<-c("cylMesh",(seq(1,ncol(tempFrame)-1))*timeSnapShot)

tempMelt<-melt(tempFrame,id="cylMesh")


# define a colour palette to be interpolated - use yellow/orange/red and reverse later
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))
pall <- cols(length(unique(tempMelt$variable)))

dev.new()

# plot the data frame up, using the previously mentioned colour palette
timePlot<-ggplot(tempMelt,aes(cylMesh,y=value,color=variable))+geom_line()

timePlot<-timePlot+scale_colour_manual(values=rev(pall))+labs(x="x (mm)", y=expression(paste("Temperature (",degree,"C) ")))+theme(legend.position="none")

print(timePlot)

print(proc.time()-ptm)