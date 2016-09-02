#################################
#################################
##
## Example code for evaluating site fidelity with rASF (rAnimalSiteFidelity)
## using cougar and coyote data.
##
## Author: Peter J. Mahoney, 2016
##
#################################
#################################

################################
## Install package dependencies:
##
## install.packages('bcpa', "rgdal", "rgeos", "maptools", "foreach", 
##                 "doParallel", "lubridate", "ggplot2", "dplyr", "raster")
##
################################

################################
## Key Functions: see examples below!
## 
## Instantiation:
## ClusterXY(xy, dt, proj4string, id, PointID, Data) ~ XY data for clustering function
## ClusterACT((act, dt, id) ~ Raw activity (accelerometer) data
##
## Main Functions:
## visualize_clusters(ClusterXY, nfixes, sbuffer, tbuffer, intime, inPar, nCores) ~ finds clusters
## addDiurnal(ClusterACT, coord, phase, solarDep, localTZ) ~ adds diurnal dummy codes to activity data
## mergeActivity(ClusterACT, F64_clusters, byCluster, aVar, t_int, inPar, nCores) ~ merge activity data with cluster output (case-specific)
## getBCPAa(ClusterACT, variable, range, tau, windowsize, clusterwidth) ~ estimates BCPA for activity only
## getBCPAa(ClusterUN, variable, range, tau, windowsize, clusterwidth) ~ estimates BCPA for activity by cluster
## getBCPAt(ClusterXY, variable, range, tau, windowsize, clusterwidth) ~ estimates BCPA for trajectory only
## getBCPAt(ClusterUN, variable, range, tau, windowsize, clusterwidth) ~ estimates BCPA for trajectory by cluster
## collapseBCPA(ClusterUN, phases, sweep) ~ Collapse BCPA statistics (across clusters) into a data.frame
## addState(ClusterUN, cBCPA) ~ see below, adds state data from K-means
## summaryClusterBehavior(ClusterUN, incRaw, byPoint, byPoly) ~ summarize cluster activity
## addSummaryToSHP(ClusterUN, incRaw, byPoint, byPoly) ~ add summarized data to shp files for export
##
## Export Functions:
## exportClust(ClusterUN, dir, fn, shps = 'All') ~ exports cluster output as ESRI shp files
## plot(ClusterUN, clustids, type='smooth', threshold,   # plot cluster objects, ?plot.bcpa
##     suppressBreakLine, legend, addDiurnal, localTZ, 
##     addState, font.lab, ...)
#################################

#############################
## Integrating Activity:
## NOTE: for the getBCPA* functions, see ?bcpa::WindowSweep
##   for help with arguments.  My implementation only modifies the
##   'data' argument to take 'ClusterXY' and 'ClusterACT' class objects
#############################


######################################
## Cougar Kill Clusters             ##
## Integrating activity sensor data ##
######################################

# Import source file for rASF
# {Currently opens the file from your working directory)
source("./rASF.R")

######
## Visualize clusters
######

# Import point data and convert to SpatialPointDataframe
d <- read.csv('./Puma_Data/F64_2014_GPS.csv')
d$dt <- mdy_hms(d$DateTime, tz = 'MST')
coordinates(d) <- d[, c('Longitude', 'Latitude')]
proj4string(d) <- CRS('+proj=longlat +datum=WGS84')

# Project for accurate measurements in meters
utmProj = CRS("+proj=utm +zone=12 +datum=NAD83")
d <- spTransform(d, utmProj)
str(d)

# Instantiate a ClusterXY object
clustxy <- ClusterXY(xy = coordinates(d), dt = d$dt, 
                     proj4string = utmProj, id = 'F64_2014', 
                     PointID = 'FixID', Data = d@data)
str(clustxy)

# Define clustering parameters
nfixes = 3
sbuffer = 100
tbuffer = 72
intime = T

# In parallel to increase speed of cluster formation?
# Set inPar = FALSE if you want to run in serial
inPar = TRUE
nCores = detectCores() - 1

# Identify clusters and creater ClusterUN class object
F64_clusters <- visualize_clusters(clustxy, nfixes, sbuffer, tbuffer, intime,
                          inPar = inPar, nCores = nCores)

# Looking at the structure of cluster output
str(F64_clusters)
F64_clusters@Data        # A list of length(cluster).  Each element contains Point Data
F64_clusters@Centroid    # Pulls the Centroids as a SpatialPointsDataFrame
F64_clusters@proj4string # Checks the Coordinate Reference System (CRS)
F64_clusters@Params      # Input parameters used

# Cluster summary
summary(F64_clusters)

# Export raw cluster information
# shps = c('All', 'pointobj', 'lineobj', 'polyobj', 'Centroid')
# If 'dir' doesn't exist, it will create a new directory with that name
exportClust(F64_clusters, dir = './F64_14', fn = 'F64_2014_cluster', shps = 'All')

# Import and merge activity (Lotek) with cluster data
# Will need to be modified for logger-specific activity data.
F64_act <- read.csv('./Puma_Data/F64_2014_Activity.csv')
F64_act$dt <- mdy_hms(F64_act$dt, tz = 'MST')

  # Create composite activity metrics for BCPA
F64_act$VeDBA <- sqrt(F64_act$X^2 + F64_act$Y^2)
F64_act$ODBA <- abs(F64_act$X) + abs(F64_act$Y)
head(F64_act)

# Convert to a ClusterACT class object
F64_dAct <- ClusterACT(act = F64_act[, c('X', 'Y', 'VeDBA', 'ODBA')], 
                   dt = F64_act$dt, id = 'F64_2014')
str(F64_dAct)

# Jitter 0's in accelerometer data
# Necessary for BCPA to work properly w/accelerometer data
F64_dAct@act$VeDBA[F64_dAct@act$VeDBA==0] <- jitter(F64_dAct@act$VeDBA[F64_dAct@act$VeDBA==0], 
                                                      amount = 0.000001)
F64_dAct@act$ODBA[F64_dAct@act$ODBA==0] <- jitter(F64_dAct@act$ODBA[F64_dAct@act$ODBA==0], 
                                                      amount = 0.000001)

# Add diurnal period
  # Current implementation only allows for a single coordinate
  # If dealing with wide ranging species, consider another way of adding
  # Diurnal data to your activity data.  Otherwise, a single representative 
  # Coordinate should suffice.
coord <- F64_clusters@OrigData[1, c('Longitude', 'Latitude')]
coord <- SpatialPoints(coord, CRS('+proj=longlat +datum=WGS84'))
F64_dAct <- addDiurnal(F64_dAct, coord=coord,
                          phase = 'withCrepuscular', solarDep=18, localTZ = 'MST')

# Grab the raw activity as a single time series (will be resource intensive
# when running a BCPA if time series is long).
F64_acDat <- mergeActivity(F64_dAct, F64_clusters, byCluster=FALSE, aVar='All')

# Alternatively, grab the raw activity by cluster with some time buffer, t_int (minutes).
F64_acDat <- mergeActivity(F64_dAct, F64_clusters, byCluster=TRUE, aVar='All', 
                           t_int=1440, inPar = T, nCores = detectCores())
str(F64_acDat@rawACT)


#############################
## Perform BCPA using R package 'bcpa':
## NOTE: for the getBCPA* function arguments, see ?bcpa::WindowSweep
##   for help.  My implementation only modifies the 'data' argument to 
##   take 'ClusterXY' and 'ClusterACT' class objects
#############################

##########################################
# Run ONLY ONE of the following options: #
##########################################

# 1) Run BCPA on the raw activity data alone (will be very resource intensive
# when running a BCPA if time series is long). No cluster data involved.
F64_acDat <- getBCPAa(F64_dAct, variable = 'VeDBA', range=0.90, tau=T, 
                    windowsize=30, clusterwidth=1)

# 2) OR, if the complete raw activity data is stored in the ClusterUN object
# Assumes mergeActivity(..., byCluster=FALSE, ...)
F64_acDat <- getBCPAa(F64_acDat, variable = 'VeDBA', range=0.90, tau=T, 
                      windowsize=30, clusterwidth=1)

# 3) OR, for much greater speed, run specific to the clusters identified.
# Assumes mergeActivity(..., byCluster=TRUE,...)
F64_acDat <- getBCPAa(F64_acDat, variable = 'VeDBA', range=0.90, tau=T, 
                      windowsize=30, clusterwidth=1)

##################################
## wrapper for plot.bcpa; see ?plot.bcpa for most arguments
##
## New Plot Arguments: 
## clustids = 'All' or integer for specific cluster
## addDiurnal = logical; add diurnal shading, T/F
## localTZ = string; sometimes necessary for diurnal estimation
## addState = logical; add behavioral state if data available, T/F
##################################

# Plot activity by clusters
plot(F64_acDat, clustids='All', type='smooth', threshold=5,
     suppressBreakLine=T, legend=F, addDiurnal=T, localTZ='MST', 
     addState=F, font.lab = 2, las=2, ylab='VeDBA')
  
  # or plot a subset of clusters
plot(F64_acDat, clustids=c(18,19), type='smooth', threshold=5,
     suppressBreakLine=T, legend=F, addDiurnal=T, localTZ='MST',
     addState=F, font.lab = 2, las=2, ylab='VeDBA')


#################################################
# Estimating behavioral state using K-means
# See package NbClust for a way of identifying the number of clusters, K
# In this case, we assume K=4
#################################################

# Collapse BCPA into a single data.frame object called BCPA
BCPA <- collapseBCPA(F64_acDat, phases=T, sweep=T)
str(BCPA)

# Imputing missing rho by using average values
cBCPA <- BCPA$phases  
cBCPA$rho.hat[is.na(cBCPA$rho.hat)] <- mean(cBCPA$rho.hat, na.rm=T)

kmBCPA <- kmeans(scale(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')]), 
                 4, algorithm='Hartigan-Wong', nstart=10)
kmCenters <- kmBCPA$centers  # pulls the state-specific centers for each metric

# Undoing the standardizing that occured during kmeans estimation
cMeans <- colMeans(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')])
cSDs <- apply(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')], MARGIN=2, sd)
for (nc in 1:ncol(kmCenters)) {
  kmCenters[,nc] <- (kmCenters[,nc] * cSDs[nc]) + cMeans[nc]
}

# Save the K-means centers
write.csv(kmCenters, './F64_14/StateSpecificMetrics_F64.csv')

# Adding Behavioral state as identified by K-means
cBCPA <- cbind(cBCPA, State = kmBCPA$cluster)
F64_datState <- addState(F64_acDat, cBCPA)

# Plot clusters with K-Means state
plot(F64_datState, clustids = c(18), type='smooth', threshold=5,
     suppressBreakLine=T, legend=F, addDiurnal=T, addState=T,
     col.state=c('firebrick4','black','firebrick3','red'),
     font.lab = 2, las=2, ylab='VeDBA', pointLabel=F)

# Pull summary Stats and export as CSV, point, or polygon file
# Includes duration, number of points, proportion of points to a given class
aBCPA <- BCPA$smooth
write.csv(aBCPA, './F64_14/F64_ActivityWindowClust_Summ.csv', row.names=F)

sumList <- summaryClusterBehavior(F64_datState, incRaw=T, byPoint=T, byPoly=T)
write.csv(as.data.frame(sumList$Poly), './F64_14/F64_PredClust_Summ.csv', 
          row.names=F)
write.csv(as.data.frame(sumList$Point), './F64_14/F64_PointClust_Summ.csv', 
          row.names=F)

F64_datState <- addSummaryToSHP(F64_datState, incRaw=T, byPoint=T, byPoly=T)
exportClust(F64_datState, dir = './F64_14/', fn = 'F64_2014_clusterSummary', 
            shps = 'All')





###########################################
## Coyote - 5 minute fixes               ##
## Estimating behavior from trajectories ##
###########################################

# Source rASF
source("./rASF.R")

# Import coyote location data and convert to SpatialPointsDataFrame
d <- read.csv('./Coyote_Data/Coyote_5min.csv')
d$dt <- dmy_hm(paste(d$date, d$hour), tz = 'GMT')
coordinates(d) <- d[, c('x', 'y')]
proj4string(d) <- CRS('+proj=utm +zone=12 +datum=NAD83')

# Store as a ClusterXY class object
clustxy <- ClusterXY(xy = coordinates(d), dt = d$dt, 
                     proj4string = d@proj4string, id = 'Coyote_1', 
                     PointID = 'time', Data = d@data) 

# Add diurnal data to ClusterXY object
clustxy <- addDiurnal(clustxy, phase = 'withCrepuscular', 
                      solarDep=18, localTZ = 'MST')

# Specify cluster parameters
nfixes = 3
sbuffer = 25
tbuffer = 120

## Cluster identification
out <- visualize_clusters(clustxy, nfixes, sbuffer, tbuffer, intime = F, 
                          inPar = T, nCores = detectCores() - 1)

## Merge activity with cluster data, 
## Activity derived from trajectory and should be sorted by time.
## Stored in ClusterUN object (out).
Coy_Act <- mergeActivity(clustxy, out, byCluster=F, t_int=1440, inPar = F)

# NOT RECOMMENDED FOR MOST TRAJECTORY DATA, BUT POSSIBLE FOR SPECIFIC APPLICATIONS
Coy_clustAct <- mergeActivity(clustxy, out, byCluster=T, t_int=1440, inPar = F)

# BCPA for trajectory only
traj <- getBCPAt(clustxy, variable = "V*cos(Theta)", range=0.90, tau=T, 
                 windowsize=30, clusterwidth = 1)

# BCPA for trajectory stored in ClusterUN@rawACT as a single data.frame
traj <- getBCPAt(Coy_Act, variable = "V*cos(Theta)", range=0.90, tau=T, 
               windowsize=30, clusterwidth = 1)

# BCPA for trajectory stored in ClusterUN@rawACT by cluster (NOT RECOMMENDED)
Coy_clustAct@rawACT <- Coy_clustAct@rawACT[1:2]
traj <- getBCPAt(Coy_clustAct, variable = "V*cos(Theta)", range=0.90, tau=T, 
                 windowsize=30, clusterwidth = 1)

# Identify behavioral state from k-means clustering (same as above)
# If estimated by cluster, will need to collapseBCPA() as we do other example
cBCPA <- traj@bcpaACT@breaks$phases

# K-means estimates
kmBCPA <- kmeans(scale(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')]), 
                 3, algorithm='Hartigan-Wong', nstart=10)
kmCenters <- kmBCPA$centers

# Backtransform K-means estimates and export as CSV
cMeans <- colMeans(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')])
cSDs <- apply(cBCPA[, c('mu.hat', 's.hat', 'rho.hat')], MARGIN=2, sd)

for (nc in 1:ncol(kmCenters)) {
  kmCenters[,nc] <- (kmCenters[,nc] * cSDs[nc]) + cMeans[nc]
}
write.csv(kmCenters, paste(dir, i, '_kmCenters.csv', sep=''))

# Add state data to compiled BCPA and reincorporate in ClusterUN object
cBCPA <- cbind(cBCPA, State = kmBCPA$cluster)
traj <- addState(traj, cBCPA)

# Plot clusters with trajectory-based BCPA estimates
plot(traj, clustids = 'All', type='smooth', threshold=5,
     suppressBreakLine=T, legend=F, addDiurnal=T, localTZ = 'MST',
     addState=T, col.state=c('black','firebrick3','red'),
     font.lab = 2, las=2, ylab='Persistence Velocity', pointLabel=F, t_int=1440)

# Add behavioral state data to shape files
traj_summaryState <- addSummaryToSHP(traj, incRaw=T, byPoint=T, byPoly=T)

# Summarize cluster behavior and export
sumList <- summaryClusterBehavior(traj, incRaw=T, byPoint=T, byPoly=T)
write.csv(as.data.frame(sumList$Poly), paste(dir, i, '_PolySumm.csv', sep=''))
write.csv(as.data.frame(sumList$Point), paste(dir, i, '_PointSumm.csv', sep=''))