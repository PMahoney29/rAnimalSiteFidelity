#####################################
##
## Functions for use with rASF v0.1
##
## NOTE: Will be converted to an R package soon!!
##
## Peter Mahoney - 2016
##
## Citation: TBD
##
#####################################

# Package Dependencies
require(bcpa)
require(rgdal)
require(rgeos)
require(maptools)
require(foreach)
require(parallel)
require(doParallel)
require(lubridate)
require(ggplot2)
require(dplyr)
require(raster)

# Set new class union class definitions
setClassUnion("charCRS", c('NULL', "character", "CRS"))
setClassUnion("charFact", c("character", "factor"))
setClassUnion('nullDataList', c("NULL", 'list', 'data.frame'))
setClassUnion('vecMat', c("vector", 'matrix'))
setClassUnion('POSIX', c("POSIXct", 'POSIXlt'))
setClassUnion('datPOSIX', c("data.frame", "POSIXct", 'POSIXlt'))
setClassUnion('nullSpatial', c("NULL", 'list', 'SpatialPolygons', 'SpatialLines', 'SpatialPoints', 'SpatialCollections'))
setClassUnion('listNum', c("NULL", 'list', 'numeric'))
setClassUnion('listChar', c("NULL", 'list', 'character'))
setClassUnion('nullCharNum', c("NULL", 'numeric', 'character'))


# Define new S4 classes
ClusterXY <- setClass(
  # Set the name for the class
  "ClusterXY",

  # Define the slots
  slots = c(
    xy = 'matrix',
    dt = 'datPOSIX',
    proj4string = 'charCRS',
    id = 'charFact',
    PointID = 'numeric',
    Data = 'nullDataList'
  ),

  validity = function(object)
  {
    l <- nrow(object$xy)
    if(length(object$dt) != l || length(object$id) != l ||
       length(object$PointID) != l || nrow(object$Data) != l) {
      return("Unequal length data inputs.  Check the length of xy, dt, id, PointID, and Data.")
    }
    if (is.null(proj4string)) {
      return('Please specify a projection for the xy data following the CRS standard')
    }
    return(TRUE)
  }
)

ClusterHL <- setClass(
  # Set the name for the class
  "ClusterHL",

  # Define the slots
  slots = c(
    ClustID = 'listNum',
    Coords = 'list',
    meandist = 'listNum',
    MinT = 'POSIX',
    MaxT = 'POSIX',
    proj4string = 'charCRS',
    Data = 'nullDataList',
    Polygons = 'nullSpatial'
  )

  #validity=function(object){  }
)

BCPAtraj <- setClass(
  # Set the name for the class
  "BCPAtraj",

  # Define the slots
  slots = c(
    traj = 'nullDataList',
    bcpa = 'ANY',
    breaks = 'list',
    ActID = 'character'
  )

  #validity=function(object){  }
  # Update in getBCPAt('ClusterUN')
)

BCPAact <- setClass(
  # Set the name for the class
  "BCPAact",

  # Define the slots
  slots = c(
    Act = 'nullDataList',
    bcpa = 'ANY',
    breaks = 'list',
    ActID = 'character'
  )

  #validity=function(object){  }
)

ClusterACT <- setClass(
  # Set the name for the class
  "ClusterACT",

  # Define the slots
  slots = c(
    act = 'vecMat',
    dt = 'datPOSIX',
    id = 'nullCharNum'
  ),

  validity=function(object)
  {
    if(length(object@dt) != nrow(object@act)) { 
      return("Unequal length data inputs.  Check the length of act, dt, id, ActID, and Data.")
    }
    return(TRUE)
  }
)

setClassUnion('nullDataListAct', c("NULL", 'list', 'data.frame', 'BCPAtraj', 'BCPAact'))

ClusterUN <- setClass(
  # Set the name for the class
  "ClusterUN",
  
  # Define the slots
  slots = c(
    UClustID = 'listNum',
    IClustIDs = 'listNum',
    CummTotFixes = 'listNum',
    CummPropFixes = 'listNum',
    Area = 'listNum',
    Centroid = 'nullSpatial',
    Data = 'nullDataList',
    OrigData = 'nullDataList',
    SpatialData = 'nullSpatial',
    Params = 'nullDataList',
    proj4string = 'charCRS',
    rawACT = 'nullDataListAct',
    bcpaACT = 'nullDataListAct'
  ),
  
  validity = function(object)
  {
    l <- length(object@UClustID)
    if(length(object@IClustIDs) != l || length(object@Data) != l ||
       length(object@CummTotFixes) != l || length(object@CummPropFixes) != l ||
       length(object@Area) != l || length(object@Centroid) != l) {
      return("Unequal length data inputs.  Check the length of xy, dt, id, PointID, and Data.")
    }
    if (is.null(proj4string)) {
      return('Please specify a projection for the xy data following the CRS standard')
    }
    return(TRUE)
  }
)


## Initialize class instance 
setMethod("initialize", "ClusterXY", function(.Object, xy, dt, proj4string, id, PointID, Data){
  ind <- which(names(Data) == PointID)
  names(Data)[ind] <- 'PointID'
  .Object@xy <- xy
  .Object@dt <- dt
  .Object@proj4string <- proj4string
  .Object@id <- id
  .Object@PointID <- Data$PointID
  .Object@Data <- Data
  .Object
})
setMethod("initialize", "ClusterHL", function(.Object, clust_list){
  .Object@ClustID = clust_list$ClustID
  .Object@Coords = clust_list$Coords
  .Object@meandist = clust_list$meandist
  .Object@MinT = clust_list$MinT
  .Object@MaxT = clust_list$MaxT
  .Object@proj4string = clust_list$Coords[[1]]@proj4string
  .Object@Data = clust_list$Data
  .Object@Polygons = clust_list$Polygon
  .Object
})
setMethod("initialize", "ClusterUN", function(.Object, input){
  .Object@UClustID = input$df$Cluster
  .Object@IClustIDs = input$df$IndPolyIDs
  .Object@CummTotFixes = input$df$CummTotFixes
  .Object@CummPropFixes = input$df$CummPropFixes
  .Object@Area = input$df$area
  .Object@Centroid = input$df$centroid
  .Object@Data = input$df$Data
  .Object@OrigData = NULL
  .Object@SpatialData = input$lip
  .Object@Params = NULL
  .Object@proj4string = NULL
  .Object@rawACT = NULL
  .Object@bcpaACT = NULL
  .Object
})
setMethod("initialize", "BCPAtraj", function(.Object, traj, bcpa, breaks) {
  .Object@traj <- traj
  .Object@breaks <- breaks
  .Object@bcpa <- bcpa
  .Object
})
setMethod("initialize", "BCPAact", function(.Object, act, bcpa, breaks) {
  .Object@Act <- act
  .Object@breaks <- breaks
  .Object@bcpa <- bcpa
  .Object
})
setMethod("initialize", "ClusterACT", function(.Object, act, dt, id) { 
  .Object@act <- act
  .Object@dt <- dt
  .Object@id <- id
  .Object
})

## General Functions
vec_of_time_differences <- function (vertex, times) {
  ## Calculates a vector of temporal distances (hours) between two geographic points

  time_list <- c()
  for (t in 1:length(times)) {
    time <- times[t]
    tim <- difftime(time, vertex, units = "hours")
    tim <- abs(tim[[1]])
    time_list <- c(time_list, tim)
  }
  return(time_list)
}

nearestPoints <- function(c, xy, XY, .Object, params) {
  ## Finds clustered points and draws a convex hull
  
  o_list <- list()
  times <- .Object@dt
  nfixes <- params$nfixes
  tbuffer <- params$tbuffer
  sbuffer <- params$sbuffer
  intime <- params$intime

  # Finds indices for all points within sbuffer
  dvec <- abs(xy[c] - xy)
  ind <- which(dvec <= sbuffer)

  # Pulls out times for all points within sbuffer
  tim <- times[ind]

  # Determines difference in time from point c to
  # every point with di <= sbuffer
  difft <- vec_of_time_differences(times[c], tim)

  # Pulls out indices (or xy points) and times with difftime <= or >= tbuffer
  if (intime == TRUE) {
    ind2 <- ind[difft <= tbuffer]
    tim2 <- tim[difft <= tbuffer]
    meandist <- mean(dvec[ind2])}
  else {
    ind2 <- ind[difft >= tbuffer | difft == 0]
    tim2 <- tim[difft >= tbuffer | difft == 0]
    meandist <- mean(dvec[ind2])}

  # Export list of lists for each point in cluster (ind2), coordinates
  # for those points, mean distance from point c to points ind2,
  # and the earliest and latest time for fixes in cluster.
  if (length(ind2) >= nfixes) {
    o_list$Data <- list(.Object@Data[ind2,])
    o_list$Coords <- list(XY[ind2,])
    o_list$meandist <- meandist
    o_list$MinT <- min(tim2)
    o_list$MaxT <- max(tim2)
    o_list$Polygons <- gConvexHull(XY[ind2,], id = .Object@PointID[c])
    o_list$ClustID <- .Object@PointID[c]
  }
  return(o_list)
}

hullOverlap <- function(p, pols, minTimes, maxTimes) {
  ## Identifies hulls that overlap in space and time
  
  Tint_p = interval(minTimes[p], maxTimes[p])
  g_pols = gIntersects(pols[p], pols, returnDense = F, byid = T)[[1]]

  i_pols = c()
  for (i in g_pols) {
    Tint_i = interval(minTimes[i], maxTimes[i])
    if (int_overlaps(Tint_p, Tint_i)) i_pols <- c(i_pols, i)
  }
  return(i_pols)
}

cluster_intersecting_hulls <- function(int_pols) {
  ## Creates sets for all hulls that intersect. 
  ## The sets contain indices for polygon objects

  clust_sets = list()
  interpols = int_pols

  while (length(interpols) > 0) {
    clust_set = interpols[[1]]
    interpols = interpols[-1]
    if (length(interpols)==0) {
      clust_sets <- append(clust_sets, list(clust_set))
    }
    else {
      Is = c()
      for (i in 1:length(interpols)) {
        if (length(intersect(clust_set, interpols[[i]])) != 0) {
          clust_set = union(clust_set, interpols[[i]])
          Is <- c(Is, i)
        }
      }
      clust_sets <- append(clust_sets, list(clust_set))
      if (!is.null(Is)) interpols <- interpols[-Is]
    }
  }
  return(clust_sets)
}

unionHulls <- function(x, cll, clust_sets, clustxy) {
  ## Unionizes overlapping hulls into a single polygon object
  
  ind <- clust_sets[[x]]
  clust_pols <- cll@Polygons[ind]
  clust_ar <- unlist(lapply(clust_pols, gArea))

  fixes <- cll@Data[ind]
  clust_nums <- unlist(lapply(fixes, nrow))


  # Retrieve indices of sorting first by nums (in reverse) and second by area
  new_ind <- order(-clust_nums, clust_ar)

  ind <- ind[new_ind]
  clust_pols <- clust_pols[new_ind]
  clust_nums <- clust_nums[new_ind]
  clust_ar <- clust_ar[new_ind]

  # Unionizing polygons in order
  lip <- clust_pols[[1]]
  if (length(clust_pols)>1) {
    for (p in 2:length(clust_pols)) {
      lip <- gUnion(lip, clust_pols[[p]], id=as.character(x))
    }
  }
  else {
    lip@polygons[[1]]@ID <- as.character(x)
  }

  # Pulling fix numbers and identifying cummulative fixes
  cumm_uniq_fixes <- list(fixes[[1]][['PointID']])
  if (length(fixes)>1) {
    for (c in 2:length(fixes)) {
      ur <- union(cumm_uniq_fixes[[c-1]], fixes[[c]][['PointID']])
      cumm_uniq_fixes <- append(cumm_uniq_fixes, list(ur))
    }
  }

  #Pulling coordinates for each fix
  uniq_fixes <- unlist(cumm_uniq_fixes[length(cumm_uniq_fixes)])
  uniq_fixes <- uniq_fixes[order(uniq_fixes)]
  act_fixes <- clustxy@Data[which(clustxy@Data[, 'PointID'] %in% uniq_fixes),]
  act_Coords <- clustxy@xy[which(clustxy@Data[, 'PointID'] %in% uniq_fixes),]
  act_fixes$UnionClusterID <- x

  cumm_tot_fixes <- unlist(lapply(cumm_uniq_fixes, length))
  total_fixes <- cumm_tot_fixes[length(cumm_tot_fixes)]
  cumm_prop_fixes <- cumm_tot_fixes / total_fixes

  o <- list()
  o$Cluster <- x
  o$IndPolyIDs <- list(ind[order(ind)])
  o$Data <- list(act_fixes)
  o$CummTotFixes <- list(cumm_tot_fixes)
  o$CummPropFixes <- list(cumm_prop_fixes)

  o$points <- SpatialPointsDataFrame(SpatialPoints(act_Coords, 
                                                   proj4string = clustxy@proj4string), 
                                     data = act_fixes, match.ID=F)

  if (class(lip)[1]=="SpatialLines") {
    lines <- lip
    o$area <- NA
    o$lines <- SpatialLinesDataFrame(lines,
                                   data = data.frame(UClustID = x, Area = o$area),
                                                     match.ID=T)
    o$centroid <- gCentroid(lines)
    row.names(o$centroid) <- x
  }
  if (class(lip)[1]=="SpatialPolygons") {
    polys <- lip
    o$area <- gArea(polys)
    o$polys <- SpatialPolygonsDataFrame(polys,
                                      data = data.frame(UClustID = x, Area = o$area,
                                                        row.names = as.character(x)), match.ID=T)
    o$centroid <- SpatialPointsDataFrame(gCentroid(polys), data = data.frame(UClustID = x), match.ID=F)
    row.names(o$centroid) <- x
  }
  return(o)
}

getDiurnal <- function(coord, startDate, endDate, tzone) {
  nt <- with_tz(seq(startDate, endDate, by = 'mins'), tzone=tzone)
  
  sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
  ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
  dawn <- crepuscule(coord, nt, POSIXct.out=T, direction='dawn', solarDep=18)[,2]
  dusk <- crepuscule(coord, nt, POSIXct.out=T, direction='dusk', solarDep=18)[,2]
  
  Night <- ifelse(nt < dawn | nt > dusk, 1, 0)
  Day <- ifelse(nt >= sr & nt < ss, 1, 0)
  Dawn <- ifelse(nt >= dawn & nt < sr, 1, 0)
  Dusk <- ifelse(nt >= ss & nt < dusk, 1, 0)
  
  nt <- with_tz(nt, tzone=tz(startDate))
  
  return(data.frame(dt = nt, Dawn, Day, Dusk, Night))
}


###################
## Class methods ##
###################

  # Clustering algorithm methods
setGeneric(name='find_clusters',
           def = function(.Object, params)
           {
             standardGeneric('find_clusters')
           })
setMethod('find_clusters', 'ClusterXY', function(.Object, params) 
{
  ## Loops through each fix and lists all other points 
  ## within distance sbuffer and time tbuffer

  XY <- SpatialPoints(.Object@xy, .Object@proj4string)
  
  # Converts xy coordinates into complex numbers for fast distance calculation
  xy <- .Object@xy[,1] + (0 + (0+1i)) * .Object@xy[,2]

  # Function for combining foreach output
  cBind <- function(l1, l2) {
    l1$Data <- append(l1$Data, l2$Data)
    l1$Coords <- append(l1$Coords, l2$Coords)
    l1$meandist <- c(l1$meandist, l2$meandist)
    l1$MinT <- append(l1$MinT, l2$MinT)
    l1$MaxT <- append(l1$MaxT, l2$MaxT)
    l1$Polygon <- c(l1$Polygon, l2$Polygon)
    l1$ClustID <- c(l1$ClustID, l2$ClustID)
    return(l1)
  }

  if (params$inPar==T) {
    cl <- makeCluster(params$nCores)
    registerDoParallel(cl)
    clusterExport(cl, c('vec_of_time_differences', 
                        'nearestPoints'))

    n_list <- foreach(c=1:length(xy),
                      .combine = cBind,
                      .packages = c('rgdal', 'rgeos', 'maptools')) %dopar% nearestPoints(c=c, xy, XY, .Object, params)

    stopCluster(cl)
  }
  else {
    n_list <- foreach(c=1:length(xy),
                         .combine = cBind,
                         .packages = c('rgdal', 'rgeos', 'maptools')) %do% nearestPoints(c=c, xy, XY, .Object, params)
  }

  if (is.list(n_list) & length(n_list) == 0) 
    stop('No clusters found, reset parameters and try again')
  
  tzi <- tz(.Object@dt)
  n_list$MinT <- with_tz(n_list$MinT, tzi)
  n_list$MaxT <- with_tz(n_list$MaxT, tzi)
  return(n_list)
})

setGeneric(name='find_intersect_hulls',
           def = function(.Object, params)
           {
             standardGeneric('find_intersect_hulls')
           })
setMethod('find_intersect_hulls', 'ClusterHL',function(.Object, params) 
{
  ## Create list of lists with indices (for pol objects in clust_list)
  ##  of intersecting hulls
  pols <- SpatialPolygons(unlist(lapply(.Object@Polygons, 
                                        function (x) x@polygons)), 
                          proj4string = .Object@proj4string)

  if (params$inPar==T) {
    cl <- makeCluster(params$nCores)
    clusterExport(cl, c('.Object', 'hullOverlap'), envir = environment())
    clusterEvalQ(cl, library(rgeos))
    clusterEvalQ(cl, library(lubridate))
    int_pols <- parLapply(cl, 1:length(.Object@Polygons),
                          function(x) hullOverlap(p = x, pols = pols,
                                                  minTimes = .Object@MinT, 
                                                  maxTimes = .Object@MaxT))
    stopCluster(cl)
  }
  else {
    int_pols <- lapply(1:length(pols),
                          function(x) hullOverlap(p = x, pols = pols,
                                                  minTimes = .Object@MinT, 
                                                  maxTimes = .Object@MaxT))
  }

  return(int_pols)
})

setGeneric(name='union_hulls_by_cluster',
           def = function(.Object, clust_sets, clustxy, params)
           {
             standardGeneric('union_hulls_by_cluster')
           })
setMethod('union_hulls_by_cluster', 'ClusterHL', function(.Object, clust_sets, 
                                                          clustxy, params) 
{
  ## Sorts and unionizes hulls (most to least fixes, ties settled by smallest to 
  ## largest areas), determines the number and proportion of fixes with every 
  ## union, and determines the cummulative area and proportion with every union

  inPar <- params$inPar
  nCores <- params$nCores

  # Function for combining foreach output
  cBind <- function(l1, l2) {
    l1$Cluster <- append(l1$Cluster, l2$Cluster)
    l1$IndPolyIDs <- append(l1$IndPolyIDs, l2$IndPolyIDs)
    l1$Data <- append(l1$Data, l2$Data)
    l1$CummTotFixes <- append(l1$CummTotFixes, l2$CummTotFixes)
    l1$CummPropFixes <- append(l1$CummPropFixes, l2$CummPropFixes)

    l1$lines <- rbind(l1$lines, l2$lines)
    l1$polys <- rbind(l1$polys, l2$polys)
    l1$points <- rbind(l1$points, l2$points)

    l1$area <- append(l1$area, l2$area)
    l1$centroid <- rbind(l1$centroid, l2$centroid)

    return(l1)
  }

  if (inPar==T) {
    cl <- makeCluster(nCores)
    registerDoParallel(cl)

    output <- foreach(x=1:length(clust_sets),
                      .combine = cBind,
                      .export = c('unionHulls'),
                      .packages = c('rgdal', 'rgeos', 'maptools')) %dopar% unionHulls(x, .Object, clust_sets, clustxy)
    stopCluster(cl)

    output2 <- SpatialCollections(points=output$points, lines=output$lines, 
                                  polygons=output$polys,
                                  proj4string = clustxy@proj4string)
    return(list(df = output[c('Cluster','IndPolyIDs','Data','CummTotFixes',
                              'CummPropFixes','area','centroid')], 
                lip = output2))
  }

  else {
    for (x in 1:length(clust_sets)) {
     o <- unionHulls(x, .Object, clust_sets, clustxy, ID)
     if (x==1) output <- o
     else output <- cBind(output, o)
    }
    output2 <- SpatialCollections(points=output$points, lines=output$lines, 
                                  polygons=output$polys,
                                  proj4string = clustxy@proj4string)
    return(list(df = output[c('Cluster','IndPolyIDs','Data','CummTotFixes',
                              'CummPropFixes','area','centroid')], 
                lip = output2))
  }
})

setGeneric(name='visualize_clusters',
           def = function(.Object, nfixes, sbuffer, tbuffer, intime, inPar, nCores)
           {
             standardGeneric('visualize_clusters')
           })
setMethod('visualize_clusters', 'ClusterXY', function(.Object, nfixes, sbuffer, 
                                                      tbuffer, intime = T, 
                                                      inPar = T, nCores) 
{
  # Full wrapper function for producing clusters
  start <- Sys.time()
  cXY <- .Object
  params <- data.frame(ID = cXY@id[1], nfixes = nfixes, sbuffer = sbuffer, 
                       tbuffer = tbuffer, intime = intime,
                       inPar = inPar, nCores = nCores)

  clust_list <- find_clusters(cXY, params=params)

  if (length(clust_list$ClustID) == 0) {
    print('No clusters found.  Revise input parameters for defining clusters')
    stop
  }
  
  else {
    clust_list <- ClusterHL(clust_list)

    int_pols <- find_intersect_hulls(clust_list, params)

    clust_sets <- cluster_intersecting_hulls(int_pols)

    out <- union_hulls_by_cluster(clust_list, clust_sets, cXY, params)
    out <- ClusterUN(out)

    Date <- Sys.time()
    SimTime <- Date - start
    out@Params <- cbind(Date = Date, SimTime = SimTime, params)
    out@OrigData <- cXY@Data
    out@proj4string <- cXY@proj4string

    return(out)
  }
})

setGeneric(name='exportClust',
           def = function(.Object, dir, fn, shps = 'All')
           {
             standardGeneric('exportClust')
           })
setMethod('exportClust', 'ClusterUN', function(.Object, dir, fn, shps = 'All')
{
  if (shps == 'All') shpList <- c('pointobj', 'lineobj', 'polyobj', 
                                  'Centroid')
  else {
    if (any(shps == 'centroid')) {
      shps <- shps[-which(shps == 'centroid')]
      if (length(shps) > 0) shps <- paste(shps, 'obj', sep='')
      shpList <- c(shps, 'Centroid')
    }
    else shpList <- paste(shps, 'obj', sep='')
  }
  
  dir.create(dir, showWarnings=F)
  
  for (s in shpList) {
    if (.hasSlot(.Object@SpatialData, s)) 
      shp <- slot(.Object@SpatialData, s)
    else shp <- slot(.Object, s)
    if (!is.null(shp)) writeOGR(shp, dir, paste(fn, '_', s, sep=''), 
                                overwrite_layer = T, driver='ESRI Shapefile')
    else print(paste(s, 'is empty.  Confirm slot contains spatial information.'))
  }
  
})

setMethod('summary', 'ClusterUN', function(object)
{
  # Summarize cluster output
  t <- object
  tab <- c()
  tab$DateTime <- t@Params$Date
  tab$RunTime <- t@Params$SimTime
  
  tab$ClusterParameters <- data.frame(t@Params[3:length(t@Params)])
  
  tab$SpatialComposition <- data.frame(TotalClusters = length(t@UClustID), 
                                       TotalPolygons = length(t@SpatialData@polyobj),
                                       TotalLines = length(t@SpatialData@lineobj),
                                       TotalPntsInClusters = nrow(t@SpatialData@pointobj), 
                                       TotalOrigPoints = nrow(t@OrigData))
  return(tab)
})

setGeneric(name='subsetCluster',
           def = function(.Object, startDate, endDate, clusterIDs = NULL)
           {
             standardGeneric('subsetCluster')
           })
setMethod('subsetCluster', 'ClusterUN', function(.Object, startDate = NULL, 
                                                 endDate = NULL, 
                                                 clusterIDs = NULL) 
{
  if(!is.null(clusterIDs) | (!is.null(startDate) & !is.null(endDate))) {
    t <- .Object

    if (!is.null(clusterIDs)) clustIDs <- clusterIDs
    if (!is.null(startDate) & !is.null(endDate)) {
      if (!is.POSIXt(startDate) & !is.POSIXt(endDate)) {
        stop('startDate and/or endDate are not in POSIXt format.  Please convert to POSIXt and re-run.')
      }

      tint <- interval(startDate, endDate)
      td <- do.call('rbind', lapply(t@Data, data.frame, stringsAsFactors=F))
      clustIDs <- unique(td[which(td$datetime %within% tint), 'UnionClusterID'])
    }

      sd <- t@SpatialData
      sd@pointobj <- sd@pointobj[which(sd@pointobj@data$UnionClusterID %in% clustIDs),]
      sd@polyobj <- sd@polyobj[which(sd@polyobj@data$UClustID %in% clustIDs),]
      if (!is.null(sd@lineobj)) 
        sd@lineobj <- sd@lineobj[which(sd@lineobj@data$UClustID %in% clustIDs),]

      out <- ClusterUN(NULL)
      out@UClustID = t@UClustID[clustIDs]
      out@IClustIDs = t@IClustIDs[clustIDs]
      out@CummTotFixes = t@CummTotFixes[clustIDs]
      out@CummPropFixes = t@CummPropFixes[clustIDs]
      out@Area = t@Area[clustIDs]
      out@Centroid = t@Centroid[clustIDs,]
      out@Data = t@Data[clustIDs]
      out@OrigData = t@OrigData
      out@SpatialData <- sd
      out@Params <- c(t@Params, ShowingSubset = paste(clustIDs, collapse = ','))
      out@proj4string <- t@proj4string
  }
  return(out)
})

setGeneric(name='addDiurnal',
           def = function(.Object, phase = c('diurnal', 'withCrepuscular'), 
                          solarDep, localTZ = 'MST',...) {
             standardGeneric('addDiurnal')
           })
setMethod('addDiurnal', 'ClusterACT', function(.Object, coord, 
                                               phase = c('diurnal', 'withCrepuscular'),
                                               solarDep=18, localTZ = 'MST') 
{
  if (class(coord) != 'SpatialPoints' & class(coord) != 'SpatialPointsDataFrame')
    stop('coord must be of class SpatialPoints or SpatialPointsDataFrame, with an assigned CRS')
  if(is.null(proj4string(coord)))
    stop('Assign a non-projected coordinate reference system to coord')
  
  nt <- with_tz(.Object@dt, tzone=localTZ)
  
  if (phase == 'diurnal') {
    sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
    ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
    
    Night <- ifelse(nt < sr | nt >= ss, 1, 0)
    Day <- ifelse(nt >= sr & nt < ss, 1, 0)
    
    .Object@dt <- data.frame(dt = .Object@dt, Night, Day)
  }
  
  if (phase == 'withCrepuscular') {
    sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
    ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
    dawn <- crepuscule(coord, nt, POSIXct.out=T, direction='dawn', solarDep=solarDep)[,2]
    dusk <- crepuscule(coord, nt, POSIXct.out=T, direction='dusk', solarDep=solarDep)[,2]
    
    Night <- ifelse(nt < dawn | nt > dusk, 1, 0)
    Day <- ifelse(nt >= sr & nt < ss, 1, 0)
    Dawn <- ifelse(nt >= dawn & nt < sr, 1, 0)
    Dusk <- ifelse(nt >= ss & nt < dusk, 1, 0)
    
    .Object@dt <- data.frame(dt = .Object@dt, Dawn, Day, Dusk, Night)
  }

  return(.Object)
})
setMethod('addDiurnal', 'ClusterXY', function(.Object,  
                                               phase = c('diurnal', 'withCrepuscular'),
                                               solarDep=18, localTZ = 'MST') 
{
  nt <- with_tz(.Object@dt, tzone=localTZ)
  coord <- as.data.frame(.Object@xy)
  coordinates(coord) <- coord
  proj4string(coord) <- .Object@proj4string
  if(!compareCRS(coord, CRS('+proj=longlat +datum=WGS84')))
    coord <- spTransform(coord, CRS('+proj=longlat +datum=WGS84'))
  
  if (phase == 'diurnal') {
    sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
    ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
    
    Night <- ifelse(nt < sr | nt >= ss, 1, 0)
    Day <- ifelse(nt >= sr & nt < ss, 1, 0)
    
    #.Object@dt <- data.frame(dt = .Object@dt, Night, Day)
    .Object@Data <- cbind(.Object@Data, Night, Day)
  }
  
  if (phase == 'withCrepuscular') {
    sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
    ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
    dawn <- crepuscule(coord, nt, POSIXct.out=T, direction='dawn', solarDep=solarDep)[,2]
    dusk <- crepuscule(coord, nt, POSIXct.out=T, direction='dusk', solarDep=solarDep)[,2]
    
    Night <- ifelse(nt < dawn | nt > dusk, 1, 0)
    Day <- ifelse(nt >= sr & nt < ss, 1, 0)
    Dawn <- ifelse(nt >= dawn & nt < sr, 1, 0)
    Dusk <- ifelse(nt >= ss & nt < dusk, 1, 0)
    
    .Object@Data <- cbind(.Object@Data, Dawn, Day, Dusk, Night)
  }
  
  return(.Object)
})
setMethod('addDiurnal', 'ClusterUN', function(.Object, xyvar=c('x','y'),
                                              phase = c('diurnal', 'withCrepuscular'), 
                                              solarDep=18, localTZ = 'MST') 
{
  
  clust <- .Object@Data
  olist <- list()
  
  for (cl in 1:length(clust)) {
    nt <- with_tz(clust[[cl]]$dt, tzone=localTZ)
    coord <- clust[[cl]][, xyvar]
    coordinates(coord) <- coord[, xyvar]
    proj4string(coord) <- .Object@proj4string
    if(!compareCRS(coord, CRS('+proj=longlat +datum=WGS84')))
      coord <- spTransform(coord, CRS('+proj=longlat +datum=WGS84'))
    
    if (phase == 'diurnal') {
      sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
      ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
      
      Night <- ifelse(nt < sr | nt >= ss, 1, 0)
      Day <- ifelse(nt >= sr & nt < ss, 1, 0)
      
      out <- data.frame(dt = clust[[cl]]@dt, Night, Day)
    }
    
    if (phase == 'withCrepuscular') {
      sr <- sunriset(coord, nt, POSIXct.out=T, direction='sunrise')[,2]
      ss <- sunriset(coord, nt, POSIXct.out=T, direction='sunset')[,2]
      dawn <- crepuscule(coord, nt, POSIXct.out=T, direction='dawn', solarDep=solarDep)[,2]
      dusk <- crepuscule(coord, nt, POSIXct.out=T, direction='dusk', solarDep=solarDep)[,2]
      
      Night <- ifelse(nt < dawn | nt > dusk, 1, 0)
      Day <- ifelse(nt >= sr & nt < ss, 1, 0)
      Dawn <- ifelse(nt >= dawn & nt < sr, 1, 0)
      Dusk <- ifelse(nt >= ss & nt < dusk, 1, 0)
      
      out <- data.frame(dt = clust[[cl]]$dt, Dawn, Day, Dusk, Night)
    }
    
    .Object@Data[[cl]] <- cbind(.Object@Data[[cl]], out)
  }
  
  return(.Object)
})



  # Methods for integrating activity into cluster output
setGeneric(name='mergeActivity',
           def = function(.Object, cDat, byCluster=TRUE, t_int=1440, 
                          inPar = T, nCores = detectCores(),...)
{
  standardGeneric('mergeActivity')
})
setMethod('mergeActivity', 'ClusterACT',function(.Object, cDat, byCluster=TRUE, 
                                                aVar='All', t_int=1440, 
                                                inPar = T, 
                                                nCores = detectCores()) 
{
  ## Test class for cDat
  if (class(cDat) != 'ClusterUN') 
    stop('Cluster object must be of class "ClusterUN"')
  if (cDat@Params$ID != .Object@id)
    stop('Cluster and activity data must come from identical individuals (IDs).')
  
  if (byCluster) {
    if (class(.Object@dt)[1] == 'data.frame') {
      pullACT <- function(cluster, .Object, aVar, t_int) {
        inter <- interval(min(cluster$dt) - (t_int*60), max(cluster$dt) + (t_int*60))
        indx <- which(.Object@dt[,'dt'] %within% inter)
        
        if (aVar=='All') act <- cbind(.Object@dt[indx, ], .Object@act[indx, ])
        else act <- cbind(.Object@dt[indx, ], .Object@act[indx, aVar])
        
        return(act)
      }}
    
    else {
      pullACT <- function(cluster, .Object, aVar, t_int) {
        inter <- interval(min(cluster$dt) - (t_int*60), max(cluster$dt) + (t_int*60))
        indx <- which(.Object@dt %within% inter)
        
        if (aVar=='All') act <- cbind(dt = .Object@dt[indx], .Object@act[indx, ])
        else act <- cbind(dt = .Object@dt[indx], .Object@act[indx, aVar])
        
        return(act)
      }}
    
    if (inPar==T) {
      cl <- makeCluster(nCores)
      clusterExport(cl, ls(envir=globalenv()))
      clusterEvalQ(cl, library(lubridate))
      out <- parLapply(cl, cDat@Data, function(x) pullACT(x, .Object, aVar, t_int))
      stopCluster(cl)
    }
    
    else {
      out <- lapply(cDat@UClustID, 
                    function(x) pullACT(cDat@Data[[x]], .Object, aVar, t_int))
    }
  }

  else {
    if (aVar=='All') out <- cbind(.Object@dt, .Object@act)
    else out <- cbind(.Object@dt, .Object@act[, aVar])
  }  
  
  cDat@rawACT <- out
  return(cDat)
})

setMethod('mergeActivity', 'ClusterXY',function(.Object, cDat, byCluster=TRUE, 
                                                 t_int=1440, 
                                                 inPar = T, 
                                                 nCores = detectCores()) 
{
  ## Test class for cDat
  if (class(cDat) != 'ClusterUN') 
    stop('Cluster object must be of class "ClusterUN"')
  if (cDat@Params$ID != .Object@id)
    stop('Cluster and activity data must come from identical individuals (IDs).')
  
  if (byCluster) {
    pullACT <- function(cluster, vt, t_int) {
      inter <- interval(min(cluster$dt) - (t_int*60), max(cluster$dt) + (t_int*60))
      indx <- which(vt$T.POSIX %within% inter)
      return(vt[indx,])
    }
    
    vt <- GetVT(.Object)
    
    if (inPar==T) {
      cl <- makeCluster(nCores)
      clusterExport(cl, ls(envir=globalenv()))
      clusterEvalQ(cl, library(lubridate))
      out <- parLapply(cl, cDat@Data, function(x) pullACT(x, vt, t_int))
      stopCluster(cl)
    }
    
    else {
      out <- lapply(cDat@UClustID, 
                    function(x) pullACT(cDat@Data[[x]], vt, t_int))
    }
    
    cDat@rawACT <- out
  }

  else cDat@rawACT <- GetVT(.Object)
  
  return(cDat)
})

setGeneric(name='getBCPAa',
           def = function(.Object, variable = "VeDBA", range=0.90, tau=T, 
                          windowsize=50, clusterwidth=1,...)
{
  standardGeneric('getBCPAa')
})
setMethod('getBCPAa', 'ClusterACT', function(.Object, variable = "VeDBA", 
                                             range=0.90, tau=T, windowsize=50, 
                                             clusterwidth=1)
{
  ## Wrapper function for some common functions in Package 'bcpa'
  d.out <- WindowSweep(.Object, variable = variable, range=range, tau=tau, 
                       windowsize=windowsize)
  brs <- ChangePointSummary(d.out, clusterwidth = clusterwidth)
  
  if (class(.Object@dt)[1] == 'data.frame') {
    brs$phases$tp0 <- with_tz(c(min(.Object@dt[,'dt']) - 1, 
                                brs$breaks$middle.POSIX) + 1, 
                              tzone=tz(.Object@dt[,'dt']))
    brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, 
                                max(.Object@dt[,'dt'])), 
                              tzone=tz(.Object@dt[,'dt']))
  }
  else {
    brs$phases$tp0 <- with_tz(c(min(.Object@dt) - 1, 
                                brs$breaks$middle.POSIX) + 1, 
                              tzone=tz(.Object@dt))
    brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, max(.Object@dt)), 
                              tzone=tz(.Object@dt))
  }
  
  out <- BCPAact(act = .Object@act, bcpa = d.out, breaks = brs)#, 
  #clustact = .Object)
  out
})
setMethod('getBCPAa', 'ClusterUN', function(.Object, variable = "VeDBA", 
                                            range=0.90, tau=T, windowsize=50, 
                                            clusterwidth=1)
{
  ## Wrapper function for some common functions in Package 'bcpa'
  d.out <- WindowSweep(.Object, variable=variable, range=range, tau=tau, 
                       windowsize=windowsize)
  
  out <- lapply(1:length(d.out), function (do)
  {
    brs <- ChangePointSummary(d.out[[do]], clusterwidth = clusterwidth)
    
    brs$phases$tp0 <- with_tz(c(min(.Object@rawACT[[do]][,'dt']) - 1, 
                                brs$breaks$middle.POSIX) + 1, 
                              tzone=tz(.Object@rawACT[[do]][,'dt']))
    brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, 
                                max(.Object@rawACT[[do]][,'dt'])), 
                              tzone=tz(.Object@rawACT[[do]][,'dt']))
    
    o2 <- BCPAact(act = NULL, bcpa = d.out[[do]], breaks = brs)
    
    return(o2)
  })
  
  .Object@bcpaACT <- out
  return(.Object)
})

setGeneric(name='getBCPAt',
           def = function(.Object, variable = "V*cos(Theta)", range=0.90, tau=T, 
                          windowsize=50, clusterwidth=1,...)
{
    standardGeneric('getBCPAt')
})
setMethod('getBCPAt', 'ClusterXY', function(.Object, variable = "V*cos(Theta)", 
                                            range=0.90, tau=T, windowsize=50, 
                                            clusterwidth=1)
{
      d.vt <- GetVT(.Object)
  
      ## Wrapper function for some common functions in Package 'bcpa'
      d.out <- WindowSweep(d.vt, variable = variable, range=range, tau=tau, 
                           windowsize=windowsize)
      brs <- ChangePointSummary(d.out, clusterwidth = clusterwidth)
      
      if (class(.Object@dt)[1] == 'data.frame') {
        brs$phases$tp0 <- with_tz(c(min(.Object@dt[,'dt']) - 1, brs$breaks$middle.POSIX) + 1, 
                                  tzone=tz(.Object@dt[,'dt']))
        brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, max(.Object@dt[,'dt'])), 
                                  tzone=tz(.Object@dt[,'dt']))
      }
      else {
        brs$phases$tp0 <- with_tz(c(min(.Object@dt) - 1, brs$breaks$middle.POSIX) + 1, 
                                  tzone=tz(.Object@dt))
        brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, max(.Object@dt)), 
                                  tzone=tz(.Object@dt))
      }
      
      out <- BCPAtraj(traj = d.vt, bcpa = d.out, breaks = brs)
      out
    })
setMethod('getBCPAt', 'ClusterUN', function(.Object, 
                                            variable="V*cos(Theta)", 
                                            range=0.90, tau=T, windowsize=50, 
                                            clusterwidth=1, xyvar = c('x', 'y'))
{
  if (is.null(.Object@rawACT)) {
    dtrack <- .Object@OrigData[,c(xyvar, 'dt')]
    names(dtrack) <- c('X', 'Y', 'Time')
    dvt <- GetVT(dtrack)
    
    ## Wrapper function for some common functions in Package 'bcpa'
    d.out <- WindowSweep(dvt, variable=variable, range=range, tau=tau, 
                         windowsize=windowsize)
  
    brs <- ChangePointSummary(d.out, clusterwidth = clusterwidth)
    
    brs$phases$tp0 <- with_tz(c(min(dtrack[,'Time']) - 1, 
                                brs$breaks$middle.POSIX) + 1, 
                              tzone=tz(dtrack[,'Time']))
    brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, 
                                max(dtrack[,'Time'])), 
                              tzone=tz(dtrack[,'Time']))
    
    out <- BCPAtraj(traj = dvt, bcpa = d.out, breaks = brs) 
  
    .Object@bcpaACT <- out
    return(.Object)
  }
  
  if (is.data.frame(.Object@rawACT)) {
    
    nms <- c('Z.start', 'Z.end', 'S', 'Phi', 'Theta')
    
    if (all(nms %in% names(.Object@rawACT))) { 
      dtrack <- .Object@OrigData['dt']
      dvt <- .Object@rawACT
      
      ## Wrapper function for some common functions in Package 'bcpa'
      d.out <- WindowSweep(dvt, variable=variable, range=range, tau=tau, 
                           windowsize=windowsize)
      
      brs <- ChangePointSummary(d.out, clusterwidth = clusterwidth)
      
      brs$phases$tp0 <- with_tz(c(min(dtrack[,'dt']) - 1, 
                                  brs$breaks$middle.POSIX) + 1, 
                                tzone=tz(dtrack[,'dt']))
      brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, 
                                  max(dtrack[,'dt'])), 
                                tzone=tz(dtrack[,'dt']))
      
      out <- BCPAtraj(traj = dvt, bcpa = d.out, breaks = brs) 
      
      .Object@bcpaACT <- out
      return(.Object)
    }
    else stop('Data in @rawACT is not trajectory based. 
              Consider restoring trajectory data with mergeActivity() or 
              run a separate BCPA on trajectory data and store in @bcpaACT slot.')
  }
  
  if (is.list(.Object@rawACT)) {
    
    nms <- c('Z.start', 'Z.end', 'S', 'Phi', 'Theta')
    
    if (all(nms %in% names(.Object@rawACT[[1]]))) { 
      dtrack <- .Object@OrigData['dt']
      dvt <- .Object@rawACT
      
      ## Wrapper function for some common functions in Package 'bcpa'
      cl <- makeCluster(detectCores()-1)
      clusterEvalQ(cl, setClass(
        "BCPAtraj",
        slots = c(
          traj = 'nullDataList',
          bcpa = 'ANY',
          breaks = 'list',
          ActID = 'character'
        )
      ))
      clusterEvalQ(cl, setClassUnion('nullDataList', c("NULL", 'list', 'data.frame')))
      clusterExport(cl, c('data', 'GetBestBreak', 'GetModels', 
                          'PartitionParameters', 'dtrack', 'WindowSweep', 
                          'variable', 'range', 'tau', 'windowsize', 'clusterwidth',
                          'ChangePointSummary', 'with_tz', 'tz', 'BCPAtraj', 'dvt'), 
                    envir = environment())

      out <- parLapply(cl, dvt, function(x) {
        d.out <- WindowSweep(x, variable=variable, range=range,
                             tau=tau, windowsize=windowsize)
        
        brs <- ChangePointSummary(d.out, clusterwidth = clusterwidth)
        
        brs$phases$tp0 <- with_tz(c(min(dtrack[,'dt']) - 1, 
                                    brs$breaks$middle.POSIX) + 1, 
                                  tzone=tz(dtrack[,'dt']))
        brs$phases$tp1 <- with_tz(c(brs$breaks$middle.POSIX, 
                                    max(dtrack[,'dt'])), 
                                  tzone=tz(dtrack[,'dt']))
        
        d.out <- BCPAtraj(traj = dvt, bcpa = d.out, breaks = brs) 
        return(d.out)
      })
      stopCluster(cl)
      
      .Object@bcpaACT <- out
      return(.Object)
    }
    else stop('Data in @rawACT is not trajectory based. 
              Consider restoring trajectory data with mergeActivity() or 
              run a separate BCPA on trajectory data and store in @bcpaACT slot.')
  }
})

setGeneric(name='collapseBCPA',
           def = function(.Object, phases, sweep)
{
  standardGeneric('collapseBCPA')
})
setMethod('collapseBCPA', 'ClusterUN', function(.Object, phases=T, sweep=T)
{
  if (is.null(.Object@bcpaACT)) 
    stop('No BCPA data stored in ClusterUN object, please run getBCPA*')
  if (length(.Object@bcpaACT) == 1 & 
      length(.Object@bcpaACT) != length(.Object@UClustID))
    stop('No need to collapse, bcpaACT contains the entire trajectory/activity.
         Use name@bcpaACT@breaks$phases.')
  
  bcp <- .Object@bcpaACT

  o <- c()
  if (phases) {
    for (cp in 1:length(bcp)) {
      ph <- bcp[[cp]]@breaks$phases
      ClusterID <- .Object@UClustID[cp]
      Bout <- 1:nrow(ph)
      ph <- cbind(ClusterID, Bout, ph)
      o <- rbind(o, ph)
    }
  }

  o2 <- c()  
  if (sweep) {
    for (cp in 1:length(bcp)) {
      xVar <- bcp[[cp]]@bcpa$x
      xPosix <- bcp[[cp]]@bcpa$t.POSIX
      ClusterID <- .Object@UClustID[cp]
      ox <- cbind(ClusterID, xPosix, xVar)
      
      if ("pp.smooth" %in% names(bcp[[cp]]@bcpa)) {
            smoothPP <- bcp[[cp]]@bcpa$pp.smooth
            names(smoothPP) <- paste('smooth', names(smoothPP), sep='_')
            ox <- cbind(ox, smoothPP)     
      }
      o2 <- rbind(o2, ox)
    }
    o2$xPosix <- as.POSIXct(o2$xPosix, origin=lubridate::origin,
                            tz = tz(bcp[[cp]]@bcpa$t.POSIX))
  }
  
  return(list(ID = .Object@Params$ID, phases = o, sweep = o2))
})

setGeneric(name='addState',
           def = function(.Object, bState)
{
  standardGeneric('addState')
})
setMethod('addState', 'ClusterUN', function(.Object, bState)
{
  if (is.null(.Object@bcpaACT)) 
    stop('No BCPA data stored in ClusterUN object, please run getBCPA*')
  
  if (length(.Object@bcpaACT)==1) {
    state <- bState[, 'State']
    .Object@bcpaACT@breaks$phases$State <- state
  }
  
  else {
    for (cp in 1:length(.Object@bcpaACT)) {
      state <- bState[bState$ClusterID==cp, 'State']
      .Object@bcpaACT[[cp]]@breaks$phases$State <- state
    }
  }
  
  return(.Object)
})

setGeneric(name='intPointsRAW',
           def = function(.Object)
{
  standardGeneric('intPointsRAW')
})
setMethod('intPointsRAW', 'ClusterUN', function(.Object)
{
  if (is.null(.Object@rawACT)) 
    stop('No RAW data stored in ClusterUN object')
  dat <- .Object@Data
  raw <- .Object@rawACT
  
  if(!is.data.frame(raw) & is.list(raw) & length(dat) == length(raw)) {
    if ('T.POSIX' %in% names(raw[[1]])) tName <- 'T.POSIX'
    else tName <- 'dt'
    
    cl <- makeCluster(detectCores()-1)
    clusterExport(cl, c('dat', 'raw'), 
                  envir = environment())
    out <- parLapplyLB(cl, 1:length(dat), function (i)
    {
      di <- dat[[i]]
      bci <- raw[[i]]
      
      ou <- c()
      for (r in 1:nrow(di)) {
        rec <- di[r,]
        rdt <- rec[, 'dt']
        indx <- which(abs(bci$dt - rdt) == min(abs(bci$dt - rdt)))
        
        if (tName == 'T.POSIX') rawVariable <- bci[indx, -c(1:2)]
        else rawVariable <- bci[indx, -1]
  
        o <- cbind(rec, rawVariable)
        ou <- rbind(ou, o)
      }
      return(ou)
    })
    stopCluster(cl)
  }
  
  if(is.data.frame(raw)) {
    if ('T.POSIX' %in% names(raw)) tName <- 'T.POSIX'
    else tName <- 'dt'
    
    cl <- makeCluster(detectCores()-1)
    clusterExport(cl, c('dat', 'raw', 'tName'), 
                  envir = environment())
    out <- parLapplyLB(cl, 1:length(dat), function (i)
    {
      di <- dat[[i]]
      bci <- raw
      
      ou <- c()
      for (r in 1:nrow(di)) {
        rec <- di[r,]
        rdt <- rec[, 'dt']
        indx <- which(abs(bci[,tName] - rdt) == min(abs(bci[,tName] - rdt)))[1]
        
        if (tName == 'T.POSIX') rawVariable <- bci[indx, -c(1:2)]
        else rawVariable <- bci[indx, -1]
        
        o <- cbind(rec, rawVariable)
        ou <- rbind(ou, o)
      }
      return(ou)
    })
    stopCluster(cl)
  }
  
  return(out)
})

setGeneric(name='intPointsBCPA',
           def = function(.Object)
{
  standardGeneric('intPointsBCPA')
})
setMethod('intPointsBCPA', 'ClusterUN', function(.Object)
{
  if (is.null(.Object@bcpaACT)) 
    stop('No BCPA data stored in ClusterUN object, please run getBCPA*')
  dat <- .Object@Data
  bcp <- .Object@bcpaACT
  
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c('dat', 'bcp'), 
                envir = environment())
  out <- parLapplyLB(cl, 1:length(dat), function (i)
  {
    di <- dat[[i]]
    
    if (length(bcp)==1) bci <- bcp
    else bci <- bcp[[i]]
    
    ou <- c()
    for (r in 1:nrow(di)) {
      rec <- di[r,]
      rdt <- rec[, 'dt']
      indx <- which(abs(bci@bcpa$t.POSIX - rdt) == 
                      min(abs(bci@bcpa$t.POSIX - rdt)))[1]
      bcpaVariable <- bci@bcpa$x[indx]
      bcpaT <- bci@bcpa$t[indx]
      
      indxPh <- which(bci@breaks$phases$t0 <= bcpaT & bci@breaks$phases$t1 > bcpaT)
      bcpaState <- bci@breaks$phases[indxPh, 'State']
      
      o <- cbind(rec, bcpaVariable, bcpaState)
      if (!is.null(bci@bcpa$pp.smooth)) {
        bcpaSmooth <- bci@bcpa$pp.smooth[indx,]
        o <- cbind(o, bcpaSmooth)
      }
      ou <- rbind(ou, o)
    }
    return(ou)
    })
  stopCluster(cl)
  
  return(out)
})

setGeneric(name='summaryClusterBehavior',
           def = function(.Object, incRaw, byPoint, byPoly)
{
  standardGeneric('summaryClusterBehavior')
})
setMethod('summaryClusterBehavior', 'ClusterUN', 
          function(.Object, incRaw=T, byPoint=T, byPoly=T)
{
  if (is.null(.Object@bcpaACT)) 
    stop('No BCPA data stored in ClusterUN object, please run getBCPA*')
            
  if (is.null(.Object@rawACT) & incRaw==T) 
    stop('No raw data stored in ClusterUN object') 
  
  if (byPoint==F & byPoly==F) 
    stop('Please choose the level of summary output: point and/or polygon') 

  if (!is.null(.Object@rawACT) & incRaw==T) {
    .Object@Data <- intPointsRAW(.Object)
  }
  
  iP_bcpa <- intPointsBCPA(.Object)
            
  if(byPoint) {
    oP <- do.call("rbind", iP_bcpa)
  }
  
  if(byPoly) {
    outPolData <- c()
    
    for (cl in .Object@UClustID) {
      cent <- as.data.frame(.Object@Centroid)[.Object@Centroid$UClustID==cl,]
      
      names(cent)[2:3] <- c('Cent_x', 'Cent_y')
      
      tdur <- as.numeric(difftime(max(.Object@Data[[cl]]$dt), 
                                  min(.Object@Data[[cl]]$dt),
                                  units='hours'))
      nfixes <- nrow(.Object@Data[[cl]])
      area <- .Object@Area[[cl]]
      
      outPolData <- rbind(outPolData, 
                          cbind(cent, T_duration = tdur, 
                                NumFixes = nfixes, Area = area))
    }
    if (is.list(.Object@bcpaACT)) {
      if ('State' %in% names(.Object@bcpaACT[[1]]@breaks$phases)) {
        states <- unique(collapseBCPA(.Object, phases=T, sweep=F)$phases$State)
        states <- states[order(states)]    
        pStates <- c()
        for (ic in 1:length(iP_bcpa)) {
          cis <- iP_bcpa[[ic]]$bcpaState
          oProp <- c()
          for (s in 1:length(states)) {
            oProp <- cbind(oProp, sum(cis==states[s]) / length(cis))
          }
          oProp <- as.data.frame(oProp)
          names(oProp) <- paste('PropSt_', states, sep='')
          
          pStates <- rbind(pStates, oProp)
        }
        outPolData <- cbind(outPolData, pStates, 
                            row.names=1:nrow(outPolData))
      }
    }
    
    if (isS4(.Object@bcpaACT)) {
      if ('State' %in% names(.Object@bcpaACT@breaks$phases)) {
        states <- unique(.Object@bcpaACT@breaks$phases$State)
        states <- states[order(states)]    
        pStates <- c()
        for (ic in 1:length(iP_bcpa)) {
          cis <- iP_bcpa[[ic]]$bcpaState
          oProp <- c()
          for (s in 1:length(states)) {
            oProp <- cbind(oProp, sum(cis==states[s]) / length(cis))
          }
          oProp <- as.data.frame(oProp)
          names(oProp) <- paste('PropSt_', states, sep='')
          
          pStates <- rbind(pStates, oProp)
        }
        outPolData <- cbind(outPolData, pStates, 
                            row.names=1:nrow(outPolData))
      }
    }
    
    if ('Night' %in% names(iP_bcpa[[1]])) {
      PropNight <- c()
      for (ic in 1:length(iP_bcpa)) {
        cin <- iP_bcpa[[ic]]$Night
        PropNight <- rbind(PropNight, sum(cin==1) / length(cin))
      }
      
      Diurnal <- data.frame(PropNight = PropNight)
      
      if ('Dusk' %in% names(iP_bcpa[[1]])) {
        PropCrep <- c()
        for (ic in 1:length(iP_bcpa)) {
          cidu <- iP_bcpa[[ic]]$Dusk
          cida <- iP_bcpa[[ic]]$Dawn
          PropCrep <- rbind(PropCrep, sum(cidu==1 | cida==1) / length(cidu))
        }
        Diurnal <- cbind(Diurnal, PropCrep = PropCrep)
      }
      outPolData <- cbind(outPolData, Diurnal, 
                          row.names=1:nrow(outPolData))
    }
    }
  
  object <- list()
  object$Point <- oP
  object$Poly <- outPolData
  return(object)
})

setGeneric(name='addSummaryToSHP',
           def = function(.Object, incRaw=T, byPoint=T, byPoly=T)
{
  standardGeneric('addSummaryToSHP')
})
setMethod('addSummaryToSHP', 'ClusterUN', 
          function(.Object, incRaw=T, byPoint=T, byPoly=T)
{
  if (is.null(.Object@bcpaACT)) 
    stop('No BCPA data stored in ClusterUN object, please run getBCPA*')
            
  if (is.null(.Object@rawACT) & incRaw==T) 
    stop('No raw data stored in ClusterUN object') 
            
  if (byPoint==F & byPoly==F) 
    stop('Please choose the level of summary output: point and/or polygon') 
            
  if (!is.null(.Object@rawACT) & incRaw==T) {
    .Object@Data <- intPointsRAW(.Object)
  }
            
  iP_bcpa <- intPointsBCPA(.Object)
            
  if(byPoint) {
    .Object@SpatialData@pointobj@data <- do.call("rbind", iP_bcpa)
  }
            
  if(byPoly) {
    outPolData <- c()
              
    for (cl in .Object@UClustID) {
      cent <- as.data.frame(.Object@Centroid)[.Object@Centroid$UClustID==cl,]
                
      names(cent)[2:3] <- c('Cent_x', 'Cent_y')
                
      tdur <- as.numeric(difftime(max(.Object@Data[[cl]]$dt), 
                                  min(.Object@Data[[cl]]$dt),
                                  units='hours'))
      nfixes <- nrow(.Object@Data[[cl]])
      area <- .Object@Area[[cl]]
                
      outPolData <- rbind(outPolData,
                          cbind(cent, T_duration = tdur,
                                NumFixes = nfixes, Area = area))
    }
    if (is.list(.Object@bcpaACT)) {
      if ('State' %in% names(.Object@bcpaACT[[1]]@breaks$phases)) {
        states <- unique(collapseBCPA(.Object, phases=T, sweep=F)$phases$State)
        states <- states[order(states)]    
        pStates <- c()
        for (ic in 1:length(iP_bcpa)) {
          cis <- iP_bcpa[[ic]]$bcpaState
          oProp <- c()
          for (s in 1:length(states)) {
            oProp <- cbind(oProp, sum(cis==states[s]) / length(cis))
          }
          oProp <- as.data.frame(oProp)
          names(oProp) <- paste('PropSt_', states, sep='')
          
          pStates <- rbind(pStates, oProp)
        }
        outPolData <- cbind(outPolData, pStates, 
                            row.names=1:nrow(outPolData))
      }
    }
              
    if (isS4(.Object@bcpaACT)) {
      if ('State' %in% names(.Object@bcpaACT@breaks$phases)) {
        states <- unique(.Object@bcpaACT@breaks$phases$State)
        states <- states[order(states)]    
        pStates <- c()
        for (ic in 1:length(iP_bcpa)) {
          cis <- iP_bcpa[[ic]]$bcpaState
          oProp <- c()
          for (s in 1:length(states)) {
            oProp <- cbind(oProp, sum(cis==states[s]) / length(cis))
          }
          oProp <- as.data.frame(oProp)
          names(oProp) <- paste('PropSt_', states, sep='')
          
          pStates <- rbind(pStates, oProp)
        }
        outPolData <- cbind(outPolData, pStates, 
                            row.names=1:nrow(outPolData))
      }
    }
    
    if ('Night' %in% names(iP_bcpa[[1]])) {
      PropNight <- c()
      for (ic in 1:length(iP_bcpa)) {
        cin <- iP_bcpa[[ic]]$Night
        PropNight <- rbind(PropNight, sum(cin==1) / length(cin))
      }
      
      Diurnal <- data.frame(PropNight = PropNight)
      
      if ('Dusk' %in% names(iP_bcpa[[1]])) {
        PropCrep <- c()
        for (ic in 1:length(iP_bcpa)) {
          cidu <- iP_bcpa[[ic]]$Dusk
          cida <- iP_bcpa[[ic]]$Dawn
          PropCrep <- rbind(PropCrep, sum(cidu==1 | cida==1) / length(cidu))
        }
        Diurnal <- cbind(Diurnal, PropCrep = PropCrep)
      }
      outPolData <- cbind(outPolData, Diurnal, 
                          row.names=1:nrow(outPolData))
    }
    
    if (!is.null(.Object@SpatialData@polyobj))
      .Object@SpatialData@polyobj@data <- merge(.Object@SpatialData@polyobj@data[1], outPolData)
    if (!is.null(.Object@SpatialData@lineobj))
      .Object@SpatialData@lineobj@data <- merge(.Object@SpatialData@lineobj@data[1], outPolData)
    #if (!is.null(.Object@Centroid))
    #  .Object@Centroid@data <- merge(.Object@Centroid@data[1], outPolData)
  }
  return(.Object)
})


#################################################################
####### Wrappers for Package 'bcpa' functions, v1.1 #############
#################################################################
setMethod('GetVT', 'ClusterXY', function(Data = .Object, units = "hour", 
                                         skiplast = TRUE) 
{
  Data <- data.frame(X = Data@xy[,1], Y = Data@xy[,2], Time = Data@dt)
  
  ### Code from Package 'bcpa' v1.1 ##
  if (!"Z" %in% names(Data))
    Data$Z <- Data$X + (0 + (0+1i)) * Data$Y
  Z.start <- Data$Z[-nrow(Data)]
  Z.end <- Data$Z[-1]
  S <- Mod(diff(Data$Z))
  Phi <- Arg(diff(Data$Z))
  Theta <- c(NA, diff(Phi))
  T.POSIX <- Data$Time[-nrow(Data)] + diff(Data$Time)/2
  if (inherits(Data$Time, "POSIXt")) {
    Data$Time <- as.numeric(Data$Time - Data$Time[1])
    Data$Time <- Data$Time/ifelse(units == "sec", 1, 
                                  ifelse(units == "min", 60, 
                                         ifelse(units == "hour", 60 * 60, 
                                                ifelse(units == "day", 
                                                       60 * 60 * 24, 
                                                       stop("Invalid time unit.")))))
  }
  T.start <- Data$Time[-nrow(Data)]
  T.end <- Data$Time[-1]
  dT <- T.end - T.start
  V <- S/as.vector(dT)
  T.mid <- (T.start + T.end)/2
  VT.table <- data.frame(Z.start, Z.end, S, Phi, Theta, T.start,
                         T.end, T.mid, dT, V, T.POSIX)
  if (skiplast)
    VT.table <- VT.table[-1, ]
  return(VT.table)
  #####################################
})

setMethod('WindowSweep', 'ClusterACT', function (data = .Object, variable, 
                                                 windowsize = 50, windowstep = 1, 
                                                 K = 2, tau = TRUE, range = 0.6, 
                                                 progress = TRUE, plotme = FALSE,
                                                 ...) 
{
  ## MODIFIED FROM 'bcpa' FUNCTION WindowSweep
  
  x <- eval(parse(text = variable), data@act)
  
  if(class(data@dt) == 'data.frame')
    t.POSIX <- data@dt[,'dt']
  else 
    t.POSIX <- data@dt
  
  t <- as.numeric(difftime(t.POSIX, min(t.POSIX), units = 'mins'))
  
  low <- seq(1, (length(t) - windowsize), windowstep)
  hi <- low + windowsize
  if (progress)
    pb <- txtProgressBar(min = 0, max = length(low), style = 3)
  for (i in 1:length(low)) {
    myx <- x[low[i]:hi[i]]
    myt <- t[low[i]:hi[i]]
    myestimate <- GetBestBreak(myx, myt, range, tau = tau)
    breakpoint <- myestimate[1]
    tbreak <- myestimate[2]
    allmodels <- GetModels(myx, myt, breakpoint, K, tau)
    mymodel <- allmodels[allmodels[, 3] == min(allmodels[, 3]), ]
    mymodel <- c(mymodel, Break = tbreak)
    if (i == 1)
      estimates <- mymodel
    else estimates <- rbind(estimates, mymodel)
    if (plotme) {
      plot.ts(t, x, type = "l", col = "grey")
      lines(t, x, type = "l")
      lines(myt, myx, col = "green")
      abline(v = tbreak)
    }
    if (progress & i%%10 == 0)
      setTxtProgressBar(pb, i)
  }
  if (progress)
    close(pb)
  windowsweep <- list(ws = data.frame(estimates, row.names = 1:nrow(estimates)),
                      x = x, t = t, t.POSIX = t.POSIX, windowsize = windowsize,
                      windowstep = windowstep)
  windowsweep$pp.smooth <- PartitionParameters(windowsweep,
                                               type = "smooth", ...)
  class(windowsweep) <- "bcpa"
  return(windowsweep)
})

setMethod('WindowSweep', 'ClusterUN', function (data = .Object, variable, 
                                                windowsize = 50, windowstep = 1, 
                                                K = 2, tau = TRUE, range = 0.6, 
                                                progress = TRUE, plotme = FALSE,
                                                ...) 
{

  ## MODIFIED FROM 'bcpa' FUNCTION WindowSweep
  data <- data@rawACT
  
  if(is.null(data))
    stop('No raw data assigned to the rawAct slot.')  
  
  if (is.data.frame(data)) data <- list(data)

  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c('data', 'GetBestBreak', 'GetModels', 'PartitionParameters'), 
                envir = environment())
  out <- parLapplyLB(cl, 1:length(data), function (di=x, variable, windowsize, 
                                                   windowstep, K, tau, range, 
                                                   progress, plotme,...)
  {
    d <- data[[di]]
    x <- eval(parse(text = variable), d)
    
    t.POSIX <- d[,'dt']
    
    t <- as.numeric(difftime(t.POSIX, min(t.POSIX), units = 'mins'))
    
    low <- seq(1, (length(t) - windowsize), windowstep)
    hi <- low + windowsize
    if (progress)
      pb <- txtProgressBar(min = 0, max = length(low), style = 3)
    for (i in 1:length(low)) {
      myx <- x[low[i]:hi[i]]
      myt <- t[low[i]:hi[i]]
      myestimate <- GetBestBreak(myx, myt, range, tau = tau)
      breakpoint <- myestimate[1]
      tbreak <- myestimate[2]
      allmodels <- GetModels(myx, myt, breakpoint, K, tau)
      mymodel <- allmodels[allmodels[, 3] == min(allmodels[, 3]), ]
      mymodel <- c(mymodel, Break = tbreak)
      if (i == 1)
        estimates <- mymodel
      else estimates <- rbind(estimates, mymodel)
      if (plotme) {
        plot.ts(t, x, type = "l", col = "grey")
        lines(t, x, type = "l")
        lines(myt, myx, col = "green")
        abline(v = tbreak)
      }
      if (progress & i%%10 == 0)
        setTxtProgressBar(pb, i)
    }
    if (progress)
      close(pb)
    windowsweep <- list(ws = data.frame(estimates, row.names = 1:nrow(estimates)),
                        x = x, t = t, t.POSIX = t.POSIX, windowsize = windowsize,
                        windowstep = windowstep)
    windowsweep$pp.smooth <- PartitionParameters(windowsweep,
                                                 type = "smooth",...)
    class(windowsweep) <- "bcpa"
    return(windowsweep)
  }, variable, windowsize, windowstep, K, tau, range, progress, plotme,...)
  stopCluster(cl)
  
  return(out)
})

setMethod('plot', 'ClusterUN', function(x, clustids='All',
                                        pointLabel = T, missing = T,
                                        suppressBreakLine = T, 
                                        type = c("smooth", "flat")[1],
                                        threshold = 3, clusterwidth = 1,
                                        addDiurnal = T, localTZ = Sys.timezone(),
                                        addState = T,
                                        col.state = c('blue', 'yellow', 'orange', 'red'),
                                        col.cp = rgb(0.5, 0, 0.5, 0.5), 
                                        pt.cex = 0.5, legend = TRUE,
                                        rho.where = "topleft", mu.where = "nowhere", 
                                        col.sd = "red",
                                        col.mean = "black", t_int=1440, ...)   
{
  ###########################################################
  ################### WRAPPER FUNCTION FOR PLOT.BCPA() ######
  ###########################################################
  cluster <- x  
  
  if (length(cluster@bcpaACT) > 1) {
    if (class(cluster@bcpaACT[[1]]) != 'BCPAtraj' & 
        class(cluster@bcpaACT[[1]]) != 'BCPAact')
      stop("Activity must be of class 'BCPAtraj or BCPAact'")
    
    windowsweep <- cluster@bcpaACT
    
    if (clustids[1] == 'All') clustids <- 1:length(cluster@UClustID)
    
    for (clustid in clustids) {
      clust <- cluster@Data[[clustid]]
      acti <- windowsweep[[clustid]]
      
      x <- acti@bcpa$x
      t.POSIX <- acti@bcpa$t.POSIX
      t <- acti@bcpa$t
      
      startDate <- min(t.POSIX)
      endDate <- max(t.POSIX)
      
      if (!cluster@Params$intime) ty = 'Out Time'
      else ty = 'In Time'
      
      old.par <- par(mgp=c(4,1,0), mar=c(6.25,5,4,2)+0.1, ask=T)
      plot(t.POSIX, x, type = "n", main = paste('Cluster ', clustid, ' (', ty, ')', sep=''),
           xlab=paste('Time', ' (', tz(t.POSIX), ')', sep=''), xaxt='n', ...)
      
      if (addState) {
        for (ph in 1:nrow(acti@breaks$phases)) {
          dtseq <- seq(acti@breaks$phases[ph, 'tp0'], 
                       acti@breaks$phases[ph, 'tp1'],
                       by = 'mins')
          st <- acti@breaks$phases$State[ph]
          rug(dtseq, col=col.state[st], side=3, ticksize=-0.05)
        }
      }
      
      tdiff <- difftime(max(t.POSIX), min(t.POSIX), unit = 'sec')
      
      if (tdiff <= 3600) {
        ati <- seq(round(min(t.POSIX), 'secs'), max(t.POSIX), by = 'mins')
        
        while (length(ati) > 30) {
          ati <- ati[seq(1, length(ati), by=2)]
        }
        
        axis.POSIX(1, at=ati, format=("%H:%M:%S"), las=2)
      }
      else {
        ati <- seq(round(min(t.POSIX), 'hours'), max(t.POSIX), by = 'hours')
        
        while (length(ati) > 30) {
          ati <- ati[seq(1, length(ati), by=2)]
        }
        
        axis.POSIXct(1, at=ati, format=("%H:%M"), las=2)
      }
      
      if (addDiurnal) {
        coord <- cluster@Centroid[clustid,]
        if(!compareCRS(coord, CRS('+proj=longlat +datum=WGS84')))
          coord <- spTransform(coord, CRS('+proj=longlat +datum=WGS84'))
        
        dn <- getDiurnal(coord, startDate, endDate, tzone=localTZ)
        
        nigh <- which(dn$Night == 1)
        nnigh <- which(dn$Night == 0)
        
        while (length(nigh) > 0) {
          nnigh <- nnigh[nnigh > nigh[1]]
          
          if (length(nnigh) > 0) {
            nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nnigh[1]-1])
            nigh <- nigh[nigh > nnigh[1]-1]
          }
          else {
            nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nigh[length(nigh)]])
            nigh = NULL
          }
          rect(nout$ns, min(x), nout$ne, max(x), col = 'dark gray', border=NA)
        }
        
        da <- which(dn$Dawn == 1)
        dda <- which(dn$Dawn == 0)
        du <- which(dn$Dusk == 1)
        ddu <- which(dn$Dusk == 0)
        
        while (length(da) > 0) {
          dda <- dda[dda > da[1]]
          
          if (length(dda) > 0) {
            daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[dda[1]-1])
            da <- da[da > dda[1]-1]
          }
          
          else {
            daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[da[length(da)]])
            da <- NULL
          }
          rect(daout$das, min(x), daout$dae, max(x), col = 'light gray', border=NA)
        }
        
        while (length(du) > 0) {
          ddu <- ddu[ddu > du[1]]
          
          if (length(ddu) > 0) {
            duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[ddu[1]-1])
            du <- du[du > ddu[1]-1]
          }
          
          else {
            duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[du[length(du)]])
            du <- NULL
          }
          rect(duout$dus, min(x), duout$due, max(x), col = 'light gray', border=NA)
        }
      }
      
      ## Add Cluster point location
      abline(v = clust$dt, lwd = 2, col = 'black')
      if (pointLabel) {
        mtext(text = clust$PointID, las = 2, side = 3, outer = FALSE, at=clust$dt)
      }
      
      if (missing) {
        tint2 <- interval(min(clust$dt), max(clust$dt))
        allFixes <- cluster@OrigData$dt[which(cluster@OrigData$dt %within% tint2)]
        fixint <- min(round(diff(allFixes)))
        
        seqt <- seq(min(clust$dt), max(clust$dt), fixint)
        
        missFix <- c()
        for (st in 1:length(seqt)) {
          if (any(abs(difftime(seqt[st], allFixes, units='secs')) <= 60*5)) next
          else missFix <- c(missFix, seqt[st])
        }
        
        if (!is.null(missFix)) {
          missFix <- as.POSIXct(missFix, origin=lubridate::origin, tz = tz(cluster@OrigData$dt))
          abline(v = missFix, lwd = 2, col = 'red')
        }
      }
      
      lines(t.POSIX, x, col = "grey")
      
      ws <- acti@bcpa$ws
      
      if (type == "smooth") {
        if ("pp.smooth" %in% names(acti@bcpa))
          pp <- acti@bcpa$pp.smooth
        else pp <- PartitionParameters(acti@bcpa, type = "smooth")
        GoodBreaks <- ws$Break[ws$Model > 0]
        GoodBreaks <- as.data.frame(table(GoodBreaks))
        GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[, 1], t)])
        GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[, 1]))
        GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold, ]
        if (!suppressBreakLine)
          abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold * 2, col = col.cp)
        rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
        rho.int <- round(rho.scaled * 999 + 1)
        palette(topo.colors(1000))
        points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
               cex = pt.cex, lwd = 0.5)
        lines(t.POSIX, pp$mu.hat, lwd = 1.5, col = col.mean)
        lines(t.POSIX, pp$mu.hat + pp$s.hat, col = col.sd, lwd = 1.5)
        lines(t.POSIX, pp$mu.hat - pp$s.hat, col = col.sd, lwd = 1.5)
        rho.hat <- pp$rho.hat
      }
      
      if (type == 'flat')
        stop('Currently, only the smooth bcpa is permitted with cluster data')
      
      if (legend) {
        legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
        legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE),
                                  length = 5), 2)
        if (rho.where != "nowhere")
          legend(rho.where, bg = "white", legend = c(expression(hat(rho)),
                                                     legend.rhos), pch = 19, ncol = 3, col = c(0,
                                                                                               legend.cols), xjust = 0.5, yjust = 0.5)
        if (mu.where != "nowhere")
          legend(mu.where, bg = "white", legend = c(expression(hat(mu)),
                                                    expression(hat(mu) %+-% hat(sigma))), lty = 1,
                 lwd = 2:1, col = c("black", "red"), xjust = 0.5,
                 yjust = 0.5)
      }
      par(old.par)
    }
    palette("default")
  }
  
  if (length(cluster@bcpaACT) == 1) {
    if (class(cluster@bcpaACT) != 'BCPAtraj' & 
        class(cluster@bcpaACT) != 'BCPAact')
      stop("Activity must be of class 'BCPAtraj or BCPAact'")  
    
    acti <- cluster@bcpaACT
    windowsweep <- acti@bcpa
    t.POSIXt <- windowsweep$t.POSIX
    
    if (clustids[1] == 'All') clustids <- 1:length(cluster@UClustID)
    
    for (clustid in clustids) {
      clust <- cluster@Data[[clustid]]

      startDate <- min(clust$dt) - (t_int * 60)
      endDate <- max(clust$dt) + (t_int * 60)
      tint <- interval(startDate, endDate)
      inds <- which(t.POSIXt %within% tint)
      
      x <- windowsweep$x[inds]
      t.POSIX <- t.POSIXt[inds]
      t <- windowsweep$t[inds]
      
      if (!cluster@Params$intime) ty = 'Out Time'
      else ty = 'In Time'
      
      old.par <- par(mgp=c(4,1,0), mar=c(6.25,5,4,2)+0.1, ask=T)
      plot(t.POSIX, x, type = "n", main = paste('Cluster ', clustid, ' (', ty, ')', sep=''),
           xlab=paste('Time', ' (', tz(t.POSIX), ')', sep=''), xaxt='n', ...)
      
      if (addState) {
        act0 <- which(acti@breaks$phases[, 'tp0'] %within% tint)

        for (ph in act0) {
          dtseq <- seq(acti@breaks$phases[ph, 'tp0'], 
                       acti@breaks$phases[ph, 'tp1'],
                       by = 'mins')
          st <- acti@breaks$phases$State[ph]
          rug(dtseq, col=col.state[st], side=3, ticksize=-0.05)
        }
      }
      
      tdiff <- difftime(max(t.POSIX), min(t.POSIX), unit = 'sec')
      
      if (tdiff <= 3600) {
        ati <- seq(round(min(t.POSIX), 'secs'), max(t.POSIX), by = 'mins')
        
        while (length(ati) > 30) {
          ati <- ati[seq(1, length(ati), by=2)]
        }
        
        axis.POSIX(1, at=ati, format=("%H:%M:%S"), las=2)
      }
      else {
        ati <- seq(round(min(t.POSIX), 'hours'), max(t.POSIX), by = 'hours')
        
        while (length(ati) > 30) {
          ati <- ati[seq(1, length(ati), by=2)]
        }
        
        axis.POSIXct(1, at=ati, format=("%H:%M"), las=2)
      }
      
      if (addDiurnal) {
        coord <- cluster@Centroid[clustid,]
        if(!compareCRS(coord, CRS('+proj=longlat +datum=WGS84')))
          coord <- spTransform(coord, CRS('+proj=longlat +datum=WGS84'))
        
        dn <- getDiurnal(coord, startDate, endDate, tzone=localTZ)
        
        nigh <- which(dn$Night == 1)
        nnigh <- which(dn$Night == 0)
        
        while (length(nigh) > 0) {
          nnigh <- nnigh[nnigh > nigh[1]]
          
          if (length(nnigh) > 0) {
            nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nnigh[1]-1])
            nigh <- nigh[nigh > nnigh[1]-1]
          }
          else {
            nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nigh[length(nigh)]])
            nigh = NULL
          }
          rect(nout$ns, min(x), nout$ne, max(x), col = 'dark gray', border=NA)
        }
        
        da <- which(dn$Dawn == 1)
        dda <- which(dn$Dawn == 0)
        du <- which(dn$Dusk == 1)
        ddu <- which(dn$Dusk == 0)
        
        while (length(da) > 0) {
          dda <- dda[dda > da[1]]
          
          if (length(dda) > 0) {
            daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[dda[1]-1])
            da <- da[da > dda[1]-1]
          }
          
          else {
            daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[da[length(da)]])
            da <- NULL
          }
          rect(daout$das, min(x), daout$dae, max(x), col = 'light gray', border=NA)
        }
        
        while (length(du) > 0) {
          ddu <- ddu[ddu > du[1]]
          
          if (length(ddu) > 0) {
            duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[ddu[1]-1])
            du <- du[du > ddu[1]-1]
          }
          
          else {
            duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[du[length(du)]])
            du <- NULL
          }
          rect(duout$dus, min(x), duout$due, max(x), col = 'light gray', border=NA)
        }
      }
      
      ## Add Cluster point location
      abline(v = clust$dt, lwd = 2, col = 'black')
      if (pointLabel) {
        mtext(text = clust$PointID, las = 2, side = 3, outer = FALSE, at=clust$dt)
      }
      
      if (missing) {
        tint2 <- interval(min(clust$dt), max(clust$dt))
        allFixes <- cluster@OrigData$dt[which(cluster@OrigData$dt %within% tint2)]
        fixint <- min(round(diff(allFixes)))
        
        seqt <- seq(min(clust$dt), max(clust$dt), fixint)
        
        missFix <- c()
        for (st in 1:length(seqt)) {
          if (any(abs(difftime(seqt[st], allFixes, units='secs')) <= 60*5)) next
          else missFix <- c(missFix, seqt[st])
        }
        
        if (!is.null(missFix)) {
          missFix <- as.POSIXct(missFix, origin=lubridate::origin, tz = tz(cluster@OrigData$dt))
          abline(v = missFix, lwd = 2, col = 'red')
        }
      }
      
      lines(t.POSIX, x, col = "grey")
      
      ws <- windowsweep$ws
      
      if (type == "smooth") {
        if ("pp.smooth" %in% names(acti@bcpa))
          pp <- acti@bcpa$pp.smooth[inds,]
        else pp <- PartitionParameters(acti@bcpa, type = "smooth")[inds,]
        GoodBreaks <- ws$Break[ws$Model > 0]
        GoodBreaks <- as.data.frame(table(GoodBreaks))
        GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[, 1], t)])
        GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[, 1]))
        GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold, ]
        if (!suppressBreakLine)
          abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold * 2, col = col.cp)
        rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
        rho.int <- round(rho.scaled * 999 + 1)
        palette(topo.colors(1000))
        points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
               cex = pt.cex, lwd = 0.5)
        lines(t.POSIX, pp$mu.hat, lwd = 1.5, col = col.mean)
        lines(t.POSIX, pp$mu.hat + pp$s.hat, col = col.sd, lwd = 1.5)
        lines(t.POSIX, pp$mu.hat - pp$s.hat, col = col.sd, lwd = 1.5)
        rho.hat <- pp$rho.hat
      }
      
      if (type == 'flat')
        stop('Currently, only the smooth bcpa is permitted with cluster data')
      
      if (legend) {
        legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
        legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE),
                                  length = 5), 2)
        if (rho.where != "nowhere")
          legend(rho.where, bg = "white", legend = c(expression(hat(rho)),
                                                     legend.rhos), pch = 19, ncol = 3, col = c(0,
                                                                                               legend.cols), xjust = 0.5, yjust = 0.5)
        if (mu.where != "nowhere")
          legend(mu.where, bg = "white", legend = c(expression(hat(mu)),
                                                    expression(hat(mu) %+-% hat(sigma))), lty = 1,
                 lwd = 2:1, col = c("black", "red"), xjust = 0.5,
                 yjust = 0.5)
      }
      par(old.par)
    }
    palette("default")    
  }
})

setMethod('plot', 'BCPAtraj', function(x, startDate=min(x@bcpa$t.POSIX), 
                                       endDate=max(x@bcpa$t.POSIX), 
                                       suppressBreakLine = T,
                                       type = c("smooth", "flat")[1], 
                                       threshold = 3, clusterwidth = 1, 
                                       proj4string,
                                       addDiurnal = T, localTZ = Sys.timezone(),
                                       col.cp = rgb(0.5, 0, 0.5, 0.5), pt.cex = 0.5, 
                                       legend = TRUE,
                                       rho.where = "topleft", mu.where = "nowhere", 
                                       col.sd = "red", col.mean = "black", ...)   
{
  ###########################################################################
  ################### WRAPPER FUNCTION FOR PLOT.BCPA() in package bcpa ######
  ###########################################################################
  activity <- x
  if (class(activity) != 'BCPAact' & class(activity) != 'BCPAtraj')
    stop("Activity must be of class 'BCPAtraj or BCPAact'")
  
  windowsweep <- activity@bcpa
  
  ## NEED TO IMPLEMENT ID checks
  
  tint <- interval(startDate, endDate)
  
  t.POSIXt <- windowsweep$t.POSIX
  inds <- which(t.POSIXt %within% tint)
  
  x <- windowsweep$x[inds]
  t.POSIX <- t.POSIXt[inds]
  t <- windowsweep$t[inds]
  
  plot(x ~ t.POSIX, type = "n", ...)
  
  if (addDiurnal) {
    coord <- activity@traj$Z.start[1]
    coord <- data.frame(x=Re(coord), y=Im(coord))
    coordinates(coord) <- coord
    proj4string(coord) <- proj4string
    
    if(!compareCRS(coord, CRS('+proj=longlat +datum=WGS84')))
      coord <- spTransform(coord, CRS('+proj=longlat +datum=WGS84'))
    
    dn <- getDiurnal(coord, startDate, endDate, tzone=localTZ)
    
    nigh <- which(dn$Night == 1)
    nnigh <- which(dn$Night == 0)
    
    while (length(nigh) > 0) {
      nnigh <- nnigh[nnigh > nigh[1]]
      
      if (length(nnigh) > 0) {
        nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nnigh[1]-1])
        nigh <- nigh[nigh > nnigh[1]-1]
      }
      else {
        nout <- data.frame(ns = dn$dt[nigh[1]], ne = dn$dt[nigh[length(nigh)]])
        nigh = NULL
      }
      rect(nout$ns, min(x), nout$ne, max(x), col = 'dark gray', border=NA)
    }
    
    da <- which(dn$Dawn == 1)
    dda <- which(dn$Dawn == 0)
    du <- which(dn$Dusk == 1)
    ddu <- which(dn$Dusk == 0)
    
    while (length(da) > 0) {
      dda <- dda[dda > da[1]]
      
      if (length(dda) > 0) {
        daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[dda[1]-1])
        da <- da[da > dda[1]-1]
      }
      
      else {
        daout <- data.frame(das = dn$dt[da[1]], dae = dn$dt[da[length(da)]])
        da <- NULL
      }
      rect(daout$das, min(x), daout$dae, max(x), col = 'light gray', border=NA)
    }
    
    while (length(du) > 0) {
      ddu <- ddu[ddu > du[1]]
      
      if (length(ddu) > 0) {
        duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[ddu[1]-1])
        du <- du[du > ddu[1]-1]
      }
      
      else {
        duout <- data.frame(dus = dn$dt[du[1]], due = dn$dt[du[length(du)]])
        du <- NULL
      }
      rect(duout$dus, min(x), duout$due, max(x), col = 'light gray', border=NA)
    }
  }
  
  lines(t.POSIX, x, col = "blue")
  
  ws <- windowsweep$ws
  
  if (type == "smooth") {
    if ("pp.smooth" %in% names(windowsweep))
      pp <- windowsweep$pp.smooth[inds,]
    else pp <- PartitionParameters(windowsweep, type = "smooth")
    GoodBreaks <- ws$Break[ws$Model > 0]
    GoodBreaks <- as.data.frame(table(GoodBreaks))
    GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[, 1], windowsweep$t)])
    GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[, 1]))
    GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold, ]
    if (!suppressBreakLine)
      abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold * 2, col = col.cp)
    rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
    rho.int <- round(rho.scaled * 999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
           cex = pt.cex, lwd = 0.5)
    lines(t.POSIX, pp$mu.hat, lwd = 1.5, col = col.mean)
    lines(t.POSIX, pp$mu.hat + pp$s.hat, col = col.sd, lwd = 1.5)
    lines(t.POSIX, pp$mu.hat - pp$s.hat, col = col.sd, lwd = 1.5)
    rho.hat <- pp$rho.hat
  }
  
  if (type == "flat") {
    cpsummary <- ChangePointSummary(windowsweep, clusterwidth = clusterwidth)
    phases <- cpsummary$phases
    breaks <- cpsummary$breaks
    whichphase <- findInterval(t, phases$t0)
    rho.hat <- phases$rho.hat[whichphase]
    rho.int <- round(rho.hat/max(rho.hat, na.rm = TRUE) *
                       999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
           cex = pt.cex, lwd = 0.5)
    closematch <- rep(NA, length = nrow(phases))
    for (i in 1:nrow(phases)) closematch[i] <- which(abs(t -
                                                           phases$t0[i]) == min(abs(t - phases$t0[i])))[1]
    phases$t0.POSIX <- t.POSIX[closematch]
    phases$t1.POSIX <- t.POSIX[c(closematch[-1], length(t))]
    t.mid <- (windowsweep$t[-1] + windowsweep$t[-length(windowsweep$t)])/2
    segments(phases$t0.POSIX, phases$mu.hat, phases$t1.POSIX,
             phases$mu.hat, lwd = 1.5, col = col.mean)
    segments(phases$t0.POSIX, phases$mu.hat - phases$s.hat,
             phases$t1.POSIX, phases$mu.hat - phases$s.hat, col = col.sd,
             lwd = 1.5)
    segments(phases$t0.POSIX, phases$mu.hat + phases$s.hat,
             phases$t1.POSIX, phases$mu.hat + phases$s.hat, col = col.sd,
             lwd = 1.5)
    abline(v = phases$t0.POSIX[-1], lwd = breaks$size/max(breaks$size) *
             4, col = col.cp)
  }
  
  if (legend) {
    legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
    legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE),
                              length = 5), 2)
    if (rho.where != "nowhere")
      legend(rho.where, bg = "white", legend = c(expression(hat(rho)),
                                                 legend.rhos), pch = 19, ncol = 3, col = c(0,
                                                                                           legend.cols), xjust = 0.5, yjust = 0.5)
    if (mu.where != "nowhere")
      legend(mu.where, bg = "white", legend = c(expression(hat(mu)),
                                                expression(hat(mu) %+-% hat(sigma))), lty = 1,
             lwd = 2:1, col = c("black", "red"), xjust = 0.5,
             yjust = 0.5)
  }
  par(ask = T)
  palette("default")
  par(ask = FALSE)
})

setMethod('plot', 'BCPAact', function(x, startDate=min(x@bcpa$t.POSIX), 
                                      endDate=max(x@bcpa$t.POSIX), 
                                      suppressBreakLine = T,
                                      type = c("smooth", "flat")[1], 
                                      threshold = 3, clusterwidth = 1, 
                                      addDiurnal = T,
                                      col.cp = rgb(0.5, 0, 0.5, 0.5), pt.cex = 0.5, 
                                      legend = TRUE,
                                      rho.where = "topleft", mu.where = "nowhere", 
                                      col.sd = "red", col.mean = "black", ...)   
{
  ###########################################################################
  ################### WRAPPER FUNCTION FOR PLOT.BCPA() in package bcpa ######
  ###########################################################################
  activity <- x
  if (class(activity) != 'BCPAact' & class(activity) != 'BCPAtraj')
    stop("Activity must be of class 'BCPAtraj or BCPAact'")
  
  windowsweep <- activity@bcpa
  
  ## NEED TO IMPLEMENT ID checks
  
  tint <- interval(startDate, endDate)
  
  t.POSIXt <- windowsweep$t.POSIX
  inds <- which(t.POSIXt %within% tint)
  
  x <- windowsweep$x[inds]
  t.POSIX <- t.POSIXt[inds]
  t <- windowsweep$t[inds]
  
  plot(x ~ t.POSIX, type = "n", ...)
  
  if (addDiurnal) {
    t.posix <- activity@OrigData[which(activity@OrigData$dt %within% tint),]
    if('Day' %in% names(activity@OrigData)) {
      nigh <- which(t.posix$Night == 1)
      nnigh <- which(t.posix$Night == 0)
      
      while (length(nigh) > 0) {
        nnigh <- nnigh[nnigh > nigh[1]]
        
        if (length(nnigh) > 0) {
          nout <- data.frame(ns = t.posix$dt[nigh[1]], ne = t.posix$dt[nnigh[1]-1])
          nigh <- nigh[nigh > nnigh[1]-1]
        }
        else {
          nout <- data.frame(ns = t.posix$dt[nigh[1]], ne = t.posix$dt[nigh[length(nigh)]])
          nigh = NULL
        }
        rect(nout$ns, min(x), nout$ne, max(x), col = 'dark gray', border=NA)
      }
      
      if('Dawn' %in% names(activity@OrigData)) {
        da <- which(t.posix$Dawn == 1)
        dda <- which(t.posix$Dawn == 0)
        du <- which(t.posix$Dusk == 1)
        ddu <- which(t.posix$Dusk == 0)
        
        while (length(da) > 0) {
          dda <- dda[dda > da[1]]
          
          if (length(dda) > 0) {
            daout <- data.frame(das = t.posix$dt[da[1]], dae = t.posix$dt[dda[1]-1])
            da <- da[da > dda[1]-1]
          }
          
          else {
            daout <- data.frame(das = t.posix$dt[da[1]], dae = t.posix$dt[da[length(da)]])
            da <- NULL
          }
          rect(daout$das, min(x), daout$dae, max(x), col = 'light gray', border=NA)
        }
        
        while (length(du) > 0) {
          ddu <- ddu[ddu > du[1]]
          
          if (length(ddu) > 0) {
            duout <- data.frame(dus = t.posix$dt[du[1]], due = t.posix$dt[ddu[1]-1])
            du <- du[du > ddu[1]-1]
          }
          
          else {
            duout <- data.frame(dus = t.posix$dt[du[1]], due = t.posix$dt[ddu[length(ddu)]])
            du <- NULL
          }
          rect(duout$dus, min(x), duout$due, max(x), col = 'light gray', border=NA)
        }
      }
    }
    else {
      warning('Diurnal patterns not shown...please run addDiurnal()')
    }
  }
  
  lines(t.POSIX, x, col = "blue")
  
  ws <- windowsweep$ws
  
  if (type == "smooth") {
    if ("pp.smooth" %in% names(windowsweep))
      pp <- windowsweep$pp.smooth[inds,]
    else pp <- PartitionParameters(windowsweep, type = "smooth")
    GoodBreaks <- ws$Break[ws$Model > 0]
    GoodBreaks <- as.data.frame(table(GoodBreaks))
    GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[, 1], windowsweep$t)])
    GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[, 1]))
    GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold, ]
    if (!suppressBreakLine)
      abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold * 2, col = col.cp)
    rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
    rho.int <- round(rho.scaled * 999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
           cex = pt.cex, lwd = 0.5)
    lines(t.POSIX, pp$mu.hat, lwd = 1.5, col = col.mean)
    lines(t.POSIX, pp$mu.hat + pp$s.hat, col = col.sd, lwd = 1.5)
    lines(t.POSIX, pp$mu.hat - pp$s.hat, col = col.sd, lwd = 1.5)
    rho.hat <- pp$rho.hat
  }
  
  if (type == "flat") {
    cpsummary <- ChangePointSummary(windowsweep, clusterwidth = clusterwidth)
    phases <- cpsummary$phases
    breaks <- cpsummary$breaks
    whichphase <- findInterval(t, phases$t0)
    rho.hat <- phases$rho.hat[whichphase]
    rho.int <- round(rho.hat/max(rho.hat, na.rm = TRUE) *
                       999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int,
           cex = pt.cex, lwd = 0.5)
    closematch <- rep(NA, length = nrow(phases))
    for (i in 1:nrow(phases)) closematch[i] <- which(abs(t -
                                                           phases$t0[i]) == min(abs(t - phases$t0[i])))[1]
    phases$t0.POSIX <- t.POSIX[closematch]
    phases$t1.POSIX <- t.POSIX[c(closematch[-1], length(t))]
    t.mid <- (windowsweep$t[-1] + windowsweep$t[-length(windowsweep$t)])/2
    segments(phases$t0.POSIX, phases$mu.hat, phases$t1.POSIX,
             phases$mu.hat, lwd = 1.5, col = col.mean)
    segments(phases$t0.POSIX, phases$mu.hat - phases$s.hat,
             phases$t1.POSIX, phases$mu.hat - phases$s.hat, col = col.sd,
             lwd = 1.5)
    segments(phases$t0.POSIX, phases$mu.hat + phases$s.hat,
             phases$t1.POSIX, phases$mu.hat + phases$s.hat, col = col.sd,
             lwd = 1.5)
    abline(v = phases$t0.POSIX[-1], lwd = breaks$size/max(breaks$size) *
             4, col = col.cp)
  }
  
  if (legend) {
    legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
    legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE),
                              length = 5), 2)
    if (rho.where != "nowhere")
      legend(rho.where, bg = "white", legend = c(expression(hat(rho)),
                                                 legend.rhos), pch = 19, ncol = 3, col = c(0,
                                                                                           legend.cols), xjust = 0.5, yjust = 0.5)
    if (mu.where != "nowhere")
      legend(mu.where, bg = "white", legend = c(expression(hat(mu)),
                                                expression(hat(mu) %+-% hat(sigma))), lty = 1,
             lwd = 2:1, col = c("black", "red"), xjust = 0.5,
             yjust = 0.5)
  }
  par(ask = T)
  palette("default")
  par(ask = FALSE)
})
###########################################################

#### Deprecate
plotActivity <- function(aDat, date, time, actVar, int = 'Monthly', aggInt = 'Hourly') {
  aDat$datetime <- dmy_hms(paste(aDat[, date], aDat[, time]))
  aDat$Time <- hms(aDat[, time])
  aDat$month <- month(aDat$datetime)
  aDat$week <- week(aDat$datetime)
  aDat$hour <- hour(aDat$datetime)
  aDat$minute <- minute(aDat$datetime)
  
  if (aggInt == 'Hourly') {
    if (int =='Monthly') {
      iDat <- aggregate(list(X = aDat[, actVar]),
                        list(interval = aDat$month, hour = aDat$hour), mean)
      sdDat <- aggregate(list(sdX = aDat[, actVar]),
                         list(interval = aDat$month, hour = aDat$hour), sd)
      iDat$se <- sdDat$sdX / sqrt(iDat$X)
      iDat$Lower <- iDat$X - (2*iDat$se)
      iDat$Upper <- iDat$X + (2*iDat$se)
    }
    if (int =='Weekly') {
      iDat <- aggregate(list(X = aDat[, actVar]),
                        list(interval = aDat$week, hour = aDat$hour), mean)
      sdDat <- aggregate(list(sdX = aDat[, actVar]),
                         list(interval = aDat$week, hour = aDat$hour), sd)
      iDat$se <- sdDat$sdX / sqrt(iDat$X)
      iDat$Lower <- iDat$X - (2*iDat$se)
      iDat$Upper <- iDat$X + (2*iDat$se)
    }
  }
  else {print('Only aggregating intervals of an hour are currently implemented, please assign aggInt = "Hourly"')}
  
  #bDat <- melt(aDat, id.vars=c('month', 'hour'), measure.vars=c('VeDBA'))
  #ggplot(bDat, aes(x=as.factor(hour), y=value)) +
  #  geom_boxplot(outlier.colour="black", outlier.shape=19, notch=T,
  #               outlier.size=1.5) + facet_grid(~month)
  
  
  p <- ggplot(iDat)+
    geom_ribbon(aes(x=hour, ymin=Lower, ymax=Upper),fill='blue')+
    geom_line(aes(x=hour,y=X),colour='black') +
    facet_grid(~interval) +
    #ggtitle(paste(aggInt, ' Activity')) +
    #ggtitle('Accelerometer Activity') +
    #xlab(aggInt) +
    xlab('Hourly Activity by Month') +
    #ylab(actVar) +
    ylab('Velocity') +
    scale_x_continuous(breaks=seq(0, max(iDat$hour), 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title=element_text(size=16,face="bold"))
  
  
  return(p)
}

pullActivity <- function(cDat, aDat, byCluster=TRUE, aVar='All', t_int=1440, 
                         inPar = T, nCores = detectCores()) {
  
  ## Test class for cDat and aDat
  if (class(cDat) != 'ClusterUN') 
    stop('Cluster object must be of class "ClusterUN"')
  if (class(aDat) != 'ClusterACT')
    stop('Activity object must be of class "ClusterACT"')
  if (cDat@Params$ID != aDat@id)
    stop('Cluster and activity data must come from identical individuals (IDs).')
  
  if (!byCluster) {
    if (aVar=='All') out <- cbind(aDat@dt, aDat@act)
    else out <- cbind(aDat@dt, aDat@act[, aVar])
  }
  
  else {
    if (class(aDat@dt)[1] == 'data.frame') {
      pullACT <- function(cluster, aDat, aVar, t_int) {
        inter <- interval(min(cluster$dt) - (t_int*60), max(cluster$dt) + (t_int*60))
        indx <- which(aDat@dt[,'dt'] %within% inter)
        
        if (aVar=='All') act <- cbind(aDat@dt[indx, ], aDat@act[indx, ])
        else act <- cbind(aDat@dt[indx, ], aDat@act[indx, aVar])
        
        return(act)
      }}
    
    else {
      pullACT <- function(cluster, aDat, aVar, t_int) {
        inter <- interval(min(cluster$dt) - (t_int*60), max(cluster$dt) + (t_int*60))
        indx <- which(aDat@dt %within% inter)
        
        if (aVar=='All') act <- cbind(dt = aDat@dt[indx], aDat@act[indx, ])
        else act <- cbind(dt = aDat@dt[indx], aDat@act[indx, aVar])
        
        return(act)
      }}
    
    if (inPar==T) {
      cl <- makeCluster(nCores)
      clusterExport(cl, ls(envir=globalenv()))
      clusterEvalQ(cl, library(lubridate))
      out <- parLapply(cl, cDat@Data, function(x) pullACT(x, aDat, aVar, t_int))
      stopCluster(cl)
    }
    
    else {
      out <- list()
      for (clust in cDat@Data) {
        o <- pullACT(clust, aDat, aVar, t_int)
        out <- append(out, list(o))
      }
    }
  }
  cDat@rawACT <- out
  return(cDat)
}
