#----------------
# These function do the preprocessing to create an ssn object from files created by TAUDEM 



SSNdat_fromTAUDEM_tree <- function(tree, save.dir = NULL){
  
  #####Function uses a tree.txt file from TAUDEM and wites to files
  #neti.dat the link ID and the SSN binary representation of that link
  
  # tree is a dataframe read from tree.txt
  #tree.txt contains the following information as a table
  # Link Number (Arbitrary - will vary depending on number of processes used)
  # Start Point Number in Network coordinates (*coord.dat) file (Indexed from 0)
  # End Point Number in Network coordinates (*coord.dat) file (Indexed from 0)
  # Next (Downstream) Link Number.  -1 indicates no links downstream, i.e. a terminal link
  # First Previous (Upstream) Link Number. Points to Link Number. -1 indicates no upstream links.
  # Second Previous (Upstream) Link Numbers. Points to Link Number. -1 indicates no upstream links. Where only one previous link is -1, it indicates an internal monitoring point where the reach is logically split, but the network does not bifurcate.
  # Strahler Order of Link
  # Monitoring point identifier at downstream end of link. -1 indicates downstream end is not a monitoring point
  # Network magnitude of the link, calculated as the number of upstream sources (following Shreve). 
  
  
  if (is.null(save.dir)) {save.dir <- "./"}
  pos.out <- which(tree$DSLINK == -1) #outlet edges
  
  for (i in 1:length(pos.out)){
    pos <- pos.out[i]
    rid <- c(tree$LINKID[pos])
    binary <- c("1")
    
    paths.x <- list()
    links.up1 <- tree$USLINK1[pos]
    links.up2 <- tree$USLINK2[pos]
    unprocessed.links <- c(links.up1, links.up2)
    paths.x[[paste(rid)]] <- binary
    if (links.up1 != -1 & links.up2 != -1){
  
    while(length(unprocessed.links)>0){
      
      start.link <- unprocessed.links[1]
      #print(start.link)
      unprocessed.links <- unprocessed.links[-1]
      
      pos.down.link <- which(tree$USLINK1 == start.link)
      if (length(pos.down.link) == 0) {
        pos.down.link <- which(tree$USLINK2 == start.link)
      }
      down.link <- tree$LINKID[pos.down.link]
      uplinks <- c(tree$USLINK1[pos.down.link],tree$USLINK2[pos.down.link])
      
      other.link <- uplinks[!(uplinks %in% start.link)]
      get.source.binary <- paths.x[[paste(down.link)]]
      other.link.processed <- sum(names(paths.x) %in% paste(other.link))
      if (other.link.processed == 0){
        newbinary <- paste0(get.source.binary,"0")
      } else {
        newbinary <- paste0(get.source.binary,"1")
      }
      pos.link <- which(tree$LINKID == start.link)
      uplinks <- c(tree$USLINK1[pos.link],tree$USLINK2[pos.link])
      uplinks <- uplinks[uplinks != -1]
      unprocessed.links <- c(unprocessed.links, uplinks)
      if (other.link != -1){
        pos.link <- which(tree$LINKID == other.link)
        uplinks <- c(tree$USLINK1[pos.link],tree$USLINK2[pos.link])
        uplinks <- uplinks[uplinks != -1]
        unprocessed.links <- c(unprocessed.links, uplinks)
      }
      unprocessed.links <- unique(unprocessed.links)
      paths.x[[paste(start.link)]] <- newbinary
      
    }
    } else {
      
    }
    datx <- data.frame(
      rid = names(paths.x), 
      binaryID = trimws( unlist(paths.x)))
    fid <- paste0(save.dir,"/netID",i,".dat")
    output.file <- file(fid, "wb")
    write.csv(datx,  file = output.file,
              row.names = FALSE, quote = FALSE
              )
    close(output.file)
  }
}


makeEdges <- function(demfile, areathreshold){
  #Function uses TauDem to process a dem raster to create a streamnetwork:
 
  
  #Input
  #demfile is the name of a geotif without .tif extension in the present working directory
  #area threshold, threshold contributing area defining streams (units of pixels)
  
  #Output in the working directory
  # a) a pit filled dem, demf.tif
  # b) a D8 flow direction raster, demp.tif
  # c) a slopes along D8 flow directions, demsd8.tif
  # d) contributing area along flow directions, demad8.tif
  # e) a stream network mask based on a threshold contributing area, areathreshold (pixels), src.tif
  # f) a stream order raster, ord.tif
  # g) text file documenting the stream network topology, tree.txt file (see TauDEM)
  # h) a stream network shapefile, edges.shp
  # i) a watershed raster, wshed.tif
  
  print("pit remove")
  sys <- paste0("mpiexec -n 4 PitRemove -z ",demfile,".tif -fel demf.tif")
  system(sys)
  
  print("flow directions")
  sys <- paste0("D8Flowdir -fel demf.tif -p demp.tif -sd8 demsd8.tif")
  system(sys)
  
  print("contributing area")
  sys <- paste0("mpiexec -n 3 AreaD8 -p demp.tif -nc -ad8 demad8.tif")
  system(sys)
  
  print("stream mask from contributing area threshold")
  sys <- paste0("mpiexec -n 3 Threshold -ssa demad8.tif -src src.tif -thresh ",areathreshold)
  system(sys)
  
  print("Stream network shapefile, washeshed raster and tree.txt")
  sys <- paste("mpiexec -n 3 StreamNet -p demp.tif -fel demf.tif -ad8 demad8.tif -src src.tif -ord ord.tif -tree tree.txt -coord coord.txt -net edges.shp -netlyr edges -w wshed.tif")
  system(sys)
}


add_netID_2shp <- function(edges, path = NULL, area = 1L){
  
  #Multiple stream networks can be created so this function adds the attribute
  #netID to the SpatialLinesDataFrame edges
  
  #Input
  #edges, of class spatialLinesDataFrame (object created by StreamNet in TauDEM)
  #path is the path to the edges shapefile
  #minarea is the area of a pixel in the dem
  #vALUE
  #Modified edges object with 
  #a net id assigend to separate river networks
  #a unique river segment id rid assinged to each link
  #the local area of the link
  
  if (is.null(path)){ path = "./"}
  edges@data$rid <- edges@data$LINKNO
  edges@data$netID <- NA
  fls <- list.files(path = path, pattern = "*.dat", full.names = TRUE)
  for (i in 1:length(fls)){
    net.rid <- read.csv(fls[i], colClasses = c("integer","character"))$rid
    pos <- which( edges@data$rid %in% net.rid)
    edges@data$netID[pos] <- i
  }
  edges@data$upDist <- edges@data$DOUTSTART
  edges@data$rid <- edges@data$LINKNO
  edges@data$area <- edges@data$DSContArea - shp@data$USContArea
  edges@data$area[edges@data$area == 0 ] <- area
  
  edges
}


dists <- function(crds){
   if (nrow(crds) > 2){
   dxs  <- apply(cbind(crds[-length(crds[1,]),],crds[-1,]), 
                 MARGIN=1,
        FUN = function(x){
          ((x[1]-x[3])^2+(x[2]-x[4])^2)^0.5
        })
   cumdist <- c(0, cumsum(dxs))
   } else {
     if (nrow(crds) == 2){
       
         cumdist <- c(0, ((crds[1] - crds[2])^2+(crds[3] - crds[4])^2)^0.5)
         if (cumdist[2] == 0){
           cumdist[2] = shp@data$Length[i]
         }
     } else {
       cumdist <- c(0,0)
     }
   }
   cumdist
}


makePredPoints <- function(edges, sites){
  #Creates a SpatialPointsDataFrame from streamnetwork SpatialLinesData frame
  #locID
  #rid
  #upDist
  #x
  #y
  pnts <- coordinates(edges)
  ns <- length(sites)
  newmat <- rep(NA, 1,nrow=1, ncol=6)
  for (i in 1:length(edges)){
    crds <- pnts[[i]][[1]]
    n <- dim(crds)[1]
    crds <- crds[n:1,]
    dstart <- shp@data$DOUTSTART[i]
    rid <- shp@data$rid[i]
    cumdist <- cumsum(c(0,sapply(2:dim(crds)[1], FUN = function(y) {
      sqrt((crds[y-1,1] - crds[y,1])^2+(crds[y-1,2]-crds[y,2])^2)
    })))
    updist <- dstart - cumdist
    nid <- edges@data$netID[i]
    ord <- edges@data$strmOrder[i]
    newmat <- rbind(newmat, cbind(cbind(crds, rep(rid, n)), updist, rep(nid, n), rep(ord, n)))
  }
  newmat <- newmat[-1,]
  colnames(newmat) <- c("x","y","rid","upDist","netID","strmOrder")
  
  
  df2 <- data.frame(newmat)
  df2$pid <-   ns + 1:dim(df2)[1]
  df2$locID <- ns + 1:dim(df2)[1]
  ##remove pred points with same coordinates as sites
  #may not be necessary
  #unique location id, match sites and predpoints
  #crdsites <- coordinates(sites)
  #for (i in 1:dim(crdsites)[1]){
  #  pos <- which(df2$x == crdsites[i,1] & df2$y == crdsites[i,2])
  #  df2$locID[pos] <- i
  #}
  #df2$locID[which(is.na(df2$locID))] <- 1:sum(is.na(df2$locID)) + ns
  
  df2$E = df2$x
  df2$N = df2$y
  coordinates(df2) <- ~x+y
  proj4string(df2) <- proj4string(shp)
  
  return(df2)
}


pointcount = function(r, pnts){
  #Create a raster from the raster with values equal to the number of points in pts in each cell of r
  
  #Input
  #r raster
  #pnts SpatialPointsDataFrame
  
  pts <- sp::coordinates(pnts)
  # make a raster of zeroes like the input
  r2 = r
  r2[!is.na(r[])] <- 0
  # get the cell index for each point and make a table:
  counts = table(raster::cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] <- counts
  return(r2)
}

accumDownStream <- function(demd8,weightaccum, weight = NULL){
  #accumulate along flow directions given by demd8 the weights given by weightaccum
  
  #Input 
  #demd8 character specifying name of flow direction raster i.e. "demp.tif"
  #weightaccum, character specifying name of raster to contain the accumulated values
  #weight optional character vector specifying raster which to accumulate, if weight is null the accumulated area in pixels is returned
  
  #value 
  #A tif file raster named by weightaccum is written to the working directory
  if (!is.null(weight)){
    sys <- paste0("mpiexec -n 3 AreaD8 -p ",weightaccum," -ad8 ",demd8," -wg ",weight," -nc ")
  } else {
    sys <- paste0("mpiexec -n 3 AreaD8 -p ",weightaccum," -ad8 ",demd8," -nc ")
  }
  system(sys)
}


net2keep <- function(filepath = paste0(getwd(),"/ssn"), prednames = "predPoints" ){
  #Sometimes TauDEM can create many stream networks
  #But to do statistics you need at least a few points on a network
  #This functionis designed to drop parts of edges, and predPoints which are not 
  #on a network with at least one sites value
  
  #Input
  #filepath a character giving the path to the ssn directory
  #prednaes: a character giving the name of the prediction points shape file.
  
  #Output
  #trimmed saved edges.shp, predPoints.shp in the ssn directory (or other name provided)
  #note the netIds are renumbered from 1 to the number of valid netIds
  #unused *.dat files deleted from the ssn directory
  
  sites <- readOGR(dsn = filepath, layer = "sites")
  edges <- readOGR(dsn = filepath, layer = "edges")
  predpoints <- readOGR(dsn = filepath, layer = prednames)
  
  neids <- sort(unique(sites@data$netID))
  allneids <- sort(unique(edges@data$netID))
  remneids <- allneids[!(allneids %in% neids)]
  
  sites <- sites[sites@data$netID %in% neids,]
  edges <- edges[edges@data$netID %in% neids,]
  predpoints <- predpoints[predpoints@data$netID %in% neids,]
  
  #dat.files <- list.files(path = filepath, pattern = "*.dat")
  if (length(remneids)>0){
    for (i in 1:length(remneids)){
      fid <- paste0(filepath, "/netID",remneids[i],".dat")
      file.remove(fid)
    }
  }
  
  
  nidmap <- 1:length(neids)
  dat.files <- list.files(path = filepath, pattern = "*.dat")
  for (i in nidmap){
    sites@data$netID[sites@data$netID == neids[i]] <-  i
    edges@data$netID[edges@data$netID == neids[i]] <-  i
    predpoints@data$netID[predpoints@data$netID == neids[i]] <-  i
    file.rename(paste0(filepath,"/",dat.files[i]), paste0(filepath, "/netID",i,".dat"))
    
  }
  
  
  
  writeOGR(sites, dsn = filepath, layer = "sites", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  writeOGR(edges, dsn = filepath, layer = "edges", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  writeOGR(predpoints, dsn = filepath, layer = prednames, driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
}


library(rgdal)
library(sp)

#Set the working directory
setwd(".\\Data")

#Process the dem and make the edges shapefile, and tree.txt
makeEdges("dem", areathreshold = 1000)

#If not already done create the ssn folder

#Read the tree.txt file and create the *.dat files
tree <- read.table("tree.txt",header=FALSE)
names(tree) <- c("LINKID","STRT.PNT","END.PNT","DSLINK","USLINK1","USLINK2","STRAHLER","MONITORING","SHREVE")
save.dir = paste0(getwd(),"/ssn")
SSNdat_fromTAUDEM_tree(tree, save.dir = save.dir)

#Read, edit and save the edges shapefile in the ssn folder
edges <- readOGR(dsn = "edges.shp")
edges <- add_netID_2shp(edges, path = paste0(getwd(),"/ssn"))
writeOGR(ecdges,dsn = paste0(getwd(),"/ssn"), layer = "edges",overwrite_layer = TRUE, 
         driver = "ESRI Shapefile")

#Read the observations shapefile and save in the ssn folder
sites <- readOGR(getwd(), layer = "apha2")
mcpa <- sites@data[,grep("MCPA", names(sites@data))]
maxmcpa <- apply(mcpa, MARGIN = 1, FUN = function(x) ifelse(is.finite(max(x, na.rm=TRUE)),max(x, na.rm=TRUE), 0)) 
sites@data$maxmcpa <- maxmcpa

#Create set of prediction points on the stream network
#This simple approach here gets points at each raster cells on the stream network (edges.shp) 
pps <- makePredPoints(edges, sites)
writeOGR(pps, dsn = paste0(getwd(),"/ssn"), layer = "predPoints", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Simple approach to snap sites to the nearest stream network node
#You could use TuaDEM to move points downslope to the nearest stream network cell instead or
#simply do in manally in a GIS (using sites.shp and src.tif
xy <- coordinates(pps)
csites <- coordinates(sites)
wichmin <- apply(csites, MARGIN = 1, FUN = function(x){
  mindists <- apply(xy, MARGIN = 1, FUN = function(y){
    min(sqrt((x[1] - y[1])^2+(x[2]-y[2])^2))
  })
  which.min(mindists)
})
sites@data$rid <- pps@data$rid[wichmin]
sites@data$upDist <- pps@data$upDist[wichmin]
sites@data$netID <- pps@data$netID[wichmin]
sites@data$strmOrder <- pps@data$strmOrd[wichmin]
sites@data$locID <- pps@data$locID[wichmin]
sites@data$elev <- raster::extract(raster::raster("dem.tif"), sites)
sites@data$pid = 1:length(sites@data[,1])
writeOGR(sites, dsn = paste0(getwd(),"/ssn"), layer = "sites", driver = "ESRI Shapefile", overwrite_layer = TRUE)



#Clean up and keep only parts of network with observations (sites) 
net2keep(filepath = paste0(getwd(),"/ssn"), prednames = "predPoints" )

#Read in ssn object using ssn package
library(SSN)
F01.ssn <- importSSN(filepath = paste0(getwd(),"/ssn"), predpts = "predPoints", o.write = TRUE)
createDistMat(F01.ssn, predpts = "predPoints", o.write = TRUE, amongpreds = TRUE)
F01.ssn <- importSSN(filepath = paste0(getwd(),"/ssn"), predpts = "predPoints", o.write = TRUE)
plot(F01.ssn, VariableName = "maxmcpa")

#Create an additive function value
F01.ssn <- additive.function(F01.ssn, "area", "computed.afv")

#x <- readRDS("D:\\GavanMcGrath_14062018\\Data\\AHPC\\SSN\\LoughForbes\\ssn\\distance\\obs\\dist.net3.RData")
#x2 <- readRDS("D:\\GavanMcGrath_14062018\\Data\\AHPC\\SSN\\MiddleFork04.ssn\\distance\\obs\\dist.net1.RData")


#Exponential.tailup and taildown work but for some reason euclidian productes errors
fit <- glmssn(maxmcpa ~ upDist, 
              ssn.object = F01.ssn, 
              family = "Gaussian", 
              CorModels = c("Exponential.tailup"),
              use.nugget = TRUE,
              use.anisotropy = FALSE,
              addfunccol = "computed.afv",
              #trialscol = NULL, 
              EstMeth = "ML",
              #useTailDownWeight = FALSE, trans.power = NULL,trans.shift = 0,
              control = list(max.range.factor = 4, 
                             trunc.pseudo = NULL,
                             maxiter.pseudo = 20, 
                             beta.converge = 1e-05)
)

F01.pred <- predict(fit, "predPoints")

summary(fit)
covparms(fit)
varcomp(fit)
plot(fit)
GR2(fit)


plot(ssnpred)
plot(ssnpred, VarPlot = "Predictions", breaktype = "quantile", lwd = 2)

#residuals
resid <- residuals(fit, cross.validation=TRUE)
plot(resid)
hist(resid, xlab = "Raw Residuals")
qqnorm(resid)

ESVF <- Torgegram(F01.ssn, "maxmcpa",nlagcutoff = 6,maxlag = 30000)
plot(ESVF, xlim = c(0,30000))


plot(ESVF, sp.relationship = "fu", min.cex = .4, max.cex = 4,
     main = "Flow-unconnected Torgegram")

plot(ESVF, min.cex = .4, max.cex = 4, col = c("darkgray", "black"),
     main = "", xlab = "Stream Distance (m)",xlim = c(0,50000), ylim = c(0,1))


ESVF <- Torgegram(resid, "_resid_", nlagcutoff = 4,maxlag = 30000)
plot(ESVF, xlim = c(0,20000))
