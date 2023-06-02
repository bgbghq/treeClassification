##################################################################################################################
# Program   treeSpeciesClassification.R
# Purpose   Classify trees by species using NEON lidar, hyperspectral, and field data
# Person    Andy Whelan
# Date      February 20, 2023
# Modified  February 20, 2023
################################################################################################################


#---> setup
# Libraries
packages = c("httr", "jsonlite","ggplot2","neonUtilities","neonOS","devtools","sf","terra","lidR","data.table",
             "hsdar")
for(package in packages) {
  if(suppressMessages(!require(package,character.only = T))){
    install.packages(package)
    suppressMessages(library(package,character.only = T))
  }
}

devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
library(geoNEON)
library(rhdf5)

options(stringsAsFactors = F)



#--------------- Setup data -----------------------------------------------------------------
# load veg structure data
vst <- loadByProduct(dpID="DP1.10098.001", 
                     site="JERC",
                     check.size=F)




# Get available lidar data for the JERC 
meta = content(GET("http://data.neonscience.org/api/v0/products/DP1.30003.001"), as="text")
avail = fromJSON(meta, simplifyDataFrame = T, flatten=T)
avail$data$siteCodes$siteCode
avail$data$siteCodes$availableMonths[[21]] # 21st is JERC

# Get available hyperspectral data for the JERC 
hypMeta = content(GET("http://data.neonscience.org/api/v0/products/DP3.30006.001"), as="text")
hypAvail = fromJSON(hypMeta, simplifyDataFrame = T, flatten=T)
hypAvail$data$siteCodes$siteCode
hypAvail$data$siteCodes$availableMonths[[21]] # 21st is JERC


# start by subsetting data to plots with trees
vst.trees <- vst$vst_perplotperyear[which(
  vst$vst_perplotperyear$treesPresent=="Y"),] # plots with trees
vstData = vst$vst_apparentindividual # dbhs and heights

vst.loc <- getLocTOS(data=vst$vst_mappingandtagging,
                     dataProd="vst_mappingandtagging") # locations

vstData = vstData[grep("2021", vstData$date),] # latest data
vstData = vstData[vstData$plotID %in% vst.trees$plotID,]
vst.loc = vst.loc[vst.loc$plotID %in% vst.trees$plotID,] 
vst.loc1 = vst.loc[!is.na(vst.loc$adjEasting),] # just get the mapped trees
vstData = vstData[vstData$individualID %in% vst.loc1$individualID,]
vstData = vstData[order(vstData$date),]
length(unique(vstData$individualID))
length(unique(vstData$individualID[1:147]))
vstData = vstData[1:147,] # the latest of the latest measurements for trees. No duplicate individualIDs
dim(vstData)
vst.loc1 = vst.loc1[vst.loc1$individualID %in% vstData$individualID,]
vst.loc1 = vst.loc1[order(vst.loc1$individualID),]
length(unique(vst.loc$individualID))
vstData = vstData[order(vstData$individualID),]
vstData$individualID == vst.loc1$individualID
vstData$adjEasting = vst.loc1$adjEasting
vstData$adjNorthing = vst.loc1$adjNorthing
vstData$taxonID = vst.loc1$taxonID

# make an sf object out of vstData
trees = st_as_sf(vstData, coords=c("adjEasting","adjNorthing"), crs=32616)
table(trees$taxonID)
ext = st_bbox(trees) # get the extent
extRound = lapply(ext, function(x) round(x, digits=-3)) # For the lidar tiles, next 3 lines.
eastings = seq(extRound[[1]], extRound[[3]], 1000)
northings = seq(extRound[[2]], extRound[[4]], 1000)

#---> Get the lidar data
for(easting in eastings) { # DP3.30003.001 is for classified lidar
  for(northing in northings) {
    byTileAOP("DP1.30003.001", site=c("JERC"), year="2021", check.size=F,
              easting=easting, northing=northing)
  }
}

#---> Get the hyperspectral data
for(easting in eastings) { # DP3.30006.001 is for hyperspectral data
  for(northing in northings) {
    byTileAOP("DP3.30006.001", site=c("JERC"), year="2021", check.size=F,
              easting=easting, northing=northing)
  }
}



#--------------------- Get crown polygons with species -----------------------------------------------------
lidarPath = "DP1.30003.001/neon-aop-products/2021/FullSite/D03/2021_JERC_6/L1/DiscreteLidar/ClassifiedPointCloud/"
ctg = readLAScatalog(lidarPath)

# function to adjust window size based on tree height
f = function(x) {x * 0.1 + 3}

# function for tree segmentation in catalog_apply
treeSeg = function(chunk) {
  las = readLAS(chunk)
  if (is.empty(las)) return(NULL) 
  las1 = normalize_height(las, knnidw(k=8, p=2))
  las2 = filter_poi(las1, Classification != 2 & Z>-5 & Z<35)
  
  # chm
  chm_p2r_05 <- rasterize_canopy(las2, 0.5, p2r(subcircle = 0.2), pkg = "terra")
  kernel <- matrix(1,3,3)
  chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
  
  # find tree tops
  ttops = locate_trees(chm_p2r_05_smoothed, lmf(f))
  
  # segment
  algo = dalponte2016(chm_p2r_05_smoothed, ttops)
  las3 = segment_trees(las2, algo)
  crowns = tryCatch({crown_metrics(las3, func=.stdtreemetrics, geom="convex")},
           error=function(e) {
             message("Crap!")
             print(e)
             return(NA)
           })
  return(crowns)
}


opt_chunk_buffer(ctg) = 15
opt_chunk_size(ctg) = 200
opt    <- list(need_buffer = TRUE,   # catalog_apply will throw an error if buffer = 0
               automerge   = TRUE)   # catalog_apply will merge the outputs into a single obj
crowns = catalog_apply(ctg, treeSeg)
crowns = crowns[!is.na(crowns)]
first = crowns[[1]]
for(i in 2:length(crowns)) {
  first = rbind(first, crowns[[i]])
}
allCrowns = first # just rename it.
allCrowns$crownID = seq(1, nrow(allCrowns), 1) # kind of a new crownID

# intersect trees with crown segmentation
t = st_intersects(allCrowns, trees) # find which trees are located under which crowns
t1 = st_intersection(trees, allCrowns) # merge tree info with crown info
t1 = t1[!is.na(t1$height),]
t1 = as.data.table(t1) # make it a data table so the next part will work.
t2 = t1[t1[, .I[height == max(height)], by=crownID]$V1] # data from the tallest tree beneath each crown
c1 = allCrowns[lengths(t) > 0,] # Subset to crowns with measured trees beneath
c2 = c1[c1$treeID %in% t2$treeID,] # Subset again to crowns with trees with measured heights beneath
c2 = c2[order(t2$treeID),]
c2$taxonID = t2$taxonID # add taxonID
st_write(c2, "crowns.shp")
c2 = st_read("crowns.shp")
c3 = vect(c2)


# ----------------------- Hyperspectral data -----------------------------------------------------
# mask hyperspectral data
hypPath = "../firstTry/DP3.30006.001/neon-aop-products/2021/FullSite/D03/2021_JERC_6/L3/Spectrometer/Reflectance/"
hypList = list.files(hypPath)
f = paste0(hypPath, hypList[[1]])
myEPSG = h5read(f, "/JERC/Reflectance/Metadata/Coordinate_System/EPSG Code")
refInfo = h5readAttributes(f, "/JERC/Reflectance/Reflectance_Data")
myCRS = paste0("+init=epsg:", myEPSG)
myNoDataValue = as.numeric(refInfo$Data_Ignore_Value)


#---> Make a raster stack (pile) out of all the spectral bands
band2Rast = function(file, band, noDataValue, extent, CRS){
  out = h5read(file, "/JERC/Reflectance/Reflectance_Data", index=list(band,NULL,NULL))
  out = aperm(out, c(3,2,1)) # reshape the array to ('rows','columns','bands') or something like that.
  out[out == noDataValue] = NA
  outr = rast(out,crs=CRS)
  ext(outr) = extent
  return(outr)
}


tiles = list()
for(tile in hypList) {
  f = paste0(hypPath, tile)
  refInfo = h5readAttributes(f, "/JERC/Reflectance/Reflectance_Data")
  
  # Grab the UTM coordinates of the spatial extent
  xMin <- refInfo$Spatial_Extent_meters[1]
  xMax <- refInfo$Spatial_Extent_meters[2]
  yMin <- refInfo$Spatial_Extent_meters[3]
  yMax <- refInfo$Spatial_Extent_meters[4]
  
  # define the extent (left, right, top, bottom)
  rasExt <- ext(xMin,xMax,yMin,yMax)
  
  # make a list of the band numbers
  bList = seq(1, 426, 1)
  
  spat = band2Rast(file=f, band=bList, noDataValue=myNoDataValue, extent=rasExt, CRS=myCRS)
  bandNames = paste("Band_",unlist(bList),sep="")
  names(spat) = bandNames
  tiles[length(tiles) + 1] = mask(spat, c3)
}

m1 = tiles[[1]]
for(i in 2:length(tiles)) {
  m1 = merge(m1, tiles[[i]])
}
