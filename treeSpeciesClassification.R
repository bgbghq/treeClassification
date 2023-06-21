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
             "stringr","Boruta","randomForest")
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
c2 = st_read("../firstTry/crowns.shp")
c3 = vect(c2)
crowns = c3

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
  if(is.character(file)) {
    out = h5read(file, "/JERC/Reflectance/Reflectance_Data", index=list(band,NULL,NULL))
  }else{
    out = file
  }
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

# make a list of the band numbers
bList = seq(1, 426, 1)
tileLons = as.numeric(unlist(lapply(hypList, function(x) substr(x, 19, 24))))
tileLats = as.numeric(unlist(lapply(hypList, function(x) substr(x, 26,32))))
tileCoords = data.frame(lon = tileLons, lat = tileLats)


#----------------------------------- Try to extract spectra without using so much memory -----------------
# function to get the spectra in the area of each crown
# "crowns" is an sf file with crown polygons. origin: lower left
# origin of NEON spectral data is upper left


bigMatrix = list()
for(i in 1:nrow(crowns)) { # for each crown
  
  # Figure out which hyperspectral tile to load and load it as a raster
  extTmp = ext(crowns[i,]) # bounding box of the crown
  
  lonTmp = floor(extTmp[1]/1000)*1000
  latTmp = floor(extTmp[3]/1000)*1000
  
  paths = c(hypList[grep(paste0(lonTmp,"_",latTmp), hypList)], lonTmp, latTmp)
  
  tTmp = H5Fopen(paste0(hypPath,paths[1]))
  t2Tmp = tTmp&'JERC/Reflectance/Reflectance_Data'
  
  tmpMatrix = list()
  for(j in 1:426) { # for each wavelength
    thing = t2Tmp[j,,] # [bands,columns,rows], origin: upper left
    thing[thing == -9999] = NA
    thing = aperm(thing, c(2,1))
    
    # create a matrix of pixel values
    thing2 = rast(thing, crs=myCRS, extent=c(lonTmp,lonTmp+1000,latTmp,latTmp+1000))
    thing3 = crop(thing2, crowns[i])
    thing3 = mask(thing3, crowns[i], touches=F)
    thing4 = as.vector(thing3)
    
    # Change the matrix into a vector
    tmpMatrix[[length(tmpMatrix)+1]] = as.numeric(na.omit(thing4))
  }
  
  # Close the dataset
  h5closeAll()
  
  # Put the vectors into a matrix as columns and add a column that specifies the tree species
  tmpMatrix = do.call(cbind, tmpMatrix)
  tmpDf = data.frame(tmpMatrix, taxonID=crowns[i]$taxonID)
  
  # Put the matrix into the bigMatrix list
  bigMatrix[[length(bigMatrix)+1]] = tmpDf
}
# rbind the matrices to make a longer matrix
outFile = do.call(rbind, bigMatrix)

write.csv(outFile, "SpecModelData.csv")
outFile = read.csv("SpecModelData.csv")
head(outFile)
outFile = outFile[,2:428]
dim(outFile)
boruta = Boruta(x=outFile[,1:426], y=as.factor(outFile[,427]), doTrace=0, pValue=0.05)

finalvars = getSelectedAttributes(boruta, withTentative = F)

(modelRF_class = randomForest(x=outFile[ , finalvars], y=as.factor(outFile$taxonID), # This get's the response 
                           # variable to be BA per 20m pixel
                           importance = TRUE))
    
  crownChunks = list()
    extTmp = ext(crowns[i,]) # bounding box of the crown
    
    lonTmp = floor(extTmp[1]/1000)*1000
    latTmp = floor(extTmp[3]/1000)*1000
    lonMin = floor(extTmp[1,]) # round down
    lonMax = ceiling(extTmp[2]) # round up
    latMin = floor(extTmp[3,]) # round down
    latMax = ceiling(extTmp[4]) # round up
    
    # Get paths to all tiles the crown is in (up to four)
    paths = list()
    paths[[length(paths)+1]] = c(hypList[grep(paste0(lonTmp,"_",latTmp), hypList)], lonTmp, latTmp)
    if(latMax-latTmp > 1000) paths[[length(paths)+1]] = 
      c(hypList[grep(paste0(lonTmp,"_",latTmp+1000),hypList)], lonTmp, latTmp+1000)
    if(lonMax-lonTmp > 1000) paths[[length(paths)+1]] = 
      c(hypList[grep(paste0(lonTmp+1000,"_",latTmp),hypList)], lonTmp+1000, latTmp)
    if(lonMax-lonTmp > 1000 & latMax-latTmp > 1000) {
      paths[[length(paths)+1]] = c(hypList[grep(paste0(lonTmp,"_",latTmp+1000),hypList)], lonTmp, latTmp+1000)
      paths[[length(paths)+1]] = c(hypList[grep(paste0(lonTmp+1000,"_",latTmp),hypList)], lonTmp+1000, latTmp)
      paths[[length(paths)+1]] = c(hypList[grep(paste0(lonTmp+1000,"_",latTmp+1000),hypList)], lonTmp+1000, latTmp+1000)
    }
    
    for(j in 1:length(paths)) {
      
      # find the relative coordinates of the crown to subset
      xmin = round(extTmp[1])-as.numeric(paths[[j]][2]) # x extents in unprojected hyperspectral array
      if(xmin < 0) xmin=1
      xmax = round(extTmp[2])-as.numeric(paths[[j]][2])
      if(xmax > 1000) xmax=1000
      ymin = round(extTmp[3])-as.numeric(paths[[j]][3]) # y extents in unprojected hyperspectral array
      if(ymin < 0) ymin=1
      ymax = round(extTmp[4])-as.numeric(paths[[j]][3])
      if(ymax > 1000) ymax=1000
      
      # Open the tile the crown starts in and mask it
      tTmp = H5Fopen(paste0(hypPath,paths[[j]][1]))
      t2Tmp = tTmp&'JERC/Reflectance/Reflectance_Data'
      thing = t2Tmp[,xmin:xmax,ymin:ymax] # [bands,columns,rows], origin: upper left
      thing[thing == -9999] = NA
      thing = aperm(thing, c(3,2,1))
      thing2 = rast(thing, crs=myCRS, extent=c(lonTmp,lonTmp+1000,latTmp,latTmp+1000))
      ext(thing2) = c(xmin,xmax,ymin,ymax)
      res(thing2) = c(1,1)
      crownChunks[[length(crownChunks)+1]] = thing2
    }
    crownChunk = do.call(merge, crownChunks)
    ext(crownChunk) = c(lonMin,lonMax,latMin,latMax)
  }
  
}
    
barf = 


    # Identify the path(s) to the hyperspectral data tile(s).
    if(extTmp[4]-latTmp > 1000) { # if the crown strays latitudinally into another tile
      paths = paste0(paths, ",",hypList[grep(paste0(lonTmp+1000,"_",latTmp), hypList)])
      side = 1
    }
    if(extTmp[2]-lonTmp > 1000) { # If the crown strays longitudinally into another tile
      paths = paste0(paths,",", hypList[grep(paste0(lonTmp,"_",latTmp+1000), hypList)])
      side = 2
    }
    if(extTmp[2]-lonTmp > 1000 & extTmp[4]-latTmp > 1000) { # If the crown strays into 3 other tiles
      paths = paste0(paths,",", hypList[grep(paste0(lonTmp+1000,"_",latTmp+1000), hypList)])
      side = 3
    }
    
    # if there are more than one path, split them
    paths = strsplit(paths, ",")
    
    
    
    tmpSpectra = list()
    for(j in 1:length(paths)) {
      
      # got to get the right part of each tile
      if(i == 1){
        tTmp2 = tTmp&'JERC/Reflectance/Reflectance_Data'[]
      }
      if(i == 2 & side == 1){
        tTmp2 = tTmp&'JERC/Reflectance/Reflectance_Data'[]
      }
      if(i == 2 & side == 2)
      bigTile[[length(bigTile)+1]] = tTmp2[]
    }
    
  }
  
      # subset the data
      t3 = t2[,xmin:xmax,ymin:ymax]
      
      # Get attributes (projection and crs)
      a1 = tile&'/JERC/Reflectance/Metadata/Coordinate_System/'
      tmpCRS = paste0("+init=epsg:", a1$`EPSG Code`)
      
      # make a spatial extent for the crown
      utmXmin = lonTmp+xmin
      utmXmax = lonTmp+xmax
      utmYmin = latTmp+ymin
      utmYmax = latTmp+ymax
      utmExt = ext(utmXMin, utmXmax, utmYmin, utmYmax)
      
      # make it into a raster
      spat = band2Rast(t3, bList, -9999, extent=utmExt,CRS=tmpCRS)
      bandNames = paste("Band_",unlist(bList),sep="")
      names(spat) = bandNames
      tiles[length(tiles) + 1] = mask(spat, c3)
      # close dataset 
      h5closeAll()

    }
    
    
    tile1 = H5Fopen()
    print(extTmp)
  }
}
crownSpectra(crowns)
crowns = read_sf("../firstTry/crowns.shp")
