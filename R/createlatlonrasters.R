  ##Create raster for the Latitude and Longitude and add to the covariates
  
  ## GET THE COORDINATES AS A DATAFRAME FROM ANY RASTER
  
  plot(covariates$slope)
  
  template <- covariates$slope
  
  plot(template)
  
  #template <- project(template, "epsg:4283")
  
  #crs(template, describe=T)
  
  XY <- xyFromCell(template, 1:ncell(template))
  
  # ##Create a covariate named lon
  lon <- template #Get a raster template
  lon[] <- XY[,"x"]
  names(lon) <- "lon"
  
  writeRaster(lon, "./data/covariates/lon.tif")
  
  plot(lon)
  
  # ##Create a covariate named lat
  lat <- template
  lat[] <- XY[,"y"]
  names(lat) <- "lat"
  
  writeRaster(lat, "./data/covariates/lat.tif")