packages <- c("purrr","tidyverse","sf","exactextractr","raster","fasterize")
purrr::walk(packages, library, character.only=T)

property <- st_read("E:\\01_NSW phase 2\\property_guras\\property_boundaries_GURAS.shp")


woody <- raster("E:\\01_NSW phase 2\\statelevelanalysis\\data\\woody\\nsw_state_woody_100m_majorityrule.tif")

template <- raster("E:\\01_NSW phase 2\\statelevelanalysis\\data\\covariates\\13_lscrocky.tif")
plot(template)

crs(template)

property <- st_transform(property, "+init=epsg:3577")

property <- property %>%
  mutate(totalwoodypix = exact_extract(woody, ., fun = "sum")) %>%
  mutate(totalwoodyarea = totalwoodypix * (100*100/10000)) %>% #inha
  mutate(propwoody = totalwoodyarea/AREA_H)

totalwoodyimg = fasterize(property, woody, field = "totalwoodyarea")
propwoodyimg = fasterize(property, woody, field = "propwoody")
propareaimg = fasterize(property, woody, field = "AREA_H")

writeRaster(totalwoodyimg, "E:\\01_NSW phase 2\\statelevelanalysis\\data\\covariates\\totalwoody.tif")
writeRaster(propwoodyimg, "E:\\01_NSW phase 2\\statelevelanalysis\\data\\covariates\\propwoody.tif")
writeRaster(propareaimg, "E:\\01_NSW phase 2\\statelevelanalysis\\data\\covariates\\propertyarea.tif")
