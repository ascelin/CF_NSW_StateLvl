library(terra)

woody <- rast("../data/woody/nsw_state_woody_100m_majorityrule.tif")

privateproperty <- rast("../data/propertydata/property_privateland.tif")


woody.on.private <- woody %>%
                      mask(., privateproperty)

writeRaster(woody.on.private, "../data/woodyonprivateland/woodyonprivate.tif")
