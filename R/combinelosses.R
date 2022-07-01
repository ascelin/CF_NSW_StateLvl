library(terra)
library(purrr)

path <- ("../data/loss")

loss.data <- dir(path, full.names=T, pattern = "majorityrule")

loss.data

lossrasters <- loss.data %>%
                 map(., rast)

mosaiced <- do.call(mosaic, lossrasters)

writeRaster(mosaiced, "../data/loss/nsw_state_allagents_resampled100m_majorityrule.tif")
