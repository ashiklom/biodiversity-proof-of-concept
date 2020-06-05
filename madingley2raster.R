library(ncdf4)
library(raster)

ncfiles <- list.files("data/madingley", "*.nc$", full.names = TRUE)

rasterize <- function(x, lat, lon) {
  flip(raster(
    apply(x, c(3, 2), mean, na.rm = TRUE),
    xmn = min(lon), xmx = max(lon),
    ymn = min(lat), ymx = max(lat),
    crs = "+init=epsg:4326"
  ), 2)
}

for (f in ncfiles) {
  print(f)
  ifile <- regmatches(f, regexpr("([[:digit:]])(?=.nc)", f, perl = TRUE))
  nc <- nc_open(f)
  lat <- ncvar_get(nc, "Latitude")
  lon <- ncvar_get(nc, "Longitude")
  time <- ncvar_get(nc, "Time step")
  decid <- ncvar_get(nc, "deciduousbiomass density", start = c(tstart, 1, 1))
  ever <- ncvar_get(nc, "evergreenbiomass density", start = c(tstart, 1, 1))
  nc_close(nc)
  tstart <- which(time > 800)[1]
  decid_r <- rasterize(decid, lat, lon)
  ever_r <- rasterize(ever, lat, lon)
  tot_biomass <- decid_r + ever_r
  writeRaster(tot_biomass, sprintf("data/madingley/biomass_%s.tif",
                                   ifile))
}
