rasterize <- function(x, lat, lon) {
  raster::flip(raster::raster(
    apply(x, c(3, 2), mean, na.rm = TRUE),
    xmn = min(lon), xmx = max(lon),
    ymn = min(lat), ymx = max(lat),
    crs = "+init=epsg:4326"
  ), 2)
}

madingley_get_biomass <- function(ncfile, outfile) {
  nc <- ncdf4::nc_open(ncfile)
  lat <- ncdf4::ncvar_get(nc, "Latitude")
  lon <- ncdf4::ncvar_get(nc, "Longitude")
  time <- ncdf4::ncvar_get(nc, "Time step")
  tstart <- which(time > 800)[1]
  biomass <- ncdf4::ncvar_get(nc, "autotrophbiomass density",
                              start = c(tstart, 1, 1))
  ncdf4::nc_close(nc)
  biomass_r <- rasterize(biomass, lat, lon)
  raster::writeRaster(biomass_r, outfile, overwrite = TRUE)
  invisible(outfile)
}
