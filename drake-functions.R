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

VSEM_random <- function(PAR) {
  # TODO: Need to save the parameter values
  pars <- c(
    KEXT = rnorm(1, 0.5, 0.05),
    LAR = rnorm(1, 1.5, 0.1),
    LUE = rnorm(1, 0.002, 0.0001),
    GAMMA = rnorm(1, 0.4, 0.05),
    tauV = rnorm(1, 1440, 10),
    tauS = rnorm(1, 27370, 100),
    tauR = rnorm(1, 1440, 10),
    Av = rnorm(1, 0.5, 0.05),
    Cv = rnorm(1, 3, 0.5),
    Cs = rnorm(1, 15, 1),
    Cr = rnorm(1, 3, 0.5)
  )
  list(
    vsem = BayesianTools::VSEM(pars, PAR = PAR),
    params = pars
  )
}

prepare_merra_par <- function(merrafiles) {
  # Country boundaries polygon
  land <- rnaturalearth::ne_countries(
    continent = "North America",
    returnclass = "sf"
  )

  # Madingley maximum is 64.5
  lat_max <- 64.5
  lat_min <- 35
  lon_min <- -165
  lon_max <- -50
  bbox <- sf::st_bbox(c(xmin = lon_min, xmax = lon_max,
                        ymin = lat_min, ymax = lat_max),
                      crs = 4326)
  # Subset the land mask to the Madingley polygon
  land_sub <- sf::st_crop(land, bbox)

  # Read MERRA
  s <- stars::read_stars(merrafiles, sub = "SWGDN")
  sf::st_crs(s) <- sf::st_crs(4326)
  # Subset to just the North America land
  s_sub <- sf::st_crop(s, land_sub)
  # Extract the SW data
  list(
    data = dplyr::pull(s_sub),
    bbox = sf::st_bbox(s_sub)
  )
}

vsem_grid_ensemble <- function(merra_par, nens = 100) {
  # HACK: X vectors for interpolating from monthly to daily
  x1 <- seq(1, 12)
  x2 <- seq(1, 12, length.out = 365)
  ny1 <- dim(merra_par)[1]
  ny2 <- dim(merra_par)[2]

  vsem_test <- VSEM_random(1)
  bout <- array(numeric(), c(ny1, ny2, nens, 4))
  bparam <- matrix(numeric(), nens, length(vsem_test$params))
  colnames(bparam) <- names(vsem_test$params)
  pb <- progress::progress_bar$new(total = ny1 * ny2)
  for (i in seq_len(ny1)) {
    for (j in seq_len(ny2)) {
      if (!all(is.na(merra_par[i,j,]))) {
        par <- spline(x1, xout = x2, y = merra_par[i,j,])$y / 2
        for (s in seq_len(nens)) {
          vsem <- VSEM_random(par)
          bparam[s,] <- vsem$params
          bout[i,j,s,] <- colMeans(vsem$vsem)
        }
      }
      pb$tick()
    }
  }
  list(bout = bout, bpaam = bparam)
}
