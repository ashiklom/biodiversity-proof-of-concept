makeraster <- function(x, lat, lon) {
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
  biomass_r <- makeraster(biomass, lat, lon)
  raster::writeRaster(biomass_r, outfile, overwrite = TRUE)
  invisible(outfile)
}

rnorm_gt0 <- function(...) {
  x <- -1
  while (any(x < 0)) {
    x <- rnorm(...)
  }
  x
}

VSEM_random <- function(PAR, Cv0 = 3) {
  # TODO: Need to save the parameter values
  pars <- c(
    KEXT = rnorm_gt0(1, 0.5, 0.15),
    LAR = rnorm_gt0(1, 1.5, 0.3),
    LUE = rnorm_gt0(1, 0.002, 0.001),
    GAMMA = rnorm(1, 0.4, 0.1),
    tauV = rnorm(1, 1440, 10),
    tauS = rnorm(1, 27370, 100),
    tauR = rnorm(1, 1440, 10),
    Av = rnorm(1, 0.5, 0.05),
    Cv = rnorm(1, Cv0, 0.5),
    Cs = rnorm(1, Cv0*5, 1),
    Cr = rnorm(1, Cv0, 0.5)
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
  sf::st_crop(s, land_sub)
}

vsem_grid_ensemble <- function(merra_par, vsem_Cv0 = NULL, nens = 100) {
  # HACK: X vectors for interpolating from monthly to daily
  x1 <- seq(1, 12)
  x2 <- seq(1, 12, length.out = 365)
  ny1 <- dim(merra_par)[1]
  ny2 <- dim(merra_par)[2]

  vsem_test <- VSEM_random(1)
  bout <- array(numeric(), c(ny1, ny2, nens, 4))
  bparam <- matrix(numeric(), nens, length(vsem_test$params))
  colnames(bparam) <- names(vsem_test$params)
  for (i in seq_len(ny1)) {
    for (j in seq_len(ny2)) {
      if (!all(is.na(merra_par[i,j,]))) {
        if (!is.null(vsem_Cv0)) {
          Cv0 <- vsem_Cv0[j,i]
        } else {
          Cv0 <- 3
        }
        if (!Cv0 > 0) next
        par <- spline(x1, xout = x2, y = merra_par[i,j,])$y / 2
        for (s in seq_len(nens)) {
          vsem <- VSEM_random(par, Cv0)
          bparam[s,] <- vsem$params
          bout[i,j,s,] <- colMeans(vsem$vsem)
        }
      }
    }
  }
  list(bout = bout, bpaam = bparam)
}
