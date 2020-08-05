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
  biomass_raw <- ncdf4::ncvar_get(nc, "autotrophbiomass density",
                                  start = c(tstart, 1, 1))
  # Madingley units are log ([g km-2] + 1)
  biomass_leaf <- udunits2::ud.convert(exp(biomass_raw) - 1, "g km-2", "kg m-2")
  # Foliar biomass is ~10% of total AGB
  # C.f. Bond-Lamberty et al. 2002
  biomass <- 10 * biomass_leaf
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

  bparam <- cbind(
    KEXT = rnorm_gt0(nens, 0.5, 0.15),
    LAR = rnorm_gt0(nens, 1.5, 0.3),
    # Default LUE: 0.002
    ## LUE = rlnorm(1, log(0.002), 0.4),
    LUE = rnorm_gt0(nens, 0.002, 0.0005),
    GAMMA = rnorm_gt0(nens, 0.4, 0.05),
    tauV = rnorm(nens, 1440, 10),
    tauS = rnorm(nens, 27370, 100),
    tauR = rnorm(nens, 1440, 10),
    Av = rnorm(nens, 0.5, 0.05),
    Cv = rnorm(nens, 3, 0.5),
    Cs = rnorm(nens, 15, 1),
    Cr = rnorm(nens, 3, 0.5)
  )

  bout <- array(numeric(), c(ny1, ny2, nens, 4))
  for (i in seq_len(ny1)) {
    message("Progress: ", i, " of ", ny1)
    for (j in seq_len(ny2)) {
      if (!all(is.na(merra_par[i,j,]))) {
        if (!is.null(vsem_Cv0)) {
          Cv0 <- vsem_Cv0[j,i]
        } else {
          Cv0 <- 3
        }
        if (!is.na(Cv0) & !Cv0 > 0) next
        par <- spline(x1, xout = x2, y = merra_par[i,j,])$y / 2
        for (s in seq_len(nens)) {
          vsem <- VSEM_random(par, Cv0)
          bparam[s,] <- vsem$params
          bout[i,j,s,] <- colMeans(vsem$vsem)
        }
      }
    }
  }
  list(bout = bout, bparam = bparam)
}

distribution_overlap <- function(r1, r2) {
  stopifnot(nrow(r1) == nrow(r2),
            ncol(r1) == ncol(r2))
  r1_mean <- r1[[1]]
  r1_sd <- r1[[2]]
  r2_mean <- r2[[1]]
  r2_sd <- r2[[2]]

  # This is currently numerical. But it can be done analytically, because these
  # are both normal.
  # Prec(x) = function(x) 1 / var(x)
  # Mean = (Prec1 * Mean1 + Prec2 * Mean2) / (Prec1 + Prec2)
  # Var = 1/(Prec1 + Prec2)

  f <- function(x, m1, s1, m2, s2) {
    f1 <- dnorm(x, m1, s1)
    f2 <- dnorm(x, m2, s2)
    pmin(f1, f2)
  }

  f_int <- function(VM, VS, MM, MS) {
    tryCatch(
      integrate(f, -Inf, Inf, VM, VS, MM, MS)$value,
      error = function(e) {
        warning(conditionMessage(e))
        return(NA)
      }
    )
  }

  intmat <- matrix(numeric(), nrow(r1), ncol(r1))
  for (i in seq_len(nrow(r1))) {
    for (j in seq_len(ncol(r1))) {
      if (is.na(r1_mean[i,j]) | is.na(r2_mean[i,j])) next
      intmat[i,j] <- f_int(r1_mean[i,j], r1_sd[i,j],
                           r2_mean[i,j], r2_sd[i,j])
    }
  }
  intrast <- raster(r1)
  intrast[] <- intmat
  intrast
}

joint_raster <- function(r1, r2) {
  stopifnot(nrow(r1) == nrow(r2), ncol(r1) == ncol(r2))
  m1 <- r1[[1]]
  s1 <- r1[[2]]
  m2 <- r2[[1]]
  s2 <- r2[[2]]
  p1 <- s1^-2
  p2 <- s2^-2
  pz <- p1 + p2
  mz <- (m1 * p1 + m2 * p2) / pz
  sz <- pz^-0.5
  raster::stack(list(Mean = mz, SD = sz))
}

raster_likelihood <- function(r, target, f = mean) {
  rsub <- raster::resample(r, target)
  pr <- rsub
  pr[] <- dnorm(rsub[], target[["Mean"]][], target[["SD"]][], log = TRUE)
  pr[pr < quantile(pr, 0.05)] <- NA
  f(pr[], na.rm = TRUE)
}

madingley_moose_raster <- function(f, base) {
  .datatable.aware <- TRUE #nolint
  # Area of a 1 degree grid cell
  # HACK: Assuming it's constant, but should really vary with latitude
  ft2km <- 30.48 / (100 * 1000)
  latsize <- 364000 * ft2km
  lonsize <- 288200 * ft2km
  mad_moose <- data.table::fread(f)
  mad_moose[, degree_area := latsize * lonsize]
  # BodyMass is in g -- here, restrict to mass > 300 kg
  # Individual biomass == Adult biomass
  dat <- mad_moose[IndividualBodyMass == AdultMass &
                     AdultMass > 400000,
                     .(density = sum(CohortAbundance) / degree_area),
                     .(Latitude, Longitude)]
  mad_moose_sf <- sf::st_as_sf(dat, coords = c("Longitude", "Latitude"),
                               crs = sf::st_crs(4326))
  rasterize(mad_moose_sf, base, field = "density")
}

land_shape <- function() {
  na <- rnaturalearth::ne_countries(
    continent = "North America",
    scale = "medium",
    returnclass = "sf"
  )
  lakes <- rnaturalearth::ne_load(
    scale = 50,
    type = "lakes",
    category = "physical",
    destdir = "data/naturalearth/ne_50m_lakes/",
    returnclass = "sf"
  ) %>%
    sf::st_crop(na)
  tm_shape(na) + tm_borders(lwd = 1) +
    tm_shape(lakes) + tm_borders(lwd = 1)
}
