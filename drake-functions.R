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
  ## rsub <- raster::resample(r, target)
  rsub <- r
  pr <- rsub
  pr[] <- dnorm(rsub[], target[["Mean"]][], target[["SD"]][], log = TRUE)
  pr[pr < quantile(pr, 0.05)] <- NA
  f(pr[], na.rm = TRUE)
}

madingley_moose_raster <- function(f, base) {
  .datatable.aware <- TRUE #nolint
  mad_moose <- data.table::fread(f)
  # BodyMass is in g -- here, restrict to mass > 400 kg
  # Individual biomass == Adult biomass
  dat <- mad_moose[IndividualBodyMass == AdultMass &
                     AdultMass > 400000,
                     .(abundance = sum(CohortAbundance)),
                     .(Latitude, Longitude)]
  mad_moose_sf <- sf::st_as_sf(dat, coords = c("Longitude", "Latitude"),
                               crs = sf::st_crs(4326))
  rasterize(mad_moose_sf, base, field = "abundance")
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
