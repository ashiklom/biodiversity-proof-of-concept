library(drake)
library(raster)

ncfiles <- list.files("data/madingley", "*.nc$", full.names = TRUE)
nc_id <- regmatches(ncfiles, regexpr("([[:digit:]])(?=.nc)", ncfiles,
                                        perl = TRUE))
outfiles <- sprintf("data/madingley/biomass_%s.tif", nc_id)
merrafiles <- list.files("data/merra-monthly", "*.nc", full.names = TRUE)

dir.create("figures", showWarnings = FALSE)

plan <- drake_plan(
  madingley_biomass = target(
    madingley_get_biomass(file_in(ncfile), file_out(outfile)),
    transform = map(ncfile = !!ncfiles, outfile = !!outfiles)
  ),
  merra_stars = prepare_merra_par(merrafiles),
  merra_fpar = {
    ms <- dplyr::slice(merra_stars, time, 1)
    ff <- tempfile(fileext = ".tif")
    stars::write_stars(ms, ff)
    merra_r <- raster(ff)
    fpar_orig <- raster(file_in("data/modis_fpar.tif"))
    resample(fpar_orig, merra_r)
  },
  vsem_raw = vsem_grid_ensemble(dplyr::pull(merra_stars),
                                3 * as.matrix(merra_fpar) / 100),
  vsem_biomass = {
    bb <- sf::st_bbox(merra_stars)
    r <- raster::brick(
      vsem_raw$bout[,,,2],
      xmn = bb["xmin"], xmx = bb["xmax"],
      ymn = bb["ymin"], ymx = bb["ymax"],
      transpose = TRUE,
      crs = "+init=epsg:4326")
    r[r <= 0] <- 0
    r
  },
  mad_biomass_all = raster::stack(file_in(!!outfiles)),
  mad_biomass = {
    mad_biomass_sub <- raster::crop(mad_biomass_all, vsem_biomass)
    mad_biomass_sub[mad_biomass_sub <= 0] <- NA
    mad_biomass_sub
  },
  vsem_reproj = raster::resample(vsem_biomass, mad_biomass),
  mad_range = range(mad_biomass),
  vsem_range = range(vsem_reproj),
  rangeplot = {
    png(file_out("figures/biomass-range-compare.png"),
        width = 7, height = 7, units = "in", res = 300)
    par(mfrow = c(2, 2))
    plot(mad_range[[1]], main = "Madingley min")
    plot(mad_range[[2]], main = "Madingley max")
    plot(vsem_range[[1]], main = "VSEM min")
    plot(vsem_range[[2]], main = "VSEM max")
    dev.off()
  },
  overlapplot = {
    mad_vsem_overlap <- (mad_range[["range_max"]] > vsem_range[["range_min"]]) &
      (mad_range[["range_min"]] < vsem_range[["range_max"]])
    png(file_out("figures/biomass-range-intersect.png"),
        width = 5, height = 5, units = "in", res = 300)
    plot(mad_vsem_overlap)
    dev.off()
  }
)
