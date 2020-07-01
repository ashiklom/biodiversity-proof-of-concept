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
  mad_biomass_all = raster::stack(file_in(!!outfiles)),
  mad_biomass = {
    mad_biomass_sub <- raster::crop(mad_biomass_all, merra_fpar)
    mad_biomass_sub[mad_biomass_sub <= 0] <- NA
    mad_biomass_sub
  },
  mstmip_stats = raster::resample(
    raster::brick(file_in("data/mstmip-stats.tif")),
    mad_biomass
  ),
  mad_stats = raster::stack(list(
    Mean = mean(mad_biomass),
    SD = calc(mad_biomass, sd)
  )),
  distribution_raster = distribution_overlap(mstmip_stats, mad_stats),
  distribution_plot = {
    png(file_out("figures/biomass-overlap-distribution.png"),
        width = 7, height = 7, units = "in", res = 300)
    plot(distribution_raster, main = "Fraction distribution overlap")
    dev.off()
  },
  mad_veg_compare = {
    png(file_out("figures/mad-veg-compare.png"),
        width = 10, height = 7, units = "in", res = 300)
    s <- raster::stack(list(
      Madingley = mad_stats[["Mean"]],
      MstMIP = mstmip_stats[[1]]
    ))
    plot(s)
    dev.off()
  },
  biomass_joint_raster = joint_raster(mstmip_stats, mad_stats),
  mad_likelihoods = tibble::tibble(
    file = file_in(!!outfiles),
    biomass_likelihood = purrr::map_dbl(
      file,
      biomass_likelihood,
      biomass = biomass_joint_raster
    )
  )
)
