library(drake)
library(raster)
library(tmap)
library(ggplot2)
library(magrittr, include.only = "%>%")

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
      purrr::map(file, raster::raster),
      raster_likelihood,
      target = biomass_joint_raster
    ),
    moose_likelihood = purrr::map_dbl(
      as.list(mad_moose_all),
      raster_likelihood,
      target = moose_joint_raster
    )
  ),
  moose_density = {
    # Units are individuals per sq km
    moose <- sf::read_sf(file_in("data/MooseDensity2010/MooseDensity2010.shp"))
    rasterize(moose, biomass_joint_raster, field = "Density", fun = mean)
  },
  moose_stats = stack(list(
    Mean = moose_density,
    # HACK: Assume SD = 10% of mean (so CI is +/- 20%)
    SD = moose_density * 0.1
  )),
  mad_moose_all = {
    madfiles <- file_in(!!fs::dir_ls("data/madingley/", glob = "*Moose*.csv"))
    stack(lapply(madfiles, madingley_moose_raster, base = mad_stats))
  },
  mad_moose_stats = stack(list(
    Mean = mean(mad_moose_all),
    SD = calc(mad_moose_all, sd)
  )),
  moose_overlap_probability = distribution_overlap(
    moose_stats, mad_moose_stats
  ),
  moose_joint_raster = joint_raster(moose_stats, mad_moose_stats),
  biomass_plot = {
    land <- land_shape()
    brks <- seq(0, 20, 5)
    biomass_layer <- tm_raster(
      style = "cont",
      palette = "Greens",
      title = expression("Mean biomass" ~ (kg ~ m^{-2})),
      breaks = brks
    )
    p1 <- tm_shape(mstmip_stats[[1]]) +
      biomass_layer +
      land +
      tm_layout(
        title = "MstMIP",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    p2 <- tm_shape(mad_stats[["Mean"]]) +
      biomass_layer +
      land +
      tm_layout(
        title = "Madingley",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    p3 <- tm_shape(distribution_raster) +
      tm_raster(
        style = "cont",
        palette = "BuPu",
        title = "Probability (0-1)"
      ) +
      land +
      tm_layout(
        title = "Biomass overlap",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    tmap_save(
      tmap_arrange(p1, p2, p3, ncol = 1),
      file_out("figures/biomass-compare.png"),
      width = 7.3, height = 9, units = "in", dpi = 300
    )
  },
  moose_plot = {
    land <- rnaturalearth::ne_countries(
      continent = "North America",
      scale = "medium"
    )
    land <- land_shape()
    density_layer <- tm_raster(
      style = "cont",
      palette = "YlOrRd",
      title = expression("Moose density" ~ (individuals ~ km^{-2}))
      ## breaks = brks
    )
    p1 <- tm_shape(moose_stats[["Mean"]]) +
      density_layer +
      land +
      tm_layout(
        title = "Observed",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    p2 <- tm_shape(mad_moose_stats[["Mean"]]) +
      density_layer +
      land +
      tm_layout(
        title = "Madingley",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    p3 <- tm_shape(moose_overlap_probability) +
      tm_raster(
        style = "cont",
        palette = "BuPu",
        title = "Probability (0-1)"
      ) +
      land +
      tm_layout(
        title = "Moose density overlap",
        title.position = c("left", "bottom"),
        legend.position = c("left", "bottom")
      )
    tmap_save(
      tmap_arrange(p1, p2, p3, ncol = 1),
      file_out("figures/moose-compare.png"),
      width = 7.3, height = 9, units = "in", dpi = 300
    )
  },
  mad_likelihood_plot = {
    dat <- mad_likelihoods
    dat$ensemble <- regmatches(
      basename(dat$file),
      regexpr("[[:digit:]]+", basename(dat$file))
    )
    plt <- ggplot(dat) +
      aes(x = biomass_likelihood, y = moose_likelihood,
          label = ensemble) +
      geom_label() +
      ## geom_segment(aes(x = -7350, xend = -7200,
      ##                  y = -6.5e9, yend = -5.5e9),
      ##              arrow = arrow(ends = "both")) +
      ## geom_text(aes(x = -7200, y = -5.5e9, label = "More likely"),
      ##           nudge_x = -50) +
      ## geom_text(aes(x = -7350, y = -6.5e9, label = "Less likely"),
      ##           nudge_x = 50) +
      theme_bw() +
      labs(x = "Biomass log(likelihood)",
           y = "Moose density log(likelihood)")
    ggsave(file_out("figures/madingley-likelihood.png"), plt,
           width = 5.62, height = 4.63, units = "in", dpi = 300)
  }
)
