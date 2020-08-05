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
  ibest_moose = which.max(moose_overlap_probability),
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
    land <- land_shape()
    density_layer <- tm_raster(
      style = "cont",
      palette = "YlOrRd",
      title = expression("Moose density" ~ (individuals ~ km^{-2}))
      ## breaks = seq(0, 5, 1)
    )
    coords_best <- raster::coordinates(moose_stats)[ibest_moose, ]
    coords_best_shp <- sf::st_sfc(sf::st_point(coords_best), crs = 4326)
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
      tm_shape(coords_best_shp) +
      tm_dots(shape = 13, size = 2, border.lwd = 3) +
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
      theme_bw() +
      labs(x = "Biomass log(likelihood)",
           y = "Moose density log(likelihood)")
    ggsave(file_out("figures/madingley-likelihood.png"), plt,
           width = 5.62, height = 4.63, units = "in", dpi = 300)
  },
  joint_constraint = {
    dnorm0 <- function(x, ...) {
      x[x < 0] <- 0
      dnorm(x, ...)
    }
    ibest <- ibest_moose
    obs_moose_mean <- moose_stats[["Mean"]][ibest]
    obs_moose_sd <- moose_stats[["SD"]][ibest]
    mad_moose_mean <- mad_moose_stats[["Mean"]][ibest]
    mad_moose_sd <- mad_moose_stats[["SD"]][ibest]
    mad_moose_vals <- vapply(1:6, function(x) mad_moose_all[[x]][ibest], numeric(1))
    assim_moose_mean <- moose_joint_raster[["Mean"]][ibest]
    assim_moose_sd <- moose_joint_raster[["SD"]][ibest]
    obs_veg_mean <- mstmip_stats[[1]][ibest]
    obs_veg_sd <- mstmip_stats[[2]][ibest]
    mad_veg_mean <- mad_stats[["Mean"]][ibest]
    mad_veg_sd <- mad_stats[["SD"]][ibest]
    mad_veg_all <- vapply(1:6, function(x) mad_biomass[[x]][ibest], numeric(1))
    assim_veg_mean <- biomass_joint_raster[["Mean"]][ibest]
    assim_veg_sd <- biomass_joint_raster[["SD"]][ibest]

    png(file_out("figures/joint-distributions.png"),
        width = 14, height = 7,
        units = "in", res = 300)
    par(mfrow = c(1, 2))
    curve(dnorm0(x, obs_moose_mean, obs_moose_sd), 0, 3,
          n = 1000,
          xlab = expression("Moose density" ~ (km^-2)),
          ylab = "Probability",
          main = "Moose density")
    curve(dnorm0(x, mad_moose_mean, mad_moose_sd), 0, 3,
          n = 1000,
          add = TRUE, col = "red")
    curve(dnorm0(x, assim_moose_mean, assim_moose_sd), 0, 3,
          n = 1000,
          add = TRUE, col = "green")
    rug(mad_moose_vals, col = "red")
    text(mad_moose_vals, 0.1, seq_along(mad_moose_vals))
    legend("topright", c("Observed", "Madingley", "Joint"),
           lty = 1, col = c("black", "red", "green"))
    curve(dnorm0(x, obs_veg_mean, obs_veg_sd), 0, 15,
          n = 1000, col = "black",
          main = "Vegetation biomass",
          xlab = expression("Vegetation biomass" ~ (kg ~ m^-2)),
          ylab = "Probability",
          ylim = c(0, 0.4))
    curve(dnorm0(x, mad_veg_mean, mad_veg_sd), 0, 15,
          n = 1000, col = "red", add = TRUE)
    rug(mad_veg_all, col = "red")
    text(mad_veg_all, 0.01, seq_along(mad_veg_all))
    curve(dnorm0(x, assim_veg_mean, assim_veg_sd), 0, 15,
          n = 1000, col = "green", add = TRUE)
    legend("topright", c("MstMIP", "Madingley", "Joint"),
           lty = 1, col = c("black", "red", "green"))
    dev.off()

    xrange <- c(0, 3)
    yrange <- c(0, 15)
    rx <- seq(xrange[1], xrange[2], length.out = 100)
    ry <- seq(yrange[1], yrange[2], length.out = 100)
    px <- dnorm0(rx, assim_moose_mean, assim_moose_sd)
    py <- dnorm0(ry, assim_veg_mean, assim_veg_sd)
    pxy <- px %*% t(py)
    wts_u <- dnorm0(mad_moose_vals, assim_moose_mean, assim_moose_sd) *
      dnorm0(mad_veg_all, assim_veg_mean, assim_veg_sd)
    wts <- wts_u / sum(wts_u)
    mad_moose_wtd_mean <- weighted.mean(mad_moose_vals, wts)
    mad_moose_wtd_sd <- sqrt(sum(wts * (mad_moose_vals - mad_moose_wtd_mean)^2))
    mad_veg_wtd_mean <- weighted.mean(mad_veg_all, wts)
    mad_veg_wtd_sd <- sqrt(sum(wts * (mad_veg_all - mad_veg_wtd_mean)^2))
    pxw <- dnorm0(rx, mad_moose_wtd_mean, mad_moose_wtd_sd)
    pyw <- dnorm0(ry, mad_veg_wtd_mean, mad_veg_wtd_sd)

    png(file_out("figures/joint-constraint.png"),
        width = 7.5, height = 6.5,
        units = "in", res = 300)
    layout(matrix(c(2, 1, 0, 3), 2), c(4, 1), c(1, 4))
    par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
    plot(0, 0, type = "n",
         xlab = expression("Moose density" ~ (km^-2)),
         ylab = expression("Vegetation biomass" ~ (kg ~ m^-2)),
         xlim = xrange, ylim = yrange)
    image(rx, ry, pxy, add = TRUE)
    points(mad_moose_vals, mad_veg_all, pch = 19)
    text(mad_moose_vals, mad_veg_all, seq_along(mad_moose_vals),
         pos = 3)
    par(mar = c(1, 4, 0, 0), xpd = TRUE)
    plot(rx, px, type = "l", axes = FALSE, xlab = "", ylab = "",
         ylim = range(px, pxw), col = "green")
    lines(rx, pxw, col = "blue")
    legend("right", c("Univariate constraint", "Multivariate constraint"),
           lty = 1, col = c("green", "blue"), bty = "n")
    par(mar = c(4, 1, 0, 0), xpd = FALSE)
    plot(py, ry, type = "l", axes = FALSE, xlab = "", ylab = "",
         xlim = range(py, pyw), col = "green")
    lines(pyw, ry, col = "blue")
    dev.off()
  }
)
