# Run VSEM simulations
library(sf)
library(stars)
library(rnaturalearth)
library(tidyverse)
library(raster)

# Country boundaries polygon
land <- ne_countries(continent = "North America", returnclass = "sf")

# Madingley maximum is 64.5
lat_max <- 64.5
lat_min <- 35
lon_min <- -165
lon_max <- -50
bbox <- st_bbox(c(xmin = lon_min, xmax = lon_max,
                  ymin = lat_min, ymax = lat_max),
                crs = 4326)

# Subset to the Madingley polygon
land_sub <- st_crop(land, bbox)

# Grab the monthly MERRA output, as a single "stars" object
merrafiles <- list.files("data/merra-monthly", "*.nc", full.names = TRUE)
s <- read_stars(merrafiles, sub = "SWGDN")
st_crs(s) <- st_crs(4326)

# Subset to just the North America land
s_sub <- st_crop(s, land_sub)

# Extract the SW data
y <- pull(s_sub)

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
  BayesianTools::VSEM(pars, PAR = PAR)
}

# HACK: X vectors for interpolating from monthly to daily
x1 <- seq(1, 12)
x2 <- seq(1, 12, length.out = 365)

# Aboveground, Belowground, Total (?)
outfile <- "data/vsem-output-raw.rds"
if (!file.exists(outfile)) {
  bout <- array(numeric(), c(dim(y)[1:2], 100, 3))
  pb <- progress::progress_bar$new(total = prod(dim(y)[1:2]))
  for (i in seq_len(dim(y)[1])) {
    for (j in seq_len(dim(y)[2])) {
      if (!all(is.na(y[i,j,]))) {
        # HACK: Interpolate monthly data to daily
        par <- spline(x1, xout = x2, y = y[i,j,])$y / 2
        for (s in 1:100) {
          vsem <- VSEM_random(par)
          bout[i,j,s,] <- colMeans(vsem[,-1])
        }
      }
      pb$tick()
    }
  }
  saveRDS(bout, outfile)
} else {
  bout <- readRDS(outfile)
}

s_bbox <- st_bbox(s_sub)
bout_agb <- brick(bout[,,,1],
                  xmn = s_bbox["xmin"], xmx = s_bbox["xmax"],
                  ymn = s_bbox["ymin"], ymx = s_bbox["ymax"],
                  transpose = TRUE,
                  crs = "+init=epsg:4326")
writeRaster(bout_agb, "data/vsem-agb.tif")

bout_totb <- brick(bout[,,,1] + bout[,,,2],
                  xmn = s_bbox["xmin"], xmx = s_bbox["xmax"],
                  ymn = s_bbox["ymin"], ymx = s_bbox["ymax"],
                  transpose = TRUE,
                  crs = "+init=epsg:4326")
writeRaster(bout_totb, "data/vsem-totb.tif")
