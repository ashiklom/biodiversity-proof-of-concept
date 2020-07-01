library(BayesianTools)
library(ncdf4)
library(dplyr)
library(udunits2)
library(ggplot2)
library(sf)
## library(raster)

ncfiles <- list.files("data/merra", full.names = TRUE,
                      pattern = ".*\\.nc")

# Extract coordinates (only need to do this once)
## r <- brick(ncfiles[1], varname = "SWGNT")
## coord <- coordinates(r)

read_ncfile <- function(f) {
  # Determine time string
  dt_string <- regmatches(f, regexpr("[[:digit:]]{8}", f))
  dt <- strptime(dt_string, "%Y%m%d", tz = "UTC")
  dt_seq <- seq(dt, dt + 86370, "1 h")
  nc <- nc_open(f)
  on.exit(nc_close(nc), add = TRUE)
  dat <- expand.grid(
    lon = ncvar_get(nc, "lon"),
    lat = ncvar_get(nc, "lat"),
    datetime = dt_seq
  )
  dat$sw_net <- c(ncvar_get(nc, "SWGNT"))
  dat
}

dat <- do.call(rbind, lapply(ncfiles, read_ncfile))

dat2 <- as_tibble(dat) %>%
  mutate(sw_net_mj = ud.convert(sw_net, "W m-2", "MJ m-2 day-1"))

dat_daily <- dat2 %>%
  group_by(date = as.Date(datetime), lat, lon) %>%
  summarize(PAR = mean(sw_net_mj / 2)) %>%
  ungroup()

VSEM_random <- function(PAR) {
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
  VSEM(pars, PAR = PAR)
  # TODO: All I care about is final biomass; can save a lot of memory that way
}

out1 <- VSEM_random(dat_daily[["PAR"]]) %>%
  as_tibble() %>%
  bind_cols(dat_daily)

out2 <- st_as_sf(out1, coords = c("lon", "lat"), crs = 4326)

usamap <- rnaturalearth::ne_states(country = "united states of america",
                                   returnclass = "sf")

out2 %>%
  filter(date == date[[1]]) %>%
  ggplot() +
  geom_sf(data = st_crop(usamap, st_bbox(out2))) +
  geom_sf(aes(color = NEE), size = 10) +
  scale_color_viridis_c()

out1 %>%
  filter(lat == tail(lat, 1), lon == tail(lon, 1)) %>%
  ggplot() +
  aes(x = date, y = CR) +
  geom_line()

VSEM

## z <- out2 %>%
##   filter(geometry == tail(geometry, 1))

## vsem_out_raw <- VSEM(PAR = dat_daily[["PAR"]])
## vsem_out_ens <- bind_rows(lapply(
##   1:100,
##   function(i) mutate(dat_daily, iens = i, vsem = list(VSEM_random(PAR)))
## ))

vsem_out <- dat_daily %>%
  mutate(vsem = VSEM(PAR = PAR))

ggplot(dat_daily %>% filter(lat == lat[[1]], lon == lon[[1]])) +
  aes(x = date, y = PAR) +
  geom_line()

VSEM

dat$sw_net_mj <- udunits2::ud.convert(dat$sw_net, "W m-2", "MJ m-2 day-1")

library(ggplot2)
ggplot(subset(dat, lat == x$lat[[1]] & lon == x$lon[[1]])) +
  aes(x = datetime, y = sw_net_mj) +
  geom_line()

# Read

f <- ncfiles[1]
nc <- nc_open(f)

v <- getValues(r)

result <- array(numeric(0), c(nrow(v), ncol(v), 4))
for (i in seq_len(nrow(v))) {
  result[i,,] <- VSEM(PAR = v[i,])
}

plot(result[1,,3])

out <- apply(v, 1, function(x) VSEM(PAR = x))



lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")

dat <- ncvar_get(nc, "SWGNT")

nc

ncatt_get(nc, "SWGNT")

plot(dat[1,1,], type = "l")

PAR <- dat[1,1,]
plot(PAR, type = 'l')

pars
out <- BayesianTools::VSEM(pars, PAR = dat[1,1,])

VSEM_random(1)

VSEM(PAR=dat[1,1])

##################################################
loadd(mad_biomass)
loadd(vsem_biomass)

merra_coords <- merra_par$coords
mp <- merra_par$data

pt <- sp::SpatialPoints(t(c(merra_x[1], merra_y[1])), sp::CRS("+init=epsg:4326"))
rowMeans(raster::extract(mad_biomass_all, pt))

merra_xy <- dplyr::distinct(merra_coords, x, y)
merra_sp <- sf::st_as_sf(merra_xy, coords = c("x", "y"))
biomass <- raster::extract(mad_biomass_all, merra_sp)

library(raster)
library(stars)
merra_stars
a <- dplyr::slice(merra_stars, time, 1)
b <- raster::brick(a)
merra_stars[[1]]

##################################################
library(raster)
library(stars)
merra_r <- raster(merrafiles[1], varname = "LWGEM")
fpar_orig <- raster("data/modis_fpar.tif")
fpar_merra <- resample(fpar_orig, merra_r)

m <- t(as.matrix(fpar_merra))
m <- st_as_stars(fpar_merra)
dim(dplyr::pull(st_as_stars(fpar_merra)))

x <- as.matrix(fpar_merra)
dim(dplyr::pull(merra_stars))

##################################################
vsem_d <- density(vsem_reproj[10,50])
mad_d <- density(sample(mad_biomass[10,50], 5000, replace = TRUE), n = 20)
normalize <- function(x) x / max(x, na.rm = TRUE)
plot(vsem_d$x, normalize(vsem_d$y), type = "l")
curve(dnorm(x, mean(mad_biomass[10,50]), sd(mad_biomass[10,50])),
      add = TRUE, col = "green")

vsem_mean <- calc(vsem_reproj, mean)
vsem_sd <- calc(vsem_reproj, sd)
mad_mean <- mean(mad_biomass)
mad_sd <- calc(mad_biomass, sd)

combined <- brick(list(
  VM = vsem_mean, VS = vsem_sd,
  MM = mad_mean, MS = mad_sd
))

intmat <- matrix(numeric(), nrow(vsem_mean), ncol(vsem_mean))
for (i in 1:nrow(vsem_mean)) {
  for (j in 1:ncol(vsem_mean)) {
    if (is.na(mad_mean[i,j]) | is.na(vsem_mean[i,j])) next
    intmat[i,j] <- c_int(vsem_mean[i,j], vsem_sd[i,j],
                         mad_mean[i,j], mad_sd[i,j])
  }
}

intrast <- raster(combined)
intrast[] <- intmat

plot(intrast)

vsem_sd[1,1]

c_int2 <- Vectorize(c_int)

iint <- overlay(combined, fun = c_int2)

dim(vsem_reproj)

##################################################
# How to do normal-normal data assimlation...
i <- seq(-1, 3, 0.01)
mx <- 0
sx <- 0.3
my <- 2
sy <- 0.5
x <- dnorm(i, mx, sx)
y <- dnorm(i, my, sy)
px <- 1 / sx^2
py <- 1 / sy^2
mz <- (mx * px + my * py) / (px + py)
sz <- sqrt(1 / (px + py))
z <- dnorm(i, mz, sz)
normalize <- function(x) x / max(x)
plot(i, normalize(x), type = 'l')
lines(i, normalize(y), col = "red")
lines(i, normalize(x * y), col = "blue")
lines(i, normalize(z), col = "green4")

##################################################
library(raster)
library(tmap)
ft2km <- 30.48 / (100 * 1000)
latsize <- 364000 * ft2km
lonsize <- 288200 * ft2km
degree_area <- latsize * lonsize
mad_moose_1 <- data.table::fread("data/madingley/MooseCohorts_tstep1259_1.csv")

hist(mad_moose_1[["IndividualBodyMass"]])

dat <- mad_moose_1[IndividualBodyMass > 300000,
                   .(density = sum(CohortAbundance) / degree_area),
                   .(Latitude, Longitude)]
mad_moose_sf <- sf::st_as_sf(dat, coords = c("Longitude", "Latitude"),
                             crs = sf::st_crs(4326))
mad_moose_r <- rasterize(mad_moose_sf, moose_density, field = "density")
plot(mad_moose_r)

plot(moose_density)

nc <- ncdf4::nc_open("data/madingley/GridOutputs_SSPHistoric_1.nc")

data(World)
usca <- subset(World, iso_a3 %in% c("USA", "CAN"))

hist(mad_moose_r)

mm <- mad_moose_r
mm[mm > 50] <- 50
tm_shape(mm) +
  tm_raster(style = "cont",
            breaks = c(0, 1, 2, 3, 5, 10, 20, 30, 40, 50)) +
  tm_shape(usca) +
  tm_borders() +
  tm_legend(position = c("left", "bottom"))

hist(moose$Density)

plot(mad_moose_stats)
plot(moose_density)

l <- as.list(mad_moose_all)
