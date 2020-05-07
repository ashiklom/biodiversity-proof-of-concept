library(BayesianTools)
library(ncdf4)
library(dplyr)
library(udunits2)
## library(raster)

ncfiles <- list.files("data/merra", full.names = TRUE,
                      pattern = ".*\\.nc")

# Extract coordinates (only need to do this once)
r <- brick(ncfiles[1], varname = "SWGNT")
coord <- coordinates(r)

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
  summarize(PAR = mean(sw_net_mj)) %>%
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
}

vsem_out_l <- list()
for (i in 1:10) {
  vsem_out_l[[i]] <- VSEM_random(dat_daily[["PAR"]])
}

## vsem_out_raw <- VSEM(PAR = dat_daily[["PAR"]])
vsem_out_ens <- bind_rows(lapply(
  1:100,
  function(i) mutate(dat_daily, iens = i, vsem = list(VSEM_random(PAR)))
))

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

VSEM(PAR=dat[1,1])
