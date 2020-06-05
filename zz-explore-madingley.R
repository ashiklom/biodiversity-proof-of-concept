# Create the region of interest
library(sf)
library()

library(tidyverse)
library(sf)

## library(stars)

## ncmeta::nc_meta("data/madingley/GridOutputs_SSPHistoric_1.nc")
## nc1 <- read_ncdf("data/madingley/GridOutputs_SSPHistoric_1.nc",
##                  curvilinear = c("Longitude", "Latitude"))
## nc1 <- read_stars("data/madingley/GridOutputs_SSPHistoric_1.nc",
##                   proxy = TRUE)

## b <- brick("data/madingley/GridOutputs_SSPHistoric_1.nc")
## r1 <- brick("data/madingley/GridOutputs_SSPHistoric_1.nc",
##              var = "deciduousbiomass density")

## plot(r1[[1]])



dv <- apply(decid, c(2, 3), var)

ii <- which(dv != 0, arr.ind = TRUE)[1:5,]

plot(time, decid[, ii[1,1], ii[1,2]], type = 'l')

dsub <- decid[, ii[,1], ii[,2]]
matplot(, type = 'l')
## plot(decid[,,85], type = 'l')


dat <- read_csv("data/madingley/MooseCohorts_tstep1259_1.csv")

dat %>%
  count(function)

dat_sf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  rename(gpoints = geometry) %>%
  mutate(geometry = st_make_grid(gpoints))

na <- ne_countries(continent = "North America")

ggplot(dat_sf) +
  aes(fill = CohortAbundance) +
  geom_sf(data = na) +
  geom_sf() +
  scale_fill_viridis_c()
