library(raster)

# Read biomass from VSEM
vsem <- brick("data/vsem-agb.tif")

# Read biomass from Madingley
mad_files <- list.files("data/madingley", "biomass.*\\.tif", full.names = TRUE)
mad_biomass <- stack(mad_files)
mad_biomass_sub <- crop(mad_biomass, vsem)

par(mfrow = c(1, 2))
plot(mad_biomass_sub[[1]])
plot(vsem[[1]])
