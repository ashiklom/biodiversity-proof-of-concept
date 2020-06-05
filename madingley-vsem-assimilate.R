library(raster)

# Read biomass from VSEM
vsem <- brick("data/vsem-totb.tif")

# Read biomass from Madingley
mad_files <- list.files("data/madingley", "biomass.*\\.tif", full.names = TRUE)
mad_biomass <- stack(mad_files)
mad_biomass_sub <- crop(mad_biomass, vsem)

# Mask out the ocean (zero values)
mad_biomass_sub[mad_biomass_sub == 0] <- NA

vsem_reproj <- resample(vsem, mad_biomass_sub)

mad_range <- range(mad_biomass_sub)
vsem_range <- range(vsem_reproj)

par(mfrow = c(2, 2))
plot(mad_range[[1]], main = "Madingley min")
plot(mad_range[[2]], main = "Madingley max")
plot(vsem_range[[1]], main = "VSEM min")
plot(vsem_range[[2]], main = "VSEM max")

df <- mad_range - vsem_range

mad_vsem_overlap <- (mad_range[["range_max"]] > vsem_range[["range_min"]]) &
  (mad_range[["range_min"]] < vsem_range[["range_max"]])

dev.new()
plot(mad_vsem_overlap)
