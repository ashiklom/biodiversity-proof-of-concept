library(raster)

# Read biomass from VSEM
vsem <- brick("data/vsem-totb.tif")

# Read biomass from Madingley
mad_files <- list.files("data/madingley", "biomass.*\\.tif", full.names = TRUE)
mad_biomass <- stack(mad_files)
mad_biomass_sub <- crop(mad_biomass, vsem)

vsem_reproj <- resample(vsem, mad_biomass_sub)
mad_biomass_sub

mad_range <- range(mad_biomass_range)
vsem_range <- range(vsem_reproj)

df <- mad_range - vsem_range

mad_vsem_overlap <- (mad_range[["range_max"]] > vsem_range[["range_min"]]) &
  (mad_range[["range_min"]] < vsem_range[["range_max"]])

plot(mad_vsem_overlap)
