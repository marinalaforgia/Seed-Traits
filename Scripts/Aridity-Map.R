rm(list=ls())

library(ggplot2)
library(raster)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)

# Load background raster (e.g. aridity, temp, etc.)
r <- raster("Spatial-Data/AI-Raster.tif")  # Replace with your raster file

# Convert raster to dataframe for ggplot
r_df <- as.data.frame(r, xy = TRUE)
names(r_df)[3] <- "value"

# Load state boundaries
# Get US state boundaries as an sf object
states <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf")

# Subset to Southwest states
sw_states <- states[states$name %in% c("California", "Arizona"), ]

# Load yellow points
pts <- read.csv("Data/20230530_Seeds_All-Accessions.csv")  # must have lat, lon
traits <- read.csv("Data/20230801_Seed-Traits_clean_ID.csv")
traits <- filter(traits, keep != "discard")
pts <- filter(pts, ID %in% traits$ID)
pts <- pts[!is.na(pts$lat),]
pts_sf <- st_as_sf(pts, coords = c("long", "lat"), crs = 4326)

# Get CA and AZ state shapes
states <- ne_states(country = "United States of America", returnclass = "sf")
az_ca <- states %>% filter(name %in% c("Arizona", "California"))
az_ca <- st_transform(az_ca, crs = crs(r))

# Clip raster
r_crop <- crop(r, extent(az_ca))
r_masked <- mask(r_crop, az_ca)
r_df <- as.data.frame(r_masked, xy = TRUE)
names(r_df)[3] <- "value"

# Transform points to match raster CRS
pts_sf <- st_transform(pts_sf, crs = crs(r))

# Ensure CRS match
az_ca <- st_transform(az_ca, crs = crs(r))

# Crop and mask
r_cropped <- crop(r, extent(az_ca))
r_masked <- mask(r_cropped, az_ca)

site_colors <- c(
  "Sonoran" = "#ca562c",
  "Portal" = "#da8a5c",
  "Carrizo" = "#b4c8a8",
  "Jasper Ridge" = "#80cdc1",
  "McLaughlin" = "#018571"
)

# Rescale raster from original range (100–23000) to 0–1
r_masked <- setMinMax(r_masked)
r_scaled <- (r_masked - minValue(r_masked)) / (maxValue(r_masked) - minValue(r_masked))
r_scaled[r_scaled < 0] <- 0
r_scaled[r_scaled > 1] <- 1


# Final plot
# New labels and breakpoints

# < 0.03 Hyper Arid 
# 0.03 – 0.15 Arid
# 0.15 - 0.2
# 0.2 - 0.4
# 0.4 – 0.5 Semi-Arid
# 0.5 – 0.65 Dry sub-humid
# > 0.65 Humid

# 0.12, 0.14, 0.19, 0.40, 0.50, 

r_binned <- cut(r_masked[],
                breaks = c(-Inf, 300, 1500, 2000, 4000, 5000, 6500, Inf),
                labels = c("<0.03", "0.03-0.15", "0.15-0.2", "0.2–0.4", "0.4–0.5", "0.5–0.65", ">0.65")
)

# Assign binned values back into a raster
r_class <- r_masked
values(r_class) <- r_binned

r_df <- as.data.frame(r_class, xy = TRUE)
names(r_df)[3] <- "aridity_bin"
r_df <- filter(r_df, !is.na(aridity_bin))

r_df$aridity_bin <- factor(r_df$aridity_bin, levels = c("<0.03", "0.03-0.15", "0.15-0.2", "0.2–0.4", "0.4–0.5", "0.5–0.65", ">0.65"))

# Add in sites
meta <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")

# Convert meta to sf
sites_sf <- st_as_sf(meta, coords = c("Longitude", "Latitude"), crs = 4326)
sites_sf <- st_transform(sites_sf, crs = st_crs(az_ca))  # match CRS

# Plot with stars for site locations
map <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = aridity_bin)) +
  scale_fill_manual(
    name = "Aridity Index",
    values = c(
      "<0.03" = "#7f0000",     
      "0.03-0.15" = "#ca562c",
      "0.15-0.2" = "#da8a5c",
      "0.2–0.4" = "#b4c8a8",
      "0.4–0.5" = "#80cdc1",
      "0.5–0.65" = "#018571",  
      ">0.65" = "#004d40"      
    )
  ) +
  geom_sf(data = az_ca, fill = NA, color = "black", linewidth = 0.5) +
  geom_sf(data = pts_sf, shape = 21, fill = "gold", color = "black", size = 2, stroke = 0.4) +
  #geom_sf(data = sites_sf, shape = 24, color = "black", size = 3, fill = "magenta3", stroke = 0.7, alpha = 0.8) +  # stars
  coord_sf(xlim = c(-125, -107), ylim = c(30, 43), expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) #+
  #geom_sf_text(data = sites_sf, aes(label = Dataset), size = 5, nudge_y = 1)

ggsave("Manuscript/Aridity/map-no-sites.png", map, width = 7, height = 5, units = "in")


