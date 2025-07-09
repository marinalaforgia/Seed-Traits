rm(list=ls())

# Extract and explore climate #
# 07/09/2025: This is kind of a mess at the moment...

library(raster)
library(sp)
library(tidyverse)
library(elevatr)
library(prism)
library(sf)
library(geodata)
library(SPEI)
library(infotheo)
library(dplyr)
library(tibble)
library(hydrostats)
library(tsibble)
library(entropy)

metadata <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")

#### bioclim vars ####
prism_set_dl_dir("Data/Climate/")
files <- prism_archive_ls()  # lists downloaded PRISM files

tmean_files <- prism_archive_subset("tmean", "annual", years = 2020)
r <- stack(tmean_files)
r <- worldclim_global("worldclim", var = "bio", res = 2.5)

names(r) <- c("MAT","Mean Diurnal Range", "Isothermality", "Temp Seasonality", "Max Temp", "Min Temp", "Temp Annual Range", "Mean Temp Wettest Quarter", "Mean Temp Driest Quarter", "Mean Temp Warmest Quarter", "Mean Temp Coldest Quarter", "Annual Precip", "Precip Wettest Month", "Precip Driest Month", "Precip Seasonality", "Precip Wettest Quarter", "Precip Driest Quarter" ,"Precip Warmest Quarter", "Precip Coldest Quarter")

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

seed.acc <- read.csv("Data/20230530_Seeds_All-Accessions.csv")
seed.acc <- filter(seed.acc, !is.na(lat))

coords <- data.frame(x = seed.acc$long, y = seed.acc$lat)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r, points)

df <- cbind.data.frame(seed.acc, coordinates(points), values)

#colnames(df)[4:22] <- c("MAT","Mean.Diurnal.Range", "Isothermality", "Temp.Seasonality", "Max.Temp", "Min.Temp", "Temp.Annual.Range", "Mean.Temp.Wettest.Quarter", "Mean.Temp.Driest.Quarter", "Mean.Temp.Warmest.Quarter", "Mean.Temp.Coldest.Quarter", "Annual.Precip", "Precip.Wettest.Month", "Precip.Driest.Month", "Precip.Seasonality", "Precip.Wettest.Quarter", "Precip.Driest.Quarter" ,"Precip.Warmest.Quarter", "Precip.Coldest.Quarter")

# Long term data set sites
sites <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")
site <- sites$Dataset
coords <- data.frame(x = sites$long, y = sites$lat)
points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r, points)
df2 <- cbind.data.frame(sites, coordinates(points), values)
df2 <- df2[,c(1,7,8,19:39)]
colnames(df2)[1] <- "site"

df <- merge(df, seed.acc.TW[,c(5,6)], by = "ID", all.x = T, all.y = F)

df <- df[,c(28, 2, 34:55)]
df <- df[,c(1,2,24,3:23)]
colnames(df)[3] <- "TW"
df <- rbind(df, df2)

# elevation data 
test <- get_elev_point(df[,4:5], units = "meters", prj = "EPSG:4326", src = "aws")

df <- cbind(df, test@data$elevation)
colnames(df)[28] <- "ele.m"

PCA <- prcomp(df[,6:24], scale = T)
summary(PCA)
biplot(PCA)

autoplot(PCA, x = 1, y = 2, data = df, frame = F, loadings = T, loadings.label = T, label = F, col = "AI") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = mean(df$AI))


df <- cbind(df, PCA$x[,1:3])

ggplot(df, aes(x = PC1, y = PC2, col = AI)) +
  geom_jitter(width = 0.1, height = 0.1) +
  #geom_point(aes(col = AI)) +
  #geom_point(data = df[df$site %in% site,], aes(x = PC1, y = PC2, col = AI)) +
  geom_text(data = df[df$site %in% site,], aes(label = site), col = "black") +
  theme_classic() +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = mean(df$AI))

# Thornwaite doesnt distinguish as much between sites as AI, mostly shows differences in elevation gradient
ggplot(df, aes(x = PC1, y = PC2, col = TW)) +
  geom_jitter(width = 0.1, height = 0.1) +
  #geom_point(aes(col = AI)) +
  #geom_point(data = df[df$site %in% site,], aes(x = PC1, y = PC2, col = AI)) +
  #theme_classic() +
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 80) +
    geom_text(data = df[df$site %in% site,], aes(label = site), col = "black") 
  #geom_text(aes(label = site), col = "black")



#### PRISM monthly ####
prism <- list()
sites <- c("McLaughlin", "Portal", "Carrizo-Plain", "Jasper-Ridge", "Sonoran-Desert")

for(i in sites) {
  prism[[i]] <- read.csv(paste0("Data/Long-Term-Datasets/", i ,"/Climate/PRISM_monthly.csv"))
  prism[[i]] <- prism[[i]][-c(1:10),]
  colnames(prism[[i]]) <- c("Date",	"ppt.mm",	"tmin.C",	"tmean.C", "tmax.C",	"vpdmin.hPa",	"vpdmaxhPa")
  prism[[i]]$Site <- i
  prism[[i]]$Date <- as.Date(paste0(prism[[i]]$Date,"-01"), "%Y-%m-%d")
  prism[[i]]$Year <- format(prism[[i]]$Date, "%Y")
  prism[[i]]$Month <- format(prism[[i]]$Date, "%m")
  prism[[i]] <- prism[[i]][,c(8,1,9,10,2:7)]
  prism[[i]][,5] <- as.numeric(prism[[i]][,5])
}

#turn prism into water-year annual prism data (Oct - Sept)



# Water Year rainfall
tmp <- data.frame()
bim <- c("Sonoran-Desert", "Portal")

for(i in sites) {
  month <- ddply(prism[[i]], .(Year, Month), summarize, precip = sum(ppt.mm))
  month$Month <- as.numeric(month$Month)
  month$Year <- as.numeric(month$Year)
  month <- filter(month, Year > 1982 & Year < 2020)
  sepdec.annual <- ddply(month, .(Year), summarize, precip = sum(precip[which(Month %in% 9:12)]))

    if(i %in% bim) {
    janaug.annual <- ddply(month, .(Year), summarize, precip = sum(precip[which(Month %in% 1:5)]))
  }
  
    else {
    janaug.annual <- ddply(month, .(Year), summarize, precip = sum(precip[which(Month %in% 1:8)]))
  }
  
  
  k <- data.frame(Site = i, Year = c(as.numeric(month$Year[1])+1):2019)
  k$ppt.yr <- sepdec.annual$precip[1:nrow(sepdec.annual)-1] + janaug.annual$precip[2:nrow(sepdec.annual)]
  tmp <- rbind(tmp, k)
}

site_colors <- c(
  "Sonoran-Desert" = "#ca562c",
  "Portal" = "#da8a5c",
  "Carrizo-Plain" = "#b4c8a8",
  "Jasper-Ridge" = "#80cdc1",
  "McLaughlin" = "#018571"
)

ggplot(tmp, aes(x = Year, y = ppt.yr, color = Site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = site_colors) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Year", y = "Annual precipitation (mm)")

cv <- ddply(tmp, .(Site), summarize, precip.cv = sd(ppt.yr)/mean(ppt.yr))

cv$Site <- recode_factor(cv$Site,
                         'Carrizo-Plain' = "Carrizo",
                         'Jasper-Ridge' = "Jasper Ridge",
                         'Sonoran-Desert' = "Sonoran")

meta <- merge(cv, meta, by.x = "Site", by.y = "Dataset")

# meta <- meta[,-1]
# meta <- meta[,c(2,1,3:18)]

write.csv(meta, "Data/Long-Term-Datasets/Datasets_metadata.csv", row.names = F)

# CV is much different depending on whether you include the summer monsoons or not, clear choice not to for SD where no summer census takes place, but portal is a different story


#### Colwell monthly precip ####

colwell <- data.frame(site = sites, Colwell.C = NA, Colwell.M = NA, Colwell.P = NA)
for(i in names(prism)){
  test <- prism[[i]][,c(2,5)]
  colnames(test) <- c("Date", "Q")
  test <- ts.format(test, format="%Y-%m-%d")
  colwell[colwell$site == i, ]$Colwell.P <- Colwells(test)$P # predictability
  colwell[colwell$site == i, ]$Colwell.M <- Colwells(test)$M # contingency
  colwell[colwell$site == i, ]$Colwell.C <- Colwells(test)$C # constancy
}

colwell$site <- recode_factor(colwell$site,
                         'Carrizo-Plain' = "Carrizo",
                         'Jasper-Ridge' = "Jasper Ridge",
                         'Sonoran-Desert' = "Sonoran")

meta <- merge(meta, colwell, by.y = "site", by.x = "Dataset")

#### Aridity Index ####
# Add AI
r <- raster("Spatial-Data/AI-Raster.tif")

##mask the image (two steps)
# tmpfilter <-  r < 5000
# tmpfilter <- r > 3000
# filtered_image <- mask(r, tmpfilter, maskvalue=1)
# 
# tmpfilter <- r > 4000
# filtered_image <- mask(r, tmpfilter, maskvalue=1)
# 
# 

# #copy your raster
# r.copy <- r
# 
# #set all pixels that are not contained in your vector to NA
# r.copy[r.copy > 3800] <- NA
# r.copy[r.copy < 2000] <- NA
# 
# leaflet() %>%
#   addTiles() %>%
#   addRasterImage(r.copy)


# spatial extent of the site varies so the buffer will vary
tmp <- st_as_sf(meta, coords = c("Longitude", "Latitude"), crs = "WGS84")
tmp <- st_transform(tmp, crs = 3310)

site_buffer <- st_buffer(tmp, 5000)

site_buffer <- st_transform(site_buffer, crs = "WGS84")

AI.site <- terra::extract(r, site_buffer, mean)

meta <- cbind(meta, AI.site)

# leaflet(site_buffer) %>%
#   addTiles() %>%
#   addPolygons

meta$AI.site <- meta$AI.site*0.0001

write.csv(meta, "Data/Long-Term-Datasets/Datasets_metadata.csv", row.names = F)

#### PRISM climate vars #####
sites <- metadata %>% select(Site = Dataset, lat = Latitude, lon = Longitude)

# 2. Set PRISM directory and load data
prism_set_dl_dir("Data/Climate/")

# Function to get Oct–May months across years
get_oct_may_indices <- function(years) {
  expand.grid(
    year = years,
    month = c(10:12, 1:5)
  ) %>% mutate(
    prism_year = ifelse(month >= 10, year, year + 1)
  )
}

# 3. Build raster stacks for Tmean and Precip
years <- 1981:2020
months <- sprintf("%02d", c(10:12, 1:5))

# Load PRISM tmean and ppt rasters
# Download monthly tmean for 1981–2020
# get_prism_monthlys(
#   type = "tmean",
#   years = 1981:2020,
#   mon = 1:12,
#   keepZip = FALSE
# )

tmean_files <- prism_archive_subset("tmean", "monthly", years = years)
ppt_files   <- prism_archive_subset("ppt", "monthly", years = years)

tmean_bils <- list.files(
  path = prism_get_dl_dir(), 
  pattern = "tmean.*\\.bil$", 
  full.names = TRUE, 
  recursive = TRUE
)

ppt_bils <- list.files(
  path = prism_get_dl_dir(), 
  pattern = "ppt.*\\.bil$", 
  full.names = TRUE, 
  recursive = TRUE
)

tmean_stack <- stack(tmean_bils)
ppt_stack   <- stack(ppt_bils)

# Works for both tmean_stack and ppt_stack
extract_dates <- function(stack) {
  yms <- str_extract(names(stack), "\\d{6}")
  as.Date(paste0(yms, "01"), format = "%Y%m%d")
}

t_dates <- extract_dates(tmean_stack)
p_dates <- extract_dates(ppt_stack)

# 4. Extract Oct–May values
extract_site_ts <- function(site_row) {
  lat <- site_row$lat
  lon <- site_row$lon
  site <- site_row$Site
  
  message("Processing site: ", site)
  
  temp_vals <- raster::extract(tmean_stack, matrix(c(lon, lat), ncol = 2))
  ppt_vals  <- raster::extract(ppt_stack, matrix(c(lon, lat), ncol = 2))
  
  if (all(is.na(temp_vals))) {
    message("  All temp values are NA")
    return(tibble())
  }
  if (all(is.na(ppt_vals))) {
    message("  All ppt values are NA")
    return(tibble())
  }
  
  # Date parsing
  t_dates <- extract_dates(tmean_stack)
  p_dates <- extract_dates(ppt_stack)
  
  df <- tibble(
    date = t_dates,
    tmean = as.numeric(temp_vals),
    ppt = as.numeric(ppt_vals),
    month = as.integer(format(t_dates, "%m")),
    year = as.integer(format(t_dates, "%Y"))
  ) %>%
    filter(month %in% c(10:12, 1:5)) %>%
    mutate(
      season_year = ifelse(month >= 10, year, year + 1)
    ) %>%
    group_by(season_year) %>%
    summarise(
      site = site,
      tmean_season = mean(tmean, na.rm = TRUE),
      ppt_season = sum(ppt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()
  
  message("  Extracted ", nrow(df), " rows for site")
  return(df)
}

# Apply to all sites
site_summaries <- sites %>% 
  split(.$Site) %>%
  map_dfr(extract_site_ts)

# 5. Calculate summary metrics per site
library(zoo)
library(infotheo)

  
site_colors <- c(
  "Sonoran" = "#ca562c",
  "Portal" = "#da8a5c",
  "Carrizo" = "#b4c8a8",
  "Jasper Ridge" = "#80cdc1",
  "McLaughlin" = "#018571"
)


ggplot(site_summaries, aes(y = ppt_season, x = season_year, col = site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = site_colors) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Year", y = "Annual precipitation (mm)")
  

#### Climate water deficit ####

# Build global CWD vector from all sites
all_cwd <- c()

all_cwd_df <- tibble()

for (i in 1:nrow(sites)) {
  lon <- sites$lon[i]
  lat <- as.numeric(sites$lat[i])
  
  tmean_vals <- raster::extract(tmean_stack, cbind(lon, lat))
  ppt_vals   <- raster::extract(ppt_stack,  cbind(lon, lat))
  
  if (is.null(tmean_vals) || is.null(ppt_vals) || any(is.na(tmean_vals)) || any(is.na(ppt_vals))) next
  if (length(tmean_vals) %% 12 != 0) next
  
  years <- length(tmean_vals) / 12
  tmean_mat <- matrix(tmean_vals, nrow = years, ncol = 12, byrow = TRUE)
  ppt_mat   <- matrix(ppt_vals,   nrow = years, ncol = 12, byrow = TRUE)
  
  pet_vec <- thornthwaite(as.vector(tmean_mat), lat)
  pet_mat <- matrix(pet_vec, nrow = years, ncol = 12, byrow = F)
  
  cwd_mat <- pmax(pet_mat - ppt_mat, 0)

  # Create year/month grid
  months_seq <- 1:12
  grid <- expand.grid(year = 1981:2020, month = months_seq)
  
  # Flatten CWD row-wise
  cwd_long <- as.vector(cwd_mat)  # Jan–Dec, Jan–Dec, ...
  
  # Add site, year, month, cwd to a data frame
  cwd_df <- tibble(
    site = sites$Site[i],
    year = grid$year,
    month = grid$month,
    cwd = cwd_long
  )
  
  # Store/append to master data frame
  all_cwd_df <- bind_rows(all_cwd_df, cwd_df)
}

all_cwd_df.yr <- all_cwd_df %>%
  group_by(year,site) %>%
  summarise(CWD_annual = sum(cwd, na.rm = TRUE))
  
ggplot(all_cwd_df.yr, aes(y = CWD_annual, x = year, group = site, col = site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = site_colors) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Year", y = "Climate Water Deficit")

seasonal_cwd <-all_cwd_df %>%
  filter(month %in% c(10, 11, 12, 1, 2, 3, 4, 5)) %>%  # Oct–May
  group_by(site, year) %>%
  summarise(CWD_OctMay = sum(cwd, na.rm = TRUE))

ggplot(seasonal_cwd, aes(y = CWD_OctMay, x = year, group = site, col = site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = site_colors) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Year", y = "Climate Water Deficit (Oct - May)")




#### Colwell indices on CWD ####

# M (Constancy) — How similar conditions are month to month and year to year
# 
# C (Contingency) — How strongly seasonal the variation is (e.g., rainfall always peaking in Nov)
# 
# P (Predictability) — Total regularity = M + C; how reliable cues are in predicting future environments

library(hydrostats)
library(dplyr)
library(tibble)

# List of unique sites
site_list <- unique(all_cwd_df$site)

# Store results
colwell_results <- list()

site_list <- unique(all_cwd_df$site)

library(lubridate)

colwell_results <- list()

for (s in site_list) {
  df <- all_cwd_df %>%
    filter(site == s) %>%
    arrange(year, month) %>%
    mutate(
      Date = as.POSIXct(paste(as.numeric(year), as.numeric(month), 15), format = "%Y %m %d"),
      Q = cwd
    ) %>%
    select(Date, Q)
  
  if (sum(!is.na(df$Q)) < 24) {
    colwell_results[[s]] <- tibble(site = s, Colwell.C = NA, Colwell.M = NA, Colwell.P = NA)
    next
  }
  
  df <- as.data.frame(df)
  
  col <- Colwells(df, s = 10)
  
  colwell_results[[s]] <- data.frame(
    site = s,
    Colwell.C.cwd = col[["C"]],
    Colwell.M.cwd = col[["M"]],
    Colwell.P.cwd = col[["P"]]
  )
}

colwell_df <- bind_rows(colwell_results)

#### Merge all together ####

test <- merge(meta[,-c(20:22)], colwell_df, by.x = "Dataset", by.y = "site")
test <- merge(test, final_summary, by.x = "Dataset", by.y = "site")
test <- merge(test, final_summary2, by.x = "Dataset", by.y = "site")

write.csv(test, "Data/Long-Term-Datasets/20250709_meta-new-climate.csv", row.names = F)
