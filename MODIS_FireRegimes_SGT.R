library(MODISTools)
library(sf)
library(raster)
library(terra)
library(sp)
library(tidyverse)


# Read and plot points coordinates
traps.pts <- read.table("Data/Points_Traps_Bruna_Heitor.txt", h = T)
traps.pts

traps.pts.sp <- traps.pts

coordinates(traps.pts.sp) <- c("Longitude", "Latitude")
proj4string(traps.pts.sp) <- CRS("+proj=longlat +datum=WGS84")
traps.pts.sp
traps.pts.sp <- st_as_sf(traps.pts.sp)

mt_products()
# MODIS/Terra+Aqua Burned Area (Burned Area) Monthly L3 Global 500 m SIN Grid
mt_bands(product = "MCD64A1")

# MODIS/Terra+Aqua Burned Area (Burned Area) Monthly L3 Global 500 m SIN Grid
mt_dates("MCD64A1", traps.pts$Latitude[1], traps.pts$Longitude[1])

traps.pts <- traps.pts[, c("treatment", "Latitude", "Longitude")]

SGT.fire <- mt_subset(
  product = "MCD64A1",
  lat = -10.938,
  lon = -46.756,
  band = "Burn_Date",
  start = "2000-01-01",
  end = "2022-07-31",
  site_name = "SGT",
  km_lr = 60,
  km_ab = 60,
  internal = T,
  progress = T
)

head(SGT.fire)

saveRDS(SGT.fire, "Output/SGT_fire.rds")
SGT.fire <- readRDS("Output/SGT_fire.rds")

SGT.fire_r <- mt_to_terra(df = SGT.fire, reproject = T)
SGT.fire_r
names(SGT.fire_r)
plot(SGT.fire_r[[9]])
plot(traps.pts.sp, add = T)
plot(traps.pts.sp$geometry)

terra::writeRaster(SGT.fire_r, "Output/SGT_Fire.tif", overwrite = T)
SGT.fire_r <- rast("Output/SGT_Fire.tif")

fire.traps.wide <- data.frame(
  st_coordinates(traps.pts.sp),
  plot = traps.pts.sp$point,
  regime = traps.pts.sp$treatment,
  fire = terra::extract(SGT.fire_r, traps.pts.sp)
)

fire.traps.df <- data.frame(
  long = rep(fire.traps.wide$X, 261),
  lat = rep(fire.traps.wide$Y, 261),
  plot = rep(fire.traps.wide$plot, 261),
  regime = rep(fire.traps.wide$regime, 261),
  burn_date = c(unlist(c(fire.traps.wide[, 6:266])))
)

fire.traps.df$year <- as.integer(substr(row.names(fire.traps.df), 6, 9))
fire.traps.df$month <- as.integer(c(substr(row.names(fire.traps.df), 11, 14)))

fire.traps.df$fire <- ifelse(fire.traps.df$burn_date > 1, 1, 0)
fire.traps.df$ID <- as.integer(rep(row.names(fire.traps.wide), 261))


library(lubridate)

TSLF_function <- function(df) {
  fire.time <- seq.Date(
    as.Date("2000/11/01"),
    as.Date("2022/07/30"),
    by = "month"
  )
  # fire.time <- fire.time[-c(257:259)]
  new.df <- as.data.frame(matrix(
    NA,
    nrow = length(fire.time),
    ncol = length(unique(df$ID))
  ))
  for (i in unique(df$ID)) {
    df_index <- dplyr::filter(df, ID == i)
    df_index$time <- fire.time
    last_event_index <- cumsum(df_index$fire) + 1
    last_event_index <- c(1, last_event_index[1:length(last_event_index) - 1])
    TSLF <- c(as.Date(NA), df_index[which(df_index$fire == 1), "time"])[
      last_event_index
    ]
    new.df[, as.integer(i)] <- (fire.time - TSLF) / 30
  }
  new.df$time <- fire.time
  new.df
}


TSLF.df <- TSLF_function(fire.traps.df)

TSLF.longdf <- data.frame(
  TSLF = c(unlist(c(TSLF.df[, 1:220]))),
  time = rep(TSLF.df$time, 220),
  ID = rep(1:220, each = 261)
)
TSLF.longdf$month <- as.integer(month(TSLF.longdf$time))
TSLF.longdf$year <- as.integer(year(TSLF.longdf$time))

summary(TSLF.longdf)

library(dplyr)
fire.traps.df <- full_join(fire.traps.df, TSLF.longdf)

# ID, year, month, and fire (0 or 1)

#--- Step 1: Isolate all fire events and create a proper date column ---
fire_events <- fire.traps.df %>%
  filter(fire == 1) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>%
  arrange(ID, date) # IMPORTANT: Ensure dates are in order for each site

#--- Step 2: Calculate the intervals between consecutive fires ---
fire_intervals <- fire_events %>%
  group_by(ID) %>%
  mutate(
    # Get the date of the previous fire for that ID
    previous_fire_date = lag(date),

    # Calculate the interval in days
    interval_days = as.numeric(date - previous_fire_date)
  )

#--- Step 3: Calculate the mean interval (MFRI) for each location ---
MFRI_results <- fire_intervals %>%
  group_by(ID) %>%
  summarise(
    # Calculate the mean of the intervals, removing the NA from the first event
    MFRI_days = mean(interval_days, na.rm = TRUE),

    # Also useful to count the number of fires
    n_fires = n()
  ) %>%
  # Convert to more intuitive units (e.g., months or years)
  mutate(
    MFRI_months = MFRI_days / 30.44, # Average days in a month
    MFRI_years = MFRI_days / 365.25
  )

# View the first few rows of the result
print(head(MFRI_results))

write.csv(fire.traps.df, "Data/fire_traps_df.csv")

fire.regimes.df <- fire.traps.df %>%
  group_by(ID, month) %>%
  summarise(
    fire = sum(fire, na.rm = T),
    MeanTSLF = mean(TSLF, na.rm = T)
  )

fire.regimes.df$weight <- rep(
  c(rep(1, 5), rep(2, 2), rep(3, 3), rep(1, 2)),
  220
)
fire.regimes.df$severity <- fire.regimes.df$fire * fire.regimes.df$weight

write.csv(fire.regimes.df, "Data/fire_regimes_months_df.csv")

fire.regimes.traps <- fire.regimes.df %>%
  group_by(ID) %>%
  summarise(
    freq = sum(fire, na.rm = T),
    MeanTSLF = mean(MeanTSLF, na.rm = T),
    severity = mean(severity, na.rm = T),
  )

fire.regimes.traps$MFRI <- MFRI_results$MFRI_months

fire.regimes.traps$plot <- traps.pts.sp$point
fire.regimes.traps$trap <- c(
  paste0(traps.pts.sp$fieldtrip[1:172], traps.pts.sp$trap[1:172]),
  traps.pts.sp$treatment[173:220]
)


fire.regimes.traps$treatment <- traps.pts.sp$treatment

fire.regimes.traps$TSLF <- traps.pts.sp$tslf

# Save data
write.csv(fire.regimes.traps, "Data/fire_regimes_traps_df.csv")

# Simple plots
plot(severity ~ freq, data = fire.regimes.traps)
plot(severity ~ TSLF, data = fire.regimes.traps)
plot(severity ~ MeanTSLF, data = fire.regimes.traps)
plot(MeanTSLF ~ freq, data = fire.regimes.traps)
plot(TSLF ~ freq, data = fire.regimes.traps)
plot(MeanTSLF ~ TSLF, data = fire.regimes.traps)

plot(MFRI ~ MeanTSLF, data = fire.regimes.traps)
plot(MFRI ~ freq, data = fire.regimes.traps)
plot(severity ~ MFRI, data = fire.regimes.traps)
GGally::ggpairs(fire.regimes.traps[, c(2:5, 9)])

# Check for collinear variables
usdm::vifstep(
  as.data.frame(fire.regimes.traps[, c(2:5, 9)]),
  th = 2,
  keep = "severity"
)
