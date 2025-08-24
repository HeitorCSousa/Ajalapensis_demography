library(MODISTools)
library(sf)
library(raster)
library(terra)
library(sp)
library(tidyverse)


# Read and plot points coordinates
arms.pts <- read.table("Pontos_Arms_Bruna_Heitor.txt", h = T)
arms.pts

arms.pts.sp <- arms.pts

coordinates(arms.pts.sp) <- c("Longitude", "Latitude")
proj4string(arms.pts.sp) <- CRS("+proj=longlat +datum=WGS84")
arms.pts.sp
arms.pts.sp <- st_as_sf(arms.pts.sp)

mt_products()
# MODIS/Terra+Aqua Burned Area (Burned Area) Monthly L3 Global 500 m SIN Grid
mt_bands(product = "MCD64A1")

# MODIS/Terra+Aqua Burned Area (Burned Area) Monthly L3 Global 500 m SIN Grid
mt_dates("MCD64A1", arms.pts$Latitude[1], arms.pts$Longitude[1])

arms.pts <- arms.pts[, c("tratamento", "Latitude", "Longitude")]

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

saveRDS(SGT.fire, "SGT.fire.rds")
SGT.fire <- readRDS("SGT.fire.rds")

SGT.fire_r <- mt_to_terra(df = SGT.fire, reproject = T)
SGT.fire_r
names(SGT.fire_r)
plot(SGT.fire_r[[9]])
plot(arms.pts.sp, add = T)
plot(arms.pts.sp$geometry)

terra::writeRaster(SGT.fire_r, "SGT_Fire.tif", overwrite = T)
SGT.fire_r <- rast("SGT_Fire.tif")

fire.arms.wide <- data.frame(st_coordinates(arms.pts.sp),
  plot = arms.pts.sp$ponto,
  regime = arms.pts.sp$tratamento,
  fire = terra::extract(SGT.fire_r, arms.pts.sp)
)

fire.arms.df <- data.frame(
  long = rep(fire.arms.wide$X, 261),
  lat = rep(fire.arms.wide$Y, 261),
  plot = rep(fire.arms.wide$plot, 261),
  regime = rep(fire.arms.wide$regime, 261),
  burn_date = c(unlist(c(fire.arms.wide[, 6:266])))
)

fire.arms.df$year <- as.integer(substr(row.names(fire.arms.df), 6, 9))
fire.arms.df$month <- as.integer(c(substr(row.names(fire.arms.df), 11, 14)))

fire.arms.df$fire <- ifelse(fire.arms.df$burn_date > 1, 1, 0)
fire.arms.df$ID <- as.integer(rep(row.names(fire.arms.wide), 261))


library(lubridate)

TSLF_function <- function(df) {
  fire.time <- seq.Date(as.Date("2000/11/01"), as.Date("2022/07/30"), by = "month")
  # fire.time <- fire.time[-c(257:259)]
  new.df <- as.data.frame(matrix(NA, nrow = length(fire.time), ncol = length(unique(df$ID))))
  for (i in unique(df$ID)) {
    df_index <- dplyr::filter(df, ID == i)
    df_index$time <- fire.time
    last_event_index <- cumsum(df_index$fire) + 1
    last_event_index <- c(1, last_event_index[1:length(last_event_index) - 1])
    TSLF <- c(as.Date(NA), df_index[which(df_index$fire == 1), "time"])[last_event_index]
    new.df[, as.integer(i)] <- (fire.time - TSLF) / 30
  }
  new.df$time <- fire.time
  new.df
}


TSLF.df <- TSLF_function(fire.arms.df)

TSLF.longdf <- data.frame(
  TSLF = c(unlist(c(TSLF.df[, 1:220]))),
  time = rep(TSLF.df$time, 220),
  ID = rep(1:220, each = 261)
)
TSLF.longdf$month <- as.integer(month(TSLF.longdf$time))
TSLF.longdf$year <- as.integer(year(TSLF.longdf$time))

summary(TSLF.longdf)

library(dplyr)
fire.arms.df <- full_join(fire.arms.df, TSLF.longdf)

# ID, year, month, and fire (0 or 1)

#--- Step 1: Isolate all fire events and create a proper date column ---
fire_events <- fire.arms.df %>%
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

write.csv(fire.arms.df, "fire_arms_df.csv")

fire.regimes.df <- fire.arms.df %>%
  group_by(ID, month) %>%
  summarise(
    fire = sum(fire, na.rm = T),
    MeanTSLF = mean(TSLF, na.rm = T)
  )

fire.regimes.df$weight <- rep(c(rep(1, 5), rep(2, 2), rep(3, 3), rep(1, 2)), 220)
fire.regimes.df$severity <- fire.regimes.df$fire * fire.regimes.df$weight

write.csv(fire.regimes.df, "fire_regimes_months_df.csv")

fire.regimes.arms <- fire.regimes.df %>%
  group_by(ID) %>%
  summarise(
    freq = sum(fire, na.rm = T),
    MeanTSLF = mean(MeanTSLF, na.rm = T),
    severity = mean(severity, na.rm = T),
  )

fire.regimes.arms$MFRI <- MFRI_results$MFRI_months

fire.regimes.arms$plot <- arms.pts.sp$ponto
fire.regimes.arms$trap <- c(
  paste0(arms.pts.sp$campanha[1:172], arms.pts.sp$armadilha[1:172]),
  arms.pts.sp$tratamento[173:220]
)


fire.regimes.arms$treatment <- arms.pts.sp$tratamento

fire.regimes.arms$TSLF <- arms.pts.sp$tuf


write.csv(fire.regimes.arms, "fire_regimes_arms_df.csv")

plot(severity ~ freq, data = fire.regimes.arms)
plot(severity ~ TSLF, data = fire.regimes.arms)
plot(severity ~ MeanTSLF, data = fire.regimes.arms)
plot(MeanTSLF ~ freq, data = fire.regimes.arms)
plot(TSLF ~ freq, data = fire.regimes.arms)
plot(MeanTSLF ~ TSLF, data = fire.regimes.arms)

plot(MFRI ~ MeanTSLF, data = fire.regimes.arms)
plot(MFRI ~ freq, data = fire.regimes.arms)
plot(severity ~ MFRI, data = fire.regimes.arms)
GGally::ggpairs(fire.regimes.arms[, c(2:5, 9)])
usdm::vifstep(as.data.frame(fire.regimes.arms[, c(2:5, 9)]), th = 2, keep = "severity")
