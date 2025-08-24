rm(list = ls())


# Load packages -----------------------------------------------------------

library(tidyverse)
library(missForest)
library(Boruta)
library(GGally)
library(viridis)
library(ggfortify)
library(gamm4)
library(ggplot2)


# Read environmental data -------------------------------------------------

env.vars <- readRDS("env_vars.rds")
env.vars$plot <- factor(env.vars$plot, levels = c(
  "A1",
  "A2",
  "A3",
  "A4"
))
env.vars <- as.data.frame(env.vars)

env.vars <- do.call(
  data.frame, # Replace Inf in data by NA
  lapply(
    env.vars,
    function(x) replace(x, is.infinite(x), NA)
  )
)

summary(env.vars)


# Ecophysiological variables ----------------------------------------------

## Critical temperatures ---------------------------------------------------

# Load data
ct.sgt <- read.table("CT_Data.txt", h = T)

# Filter by species
ct.Ajalapensis <- ct.sgt[ct.sgt$Species == "Ameivula_jalapensis", ]
ct.Ajalapensis

table(ct.Ajalapensis$Sex)

## Locomotor performance ---------------------------------------------------
# Load data
loc.perf <- read.table("loc_perf_SGT.txt", h = T)
str(loc.perf)
head(loc.perf)

table(loc.perf$sp, loc.perf$SGT)

summary(abs(loc.perf$v))
summary(abs(loc.perf$a))
boxplot(abs(loc.perf$v))

# Get the maximum speed per lizard run
data <- aggregate(abs(loc.perf$v),
  by = list(
    loc.perf$sp,
    loc.perf$SGT,
    loc.perf$Temp,
    loc.perf$Run,
    loc.perf$Sex,
    loc.perf$SVL
  ),
  FUN = max, na.rm = T
)

head(data)

# Rename columns
names(data) <- c(
  "species", "SGT", "temp",
  "run", "sex", "SVL",
  "Veloc"
)

head(data)

# Visualize variation among species
boxplot(Veloc ~ species, data = data)

# Organize data
data.ctmin <- data.frame(
  "species" = ct.sgt$Species,
  "SGT" = ct.sgt$ID,
  "temp" = ct.sgt$Ctmin,
  "run" = c(rep(NA, length(ct.sgt$Ctmin))),
  "sex" = ct.sgt$Sex,
  "SVL" = ct.sgt$SVL,
  "Veloc" = c(rep(0, length(ct.sgt$Ctmin)))
)

data.ctmax <- data.frame(
  "species" = ct.sgt$Species,
  "SGT" = ct.sgt$ID,
  "temp" = ct.sgt$Ctmax,
  "run" = c(rep(NA, length(ct.sgt$Ctmin))),
  "sex" = ct.sgt$Sex,
  "SVL" = ct.sgt$SVL,
  "Veloc" = c(rep(0, length(ct.sgt$Ctmin)))
)

data.complete <- rbind(
  data,
  data.ctmin,
  data.ctmax
)

head(data.complete)
tail(data.complete)
table(data.complete$species)

# Standardize species' names
data.complete$species[data.complete$species == "A_ameiva"] <- "Ameiva_ameiva"
data.complete$species[data.complete$species == "A_jalapensis"] <- "Ameivula_jalapensis"
data.complete$species[data.complete$species == "B_heathi"] <- "Brasiliscincus_heathi"
data.complete$species[data.complete$species == "B_oxyrhina"] <- "Bachia_oxyrhina"
data.complete$species[data.complete$species == "B_oyrhina"] <- "Bachia_oxyrhina"
data.complete$species[data.complete$species == "C_nigropunctatum"] <- "Copeoglossum_nigropunctatum"
data.complete$species[data.complete$species == "H_brasilianus"] <- "Hemidactylus_brasilianus"
data.complete$species[data.complete$species == "T_oreadicus"] <- "Tropidurus_oreadicus"
data.complete$species[data.complete$species == "V_savanicola"] <- "Vanzosaura_savanicola"

table(data.complete$species)

# Thermal locomotor performance curves ------------------------------------

# Ameivula jalapensis
data.Ajalapensis <- data.complete[data.complete$species == "Ameivula_jalapensis", ]

data.Ajalapensis$sex <- as.factor(data.Ajalapensis$sex)
summary(data.Ajalapensis$sex)

# Generalized additive mixed effects model (GAMM), considering individual as random factor
m.Ajalapensis.sprint <- gamm4(Veloc ~ t2(temp),
  random = ~ (1 | SGT),
  data = data.complete[data.complete$species == "Ameivula_jalapensis", ]
)
summary(m.Ajalapensis.sprint$gam)

# Fit GAMM by sex
m.Ajalapensis.sex.sprint <- gamm4(Veloc ~ t2(temp, by = sex),
  random = ~ (1 | SGT),
  data = data.Ajalapensis
)

summary(m.Ajalapensis.sex.sprint$gam)
plot(m.Ajalapensis.sex.sprint$gam)

# Compare models
AIC(m.Ajalapensis.sprint$mer)
AIC(m.Ajalapensis.sex.sprint$mer)

# Refit removing NAs
m.Ajalapensis.sprint.refit <- gamm4(Veloc ~ t2(temp),
  random = ~ (1 | SGT),
  data = data.Ajalapensis[!is.na(data.Ajalapensis$sex), ]
)

# Compare models with likelihood ratio test (LRT)
anova(m.Ajalapensis.sprint.refit$mer, m.Ajalapensis.sex.sprint$mer)

# Plot curve
plot(m.Ajalapensis.sprint$gam,
  main = expression(italic("Ameivula jalapensis")),
  ylab = "Desempenho locomotor", xlab = "Temperatura (°C)"
)

# Predict to make a better plot
preddata_tpc <- data.frame(temp = seq(10, 50, 0.1))

pred_tpc_Ajalapensis <- predict(m.Ajalapensis.sprint$gam,
  preddata_tpc,
  se.fit = T
)

# Merge predictions and SEs in data.frame
preddata_tpc <- rbind(preddata_tpc)
preddata_tpc$predicted <- c(pred_tpc_Ajalapensis$fit)
preddata_tpc$se <- c(pred_tpc_Ajalapensis$se.fit)

# Plot thermal performance curve
ggplot(preddata_tpc, aes(x = temp, y = predicted)) +
  geom_line(color = "blue", linewidth = 1.5) +
  geom_ribbon(aes(ymin = predicted - se, ymax = predicted + se), linetype = 3, alpha = .2) +
  geom_point(
    data = data.complete[data.complete$species == "Ameivula_jalapensis", ],
    aes(x = temp, y = Veloc), alpha = 0.5
  ) +
  labs(y = "Predicted Speed (m/s)", x = "Temperature (°C)") +
  lims(x = c(15, 45), y = c(0, 1.8)) +
  theme(plot.title = element_text(hjust = 0.5))

# Optimal temperature
preddata_tpc[preddata_tpc$predicted == max(pred_tpc_Ajalapensis$fit), ]

# 31.7 for A. jalapensis

# Predict locomotor performance with microclimatic data
microclim.SGT <- readRDS("microclimate_EESGT.rds")

env.var.hour <- microclim.SGT %>%
  group_by(plot, trap.int, fieldtrip, hour, day, month, year) %>%
  summarise(
    tmed = mean(temp, na.rm = T),
    rhmed = mean(rh, na.rm = T)
  )

env.var.hour$Ajalapensis_perf <- predict.gam(m.Ajalapensis.sprint$gam,
  newdata = data.frame(temp = env.var.hour$tmed)
)

# Save object
saveRDS(m.Ajalapensis.sprint, "tpc_Ajalapensis.rds")


# Preferred temperatures --------------------------------------------------
# Load data
tpref_SGT <- read.table("Data_Tpref_SGT.txt", h = T)

# Filter A. jalapensis data
tpref_Ajalapensis <- dplyr::filter(tpref_SGT, sp == "A_jalapensis")
table(tpref_Ajalapensis$SGT, tpref_Ajalapensis$sp)
nrow(table(tpref_Ajalapensis$SGT, tpref_Ajalapensis$sp))

# Function to calculate hours of activity based on quantiles
hvtFUN <- function(temp.envr, temp.lab, quantiles, radiation) {
  vtmin <- quantile(temp.lab, quantiles[1], na.rm = TRUE)
  vtmax <- quantile(temp.lab, quantiles[2], na.rm = TRUE)
  hv <- ifelse(temp.envr > vtmin & temp.envr < vtmax, 1, 0)
  hv[radiation == 0] <- 0
  hv
}

# Calculate hours of activity considering 50% and 90% percentiles
(tpref90_Ajalapensis <- quantile(tpref_Ajalapensis$temp, c(0.05, 0.95), na.rm = T))
(tpref50_Ajalapensis <- quantile(tpref_Ajalapensis$temp, c(0.25, 0.75), na.rm = T))

# Save objects
saveRDS(
  tpref90_Ajalapensis,
  "tpref90_Ajalapensis.rds"
)

saveRDS(
  tpref50_Ajalapensis,
  "tpref50_Ajalapensis.rds"
)

# Predict with microclimatic data
env.var.hour$Ajalapensis_ha90 <- hvtFUN(
  env.var.hour$tmed,
  tpref_Ajalapensis$temp,
  c(0.05, 0.95),
  rep(1, nrow(env.var.hour))
)

env.var.hour$Ajalapensis_ha50 <- hvtFUN(
  env.var.hour$tmed,
  tpref_Ajalapensis$temp,
  c(0.25, 0.75),
  rep(1, nrow(env.var.hour))
)

# Summarize by day
ecophys.day <-
  env.var.hour %>%
  dplyr::group_by(plot, trap.int, fieldtrip, day, month, year) %>%
  dplyr::summarise(
    Ajalapensis_perf = mean(Ajalapensis_perf),
    Ajalapensis_ha50 = sum(Ajalapensis_ha50),
    Ajalapensis_ha90 = sum(Ajalapensis_ha90)
  )

summary(ecophys.day)
pts.traps <- read.table("Points_Traps.txt", h = T)
pts.traps.SGT <- pts.traps[pts.traps$local == "EESGT", ]

# Consider the day length variation as a threshold (diurnal species)
daylength <- geosphere::daylength(
  lat = mean(pts.traps.SGT$lat),
  doy = yday(paste(ecophys.day$year,
    sprintf("%02d", as.numeric(ecophys.day$month)),
    sprintf("%02d", as.numeric(ecophys.day$day)),
    sep = "-"
  ))
)

ecophys.day$Ajalapensis_ha50 <- ifelse(ecophys.day$Ajalapensis_ha50 > daylength,
  daylength,
  ecophys.day$Ajalapensis_ha50
)

ecophys.day$Ajalapensis_ha90 <- ifelse(ecophys.day$Ajalapensis_ha90 > daylength,
  daylength,
  ecophys.day$Ajalapensis_ha90
)

# Summarize by month
ecophys.month <-
  ecophys.day %>%
  dplyr::group_by(plot, trap.int, fieldtrip) %>%
  dplyr::summarise(
    Ajalapensis_perf = mean(Ajalapensis_perf),
    Ajalapensis_ha50 = mean(Ajalapensis_ha50),
    Ajalapensis_ha90 = mean(Ajalapensis_ha90)
  )

ecophys.month$plot <- factor(ecophys.month$plot, levels = c(
  "A1",
  "A2",
  "A3",
  "A4"
))

# Merge environmental variables with ecophysiological variables
env.vars <- left_join(env.vars, ecophys.month, by = c("plot", "trap.int", "fieldtrip"))
summary(env.vars)

# Open cluster to speed up computation
cl <- parallel::makeCluster(ncol(env.vars) - 1,
  setup_timeout = 0.5
)
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()

env.vars$plot.int <- as.integer(env.vars$plot)

# Imputation with missforest
env.vars.imputed <- missForest::missForest(as.matrix(env.vars[, -1]),
  verbose = T,
  variablewise = T,
  maxiter = 100,
  ntree = 1000,
  parallelize = "variables"
)
env.vars.imputed$OOBerror

names(env.vars[, -1])[env.vars.imputed$OOBerror > 1] # Variables with error higher than 1
names(env.vars[, -1])[env.vars.imputed$OOBerror > 5] # Variables with error higher than 5

env.vars.imputed.df <- as.data.frame(env.vars.imputed$ximp)
env.vars.imputed.df <- cbind(env.vars[, 1], env.vars.imputed.df)

names(env.vars.imputed.df) <- c("plot", names(env.vars.imputed.df)[-1])
summary(env.vars.imputed.df)

# Save object
saveRDS(env.vars.imputed.df, "env_vars_imputed_df.rds")

# Which variables better explain the differences among plots?
set.seed(123)
Borutaplot <- Boruta(plot ~ .,
  data = env.vars.imputed.df[, -c(2, 3, 15:18)], doTrace = 2,
  maxRuns = 500
)

Borutaplot
saveRDS(Borutaplot, "Borutaplot.rds")
Borutaplot <- readRDS("Borutaplot.rds")
plot(Borutaplot, las = 2)

windows(height = 8, width = 12)
plot(Borutaplot, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(Borutaplot$ImpHistory), function(i) {
  Borutaplot$ImpHistory[is.finite(Borutaplot$ImpHistory[, i]), i]
})
names(lz) <- colnames(Borutaplot$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1, las = 2, labels = names(Labels),
  at = 1:ncol(Borutaplot$ImpHistory), cex.axis = 0.7
)

plotImpHistory(Borutaplot)
attStats(Borutaplot)
TentativeRoughFix(Borutaplot)

# Pairwise plot
ggpairs.env <- ggpairs(env.vars.imputed.df[env.vars.imputed.df$pground > 70, -c(2, 3, 18)])
quartz(height = 8, width = 12)
ggpairs.env

# Other plots
ggplot(
  env.vars.imputed.df,
  aes(x = plot, y = VARI.all)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Vegetation index")

ggplot(
  env.vars.imputed.df,
  aes(x = plot, y = tree.density)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Tree density")

# Organize fire info
# Time since last fire (TSLF)
fire.traps.df <- read.csv("fire_traps_df.csv")
fire.traps.df <- filter(fire.traps.df, ID > 172)
fire.traps.df <- fire.traps.df[, -1]


fire.traps.df$trap <- as.factor(fire.traps.df$regime)
env.vars.imputed.df$trap <- paste(env.vars.imputed.df$plot,
  env.vars.imputed.df$trap.int,
  sep = "_"
)
env.vars.imputed.df$trap <- factor(env.vars.imputed.df$trap,
  levels = levels(fire.traps.df$trap)
)
levels(fire.traps.df$trap) == levels(env.vars.imputed.df$trap)
env.vars.imputed.df <- env.vars.imputed.df[order(env.vars.imputed.df$fieldtrip), ]
env.vars.imputed.df$year <- rep(c(2021, 2021, 2022, 2022), each = 48)
env.vars.imputed.df$month <- rep(c(3, 7, 2, 6), each = 48)

env.data <- left_join(env.vars.imputed.df, fire.traps.df, by = c(
  "plot",
  "trap",
  "year",
  "month"
))
head(env.data)
summary(env.data)
tapply(env.data$TSLF, INDEX = list(env.data$trap, env.data$fieldtrip), FUN = min)

pal.plot <- turbo(4)

# MODIS made some errors about fires during our sampling. We will correct with our field information
env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 1] <- min(env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 1])
env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 2] <- min(env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 2])
env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 3] <- min(env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 3])
env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 4] <- min(env.data$TSLF[env.data$plot == "A2" & env.data$fieldtrip == 4])

# Last burn in A1 was in June 2020
env.data$TSLF[env.data$plot == "A1" & env.data$fieldtrip == 1] <- 8
env.data$TSLF[env.data$plot == "A1" & env.data$fieldtrip == 2] <- 8 + 3
env.data$TSLF[env.data$plot == "A1" & env.data$fieldtrip == 3] <- 8 + 12
env.data$TSLF[env.data$plot == "A1" & env.data$fieldtrip == 4] <- 8 + 15

# Last burn in A3 was in May 2020
env.data$TSLF[env.data$plot == "A3" & env.data$fieldtrip == 1] <- 9
env.data$TSLF[env.data$plot == "A3" & env.data$fieldtrip == 2] <- 9 + 3
env.data$TSLF[env.data$plot == "A3" & env.data$fieldtrip == 3] <- 9 + 12
env.data$TSLF[env.data$plot == "A3" & env.data$fieldtrip == 4] <- 1

# Last burn in A3 was in May 2020
env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 1] <- min(env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 1])
env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 2] <- min(env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 1]) + 3
env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 3] <- min(env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 1]) + 12
env.data$TSLF[env.data$plot == "A4" & env.data$fieldtrip == 4] <- 0.5

# Fire severity index
# Fires in mid dry season has weight = 2
# Fires in late dry season has weight = 3
# Other fires = 1
fire.regimes <- read.csv("fire_regimes_traps_df.csv")
seq.fire <- 173:220
(severity.A2 <- rep(mean(fire.regimes$severity[seq.fire[1:12]]), 48))
(severity.A1 <- rep(mean(fire.regimes$severity[seq.fire[13:24]]), 48))
(severity.A3 <- rep(mean(fire.regimes$severity[seq.fire[25:36]]), 48))
(severity.A4 <- rep(mean(fire.regimes$severity[seq.fire[37:48]]), 48))


env.data <- env.data[order(env.vars.imputed.df$plot), ]
env.data$plot

env.data$severity <- c(severity.A1, severity.A2, severity.A3, severity.A4)

# Some plots to see relationship of TSLF with environmental data
plot(t.med ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(t.max ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(t.max.abs ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(rh.max.abs ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(rh.min.abs ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(VARI.all ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(tree.density ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])
plot(zentropy ~ TSLF, data = env.data, col = pal.plot[as.numeric(as.factor(env.data$plot))])

# Principal component analysis
pca.env.all <- prcomp(
  env.data[, c(
    "TSLF", "t.med", "t.max", "t.max.abs",
    "rh.min.abs", "rh.max.abs",
    "VARI.all", "zentropy", "zkurt", "pzabovezmean", "pground", "tree.density"
  )],
  scale = T
)
summary(pca.env.all)

pca.env.all
plot(pca.env.all)

# windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c(
  "A1",
  "A2",
  "A3",
  "A4"
))
env.data$fieldtrip <- as.factor(env.data$fieldtrip)
windows(8, 8)
plot.pca.env.all <- autoplot(pca.env.all,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3
)
plot.pca.env.all

# Just with less correlated variables
pca.env <- prcomp(
  env.data[, c(
    "TSLF", "t.med", "t.max", "t.max.abs",
    "rh.min.abs", "rh.max.abs",
    "VARI.all", "zentropy", "tree.density"
  )],
  scale = T
)
summary(pca.env)

pca.env
plot(pca.env)

# windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c(
  "A1",
  "A2",
  "A3",
  "A4"
))
env.data$fieldtrip <- as.factor(env.data$fieldtrip)
windows(8, 8)
plot.pca.env <- autoplot(pca.env,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3
)
plot.pca.env

windows(height = 12, width = 8)
gridExtra::grid.arrange(plot.pca.env.all, plot.pca.env)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)
autoplot(pca.env,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3
)

windows(height = 8, width = 10)
autoplot(pca.env,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3, x = 3, y = 4
)

# Including ecophysiological variables
pca.env.ecophys <- prcomp(
  env.data[, c(
    "TSLF", "t.med", "t.max", "t.max.abs",
    "rh.min.abs", "rh.max.abs",
    "VARI.all", "zentropy", "tree.density",
    "Ajalapensis_perf", "Ajalapensis_ha90"
  )],
  scale = T
)
summary(pca.env.ecophys)

pca.env
plot(pca.env.ecophys)

# windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c(
  "A1",
  "A2",
  "A3",
  "A4"
))
env.data$fieldtrip <- as.factor(env.data$fieldtrip)
windows(8, 8)
plot.pca.env.ecophys <- autoplot(pca.env.ecophys,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3
)
plot.pca.env.ecophys

windows(height = 12, width = 8)
quartz(height = 12, width = 8)

gridExtra::grid.arrange(plot.pca.env.ecophys, plot.pca.env)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)
autoplot(pca.env.ecophys,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3
)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)

autoplot(pca.env.ecophys,
  data = env.data, colour = "plot", shape = "fieldtrip", frame = T,
  loadings = T, loadings.label = TRUE, loadings.label.size = 3, x = 3, y = 4
)

env.data <- env.data[, c(
  "plot", "trap", "fieldtrip", "severity", "TSLF", "t.med", "t.max", "t.max.abs",
  "rh.min.abs", "rh.max.abs", "VARI.all", "zentropy", "tree.density",
  "Ajalapensis_perf", "Ajalapensis_ha90"
)]

# Other plots
quartz(height = 8, width = 8)
ggplot(
  env.data,
  aes(x = plot, y = Ajalapensis_perf)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Locomotor performance (m/s)")

quartz(height = 8, width = 8)
ggplot(
  env.data,
  aes(x = plot, y = Ajalapensis_ha90)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Hours of activity")

quartz(height = 8, width = 12)
ggpairs.env <- ggpairs(env.data[, -c(2:4)])
ggpairs.env

# Save object
saveRDS(env.data, "env_data_SGT.rds")
