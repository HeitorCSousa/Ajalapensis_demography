# Set working directory and load packages ---------------------------------
# To load a pre-saved RData session with all objects and models,
# uncomment and run the line below.
# load("Ajalapensis_demography.RData")

# To run the script from a clean slate, clear all objects from the
# current environment.
rm(list = ls())


# Load packages -----------------------------------------------------------

library(ggplot2)
library(BaSTA)
library(openCR)
library(sf)
library(sp)
library(dplyr)
library(lubridate)
library(viridis)
library(tidyverse)
library(brms)
library(projpred)
library(GGally)
library(emmeans)
library(scales)

# If needed, install INLA by uncommenting the following lines
# install.packages("INLA", repos = c(getOption("repos"),
#   INLA = "https://inla.r-inla-download.org/R/testing"
# ), dep = TRUE)
library(INLA)


# Fixed plots -------------------------------------------------------------

## Load capture data  ------------------------------------------

# Load the main dataset containing all capture records.
data <- read.table("Ameivula_jalapensis_04.txt", h = T)

# Ensure the 'Plot' column is a factor with a specific order for consistent plotting.
data$Plot <- factor(data$Plot, levels = c("A1", "A2", "A3", "A4"))

# Attach the data frame to the search path for easy access to columns.
attach(data)

# Display the first and last few rows to verify data integrity.
head(data)
tail(data)

# Show the structure of the data frame, including column types.
str(data)

## Calculate capture and recapture frequencies -----------------------------

# Total number of captures (including recaptures).
length(SVL)

# Frequency table of captures by sex.
table(Sex)

# Summary of captures per plot.
summary(as.factor(Plot))

# Frequency table of captures by plot and species.
(recap.especies <- table(Plot, Species))

# Total captures and recaptures by plot.
(cap <- table(Plot))

# Captures of new individuals (not recaptures) by plot.
(cap <- table(Plot[Recapture == "N"]))

# Recaptures by plot.
(recap <- table(Plot[Recapture == "Y"]))

# Total number of unique individuals captured.
sum(cap)

# Total number of recaptures.
sum(recap)


# Average recapture rate per individual.
sum(table(Plot[Recapture != "N"])) / length(SVL)

# Maximum number of times a single individual was recaptured.
max(table(VoucherNumber[Recapture == "Y"]))

# Frequency table of recapture counts (e.g., how many lizards were recaptured once, twice, ...).
table((table(data$VoucherNumber[data$Recapture == "Y"])))

# Create a summary table of captures vs. recaptures for each plot.
table(Plot, Recapture)
a <- table(Plot, Recapture)
colnames(a)
colnames(a) <- c("Captures", "Recaptures")

# Reorder rows by plot name for consistency.
a <- (a[order(row.names(a)), ])
print(a)

# Barplot visualizing captures and recaptures side-by-side for each plot.
barplot(
  t(a),
  beside = TRUE,
  legend.text = T,
  args.legend = list(x = "topleft"),
  main = expression(italic("Ameivula jalapensis")),
  ylim = c(0, 120),
  ylab = "Frequency"
)

# Frequency table of captures by field trip and plot.
table(Fieldtrip, Plot)
b <- table(Fieldtrip, Plot)

# Reorder columns by plot name.
b <- (b[, order(colnames(b))])
print(b)

# Barplot visualizing captures per plot for each field trip.
barplot(
  t(b),
  beside = TRUE,
  legend.text = T,
  args.legend = list(x = "topleft"),
  main = expression(italic("Ameivula jalapensis")),
  ylim = c(0, 50),
  ylab = "Captures"
)

# Frequency table of captures by sex and plot.
table(Sex, Plot)
sex <- table(Sex, Plot)

# Reorder columns by plot name.
sex <- (sex[, order(colnames(sex))])

# Reorder rows for a more logical presentation (Male, Female, Indeterminate).
sex <- sex[c(2, 1, 3), ]
print(sex)

# Barplot visualizing sex distribution across plots.
barplot(
  sex,
  legend.text = T,
  args.legend = list(x = "topleft"),
  ylim = c(0, 60),
  ylab = "Frequency",
  cex.axis = 1.0,
  cex.lab = 1.0,
  cex.names = 1.0,
  cex.sub = 1,
  beside = TRUE
)


## SVL (Snout-Vent Length) Temporal and Spatial Variation ----------------------

# Basic summary statistics for SVL.
head(data)
summary(SVL)
sd(SVL, na.rm = T)

# Plot SVL variation across field trips using a half-violin, half-dot plot.
# This visualization shows the distribution, median/quartiles, and raw data points.
ggplot(
  data,
  aes(x = as.factor(Fieldtrip), y = SVL)
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
  geom_hline(yintercept = 40, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Snout-vent length (mm)")

# Plot SVL variation among plots.
ggplot(
  data,
  aes(x = as.factor(Plot), y = SVL)
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
  geom_hline(yintercept = 40, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Snout-vent length (mm)")

# Detach the data frame from the search path.
detach(data)

# Pradel-Jolly-Seber Model ------------------------------------------------

# -- Data Preparation for openCR --
# This section merges capture data with trap coordinates and formats it for analysis.

# Load trap coordinates data.
pts.traps <- read.table("Points_Traps.txt", h = T)

# Filter for the specific locality.
pts.traps.SGT <- pts.traps[pts.traps$locality == "EESGT", ]

# Convert the coordinates from Latitude/Longitude (WGS84) to UTM.
# This is crucial for accurate distance-based calculations in spatial models.
coordinates(pts.traps.SGT) <- c("long", "lat")
proj4string(pts.traps.SGT) <- CRS("+proj=longlat +datum=WGS84")
# UTM Zone 22S
pts.traps.SGT.utm <- spTransform(pts.traps.SGT, CRS = CRS("+init=epsg:32722"))

# Extract the UTM coordinates and add them back to the data frame.
pts.traps.SGT <- pts.traps[pts.traps$locality == "EESGT", ]
pts.traps.SGT$X <- coordinates(pts.traps.SGT.utm)[, 1]
pts.traps.SGT$Y <- coordinates(pts.traps.SGT.utm)[, 2]

# Rename columns for clarity.
names(pts.traps.SGT) <- c(
  "Plot",
  "Trap",
  "Lat",
  "Long",
  "Locality",
  "X",
  "Y"
)

# Join the capture data with the trap coordinates.
Ajalapensis.planilha.XY <- left_join(
  data,
  pts.traps.SGT,
  by = c("Plot", "Trap")
)
head(Ajalapensis.planilha.XY)

# -- Create Capture History Object --
# The data needs to be structured by session and occasion for the `openCR` package.
# The following steps manually adjust the 'occasion' based on the specific sampling design.
# This logic is complex and study-specific.

# Subset data for Field Trip 1 and adjust occasions.
capts.1 <- data.frame(
  Session = Ajalapensis.planilha.XY$Fieldtrip[
    Ajalapensis.planilha.XY$Fieldtrip == 1
  ],
  ID = Ajalapensis.planilha.XY$VoucherNumber[
    Ajalapensis.planilha.XY$Fieldtrip == 1
  ],
  occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Date[
    Ajalapensis.planilha.XY$Fieldtrip == 1
  ])),
  X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Fieldtrip == 1],
  Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Fieldtrip == 1],
  date = Ajalapensis.planilha.XY$Date[Ajalapensis.planilha.XY$Fieldtrip == 1],
  Plot = Ajalapensis.planilha.XY$Plot[Ajalapensis.planilha.XY$Fieldtrip == 1]
)

capts.1$occasion[4:nrow(capts.1)] <- capts.1$occasion[4:nrow(capts.1)] + 1
capts.1$occasion[19:nrow(capts.1)] <- capts.1$occasion[19:nrow(capts.1)] + 3
capts.1$occasion <- capts.1$occasion + 2

# Subset data for Field Trip 2 and adjust occasions.
capts.2 <- data.frame(
  Session = Ajalapensis.planilha.XY$Fieldtrip[
    Ajalapensis.planilha.XY$Fieldtrip == 2
  ],
  ID = Ajalapensis.planilha.XY$VoucherNumber[
    Ajalapensis.planilha.XY$Fieldtrip == 2
  ],
  occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Date[
    Ajalapensis.planilha.XY$Fieldtrip == 2
  ])),
  X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Fieldtrip == 2],
  Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Fieldtrip == 2],
  date = Ajalapensis.planilha.XY$Date[Ajalapensis.planilha.XY$Fieldtrip == 2],
  Plot = Ajalapensis.planilha.XY$Plot[Ajalapensis.planilha.XY$Fieldtrip == 2]
)

capts.2$occasion[32:nrow(capts.2)] <- capts.2$occasion[32:nrow(capts.2)] + 1

# Subset data for Field Trip 3 and adjust occasions.
capts.3 <- data.frame(
  Session = Ajalapensis.planilha.XY$Fieldtrip[
    Ajalapensis.planilha.XY$Fieldtrip == 3
  ],
  ID = Ajalapensis.planilha.XY$VoucherNumber[
    Ajalapensis.planilha.XY$Fieldtrip == 3
  ],
  occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Date[
    Ajalapensis.planilha.XY$Fieldtrip == 3
  ])),
  X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Fieldtrip == 3],
  Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Fieldtrip == 3],
  date = Ajalapensis.planilha.XY$Date[Ajalapensis.planilha.XY$Fieldtrip == 3],
  Plot = Ajalapensis.planilha.XY$Plot[Ajalapensis.planilha.XY$Fieldtrip == 3]
)

capts.3$occasion[65:nrow(capts.3)] <- capts.3$occasion[65:nrow(capts.3)] + 1

# Subset data for Field Trip 4.
capts.4 <- data.frame(
  Session = Ajalapensis.planilha.XY$Fieldtrip[
    Ajalapensis.planilha.XY$Fieldtrip == 4
  ],
  ID = Ajalapensis.planilha.XY$VoucherNumber[
    Ajalapensis.planilha.XY$Fieldtrip == 4
  ],
  occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Date[
    Ajalapensis.planilha.XY$Fieldtrip == 4
  ])),
  X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Fieldtrip == 4],
  Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Fieldtrip == 4],
  date = Ajalapensis.planilha.XY$Date[Ajalapensis.planilha.XY$Fieldtrip == 4],
  Plot = Ajalapensis.planilha.XY$Plot[Ajalapensis.planilha.XY$Fieldtrip == 4]
)

# Combine all field trips into a single data frame.
capts <- rbind(
  capts.1[, -6],
  capts.2[, -6],
  capts.3[, -6],
  capts.4[, -6]
)

# Remove rows with missing ID or coordinates.
complete <- complete.cases(capts[, c("ID", "X", "Y")])
capts <- droplevels(capts[complete, ])

# Create the trap layout object for `openCR`.
traps.xy <- data.frame(
  trapID = interaction(pts.traps.SGT$Plot, pts.traps.SGT$Trap),
  x = pts.traps.SGT$X,
  y = pts.traps.SGT$Y
)
traps.xy <- read.traps(data = traps.xy[, -1], detector = "multi")

# Create the final capture history object (`capthist`).
Ajalapensis.CHxy <- make.capthist(
  capts,
  traps.xy,
  fmt = "XY",
  noccasions = 15,
  covnames = c("Plot")
)

# Set the time intervals between primary sampling sessions (in months).
intervals(Ajalapensis.CHxy) <- c(112, 245, 84) / 30

# Save capture history object
write.capthist(Ajalapensis.CHxy)

# Creat mesh mask
mesh <- make.mask(traps.xy, spacing = 30, type = "trapbuffer", buffer = 150)

# Summary statistics
summary(mesh)
summary(Ajalapensis.CHxy)
m.array(Ajalapensis.CHxy)
JS.counts(Ajalapensis.CHxy)
JS.direct(Ajalapensis.CHxy)

summary(traps(Ajalapensis.CHxy))

# Visualize the spatio-temporal variation in captures
par(mfrow = c(2, 2))
plot(mesh)
plot(Ajalapensis.CHxy$`1`, tracks = T, add = T)

plot(mesh)
plot(Ajalapensis.CHxy$`2`, tracks = T, add = T)

plot(mesh)
plot(Ajalapensis.CHxy$`3`, tracks = T, add = T)

plot(mesh)
plot(Ajalapensis.CHxy$`4`, tracks = T, add = T)
par(mfrow = c(1, 1))

# -- Fit Non-Spatial Pradel Models --
# These models estimate:
#   phi (φ): Apparent survival probability
#   p: Capture probability
#   f: Per capita recruitment rate

# Model 0: Constant parameters (phi(.), p(.), f(.))
# Assumes survival, capture, and recruitment are constant over time and space.
pradel.fit0 <- openCR.fit(Ajalapensis.CHxy, type = "JSSAfCL")
summary(pradel.fit0)

# Model 1: Time-dependent parameters (phi(t), p(t), f(t))
# Allows parameters to vary between primary sessions (field trips).
pradel.fit1 <- openCR.fit(
  Ajalapensis.CHxy,
  type = "JSSAfCL",
  model = list(phi ~ t, p ~ t, f ~ t),
  method = "Nelder-Mead",
  details = list(control = list(maxit = 5000)),
  start = pradel.fit0
)
summary(pradel.fit1)
p.fit1 <- predict(pradel.fit1, all = T)

# Plot time-dependent survival estimates.
ggplot(p.fit1$phi, aes()) +
  geom_pointrange(
    aes(x = session, y = estimate, ymin = lcl, ymax = ucl),
    size = 1
  ) +
  labs(x = "Occasion", y = "Survival") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Plot time-dependent capture probability estimates.
ggplot(p.fit1$p, aes()) +
  geom_pointrange(
    aes(x = session, y = estimate, ymin = lcl, ymax = ucl),
    size = 1
  ) +
  labs(x = "Occasion", y = "Capture") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Plot time-dependent recruitment estimates.
ggplot(p.fit1$f, aes()) +
  geom_pointrange(
    aes(x = session, y = estimate, ymin = lcl, ymax = ucl),
    size = 1
  ) +
  labs(x = "Occasion", y = "Recruitment") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Model 2: Plot-dependent parameters (phi(Plot), p(Plot), f(Plot))
# Allows parameters to vary between the different fire severity plots.
pradel.fit2 <- openCR.fit(
  Ajalapensis.CHxy,
  type = "JSSAfCL",
  model = list(phi ~ Plot, p ~ Plot, f ~ Plot)
)
summary(pradel.fit2)

# Extract model estimates
p.fit2 <- predict(pradel.fit2, all = T)

# Standardize plot levels
p.fit2$phi$Plot <- factor(p.fit2$phi$Plot, levels = c("A1", "A2", "A3", "A4"))

p.fit2$p$Plot <- factor(p.fit2$p$Plot, levels = c("A1", "A2", "A3", "A4"))

p.fit2$f$Plot <- factor(p.fit2$f$Plot, levels = c("A1", "A2", "A3", "A4"))

# Extract, Subset, and Rename  to make pairwise contrasts
all_coefs_df <- coef(pradel.fit2)
all_vcov <- vcov(pradel.fit2)

phi_names <- grep("^phi", rownames(all_coefs_df), value = TRUE)
phi_coefs <- all_coefs_df[phi_names, "beta"]
names(phi_coefs) <- phi_names
phi_vcov <- all_vcov[phi_names, phi_names]

names(phi_coefs)[names(phi_coefs) == "phi"] <- "(Intercept)"
names(phi_coefs) <- gsub("phi.Plot", "Plot", names(phi_coefs))

# Create the Date for the Grid and Run qdrg()

# Define the levels of the factor
factor_levels <- c("A1", "A2", "A3", "A4")

# Create a data frame containing the factor
grid_data <- data.frame(Plot = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_phi <- qdrg(
  formula = ~Plot,
  data = grid_data, # <-- Add the data argument here
  coef = phi_coefs,
  vcov = phi_vcov
)

# Calculate and View the Pairwise Contrasts
contrasts_phi <- pairs(ref_grid_phi)

print(contrasts_phi)

p_names <- grep("^p", rownames(all_coefs_df), value = TRUE)[1:4]
p_coefs <- all_coefs_df[p_names, "beta"]
names(p_coefs) <- p_names
p_vcov <- all_vcov[p_names, p_names]

names(p_coefs)[names(p_coefs) == "p"] <- "(Intercept)"
names(p_coefs) <- gsub("p.Plot", "Plot", names(p_coefs))

# Create the Data for the Grid and Run qdrg()

# Define the levels of the factor
factor_levels <- c("A1", "A2", "A3", "A4")

# Create a data frame containing the factor
grid_data <- data.frame(Plot = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_p <- qdrg(
  formula = ~Plot,
  data = grid_data,
  coef = p_coefs,
  vcov = p_vcov
)

# Calculate and View the Pairwise Contrasts
contrasts_p <- pairs(ref_grid_p)

print(contrasts_p)

f_names <- grep("^f", rownames(all_coefs_df), value = TRUE)[1:4]
f_coefs <- all_coefs_df[f_names, "beta"]
names(f_coefs) <- f_names
f_vcov <- all_vcov[f_names, f_names]

names(f_coefs)[names(f_coefs) == "f"] <- "(Intercept)"
names(f_coefs) <- gsub("f.Plot", "Plot", names(f_coefs))

# Create the Data for the Grid and Run qdrg()

# Define the levels of the factor
factor_levels <- c("A1", "A2", "A3", "A4")

# Create a data frame containing the factor
grid_data <- data.frame(Plot = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_f <- qdrg(
  formula = ~Plot,
  data = grid_data,
  coef = f_coefs,
  vcov = f_vcov
)

# Calculate and View the Pairwise Contrasts
contrasts_f <- pairs(ref_grid_f)

print(contrasts_f)

# Plot plot-dependent survival estimates.
ggplot(p.fit2$phi[c(1, 5, 9, 13), ], aes()) +
  geom_pointrange(
    aes(x = Plot, y = estimate, ymin = lcl, ymax = ucl, color = Plot),
    size = 1,
    show.legend = F
  ) +
  labs(x = "Fire severity", y = "Survival") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  ) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_phi)

# Plot plot-dependent capture probability estimates.
ggplot(p.fit2$p[c(1, 5, 9, 13), ], aes()) +
  geom_pointrange(
    aes(x = Plot, y = estimate, ymin = lcl, ymax = ucl, color = Plot),
    size = 1,
    show.legend = F
  ) +
  labs(x = "Fire severity", y = "Capture probability") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  ) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_p)

# Plot plot-dependent recruitment estimates.
ggplot(p.fit2$f[c(1, 5, 9, 13), ], aes()) +
  geom_pointrange(
    aes(x = Plot, y = estimate, ymin = lcl, ymax = ucl, color = Plot),
    size = 1,
    show.legend = F
  ) +
  labs(x = "Fire severity", y = "Recruitment") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  ) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_f)

# Model 3: Interaction between Plot and time
# Allows parameters to vary by both plot and time.
pradel.fit3 <- openCR.fit(
  Ajalapensis.CHxy,
  type = "JSSAfCL",
  model = list(
    phi ~ -1 + Plot * t,
    p ~ -1 + Plot * t,
    f ~ -1 + Plot * t
  ),
  start = pradel.fit0
)
summary(pradel.fit3)
predict(pradel.fit3, all = T)

# -- Model Comparison --
# Compare the fitted models using Akaike's Information Criterion (AIC).
# The model with the lowest AIC value is considered the best fit.
AIC(pradel.fit0, pradel.fit1, pradel.fit2, pradel.fit3)

# Convergence problems. Trying to solve them, but they seem unindentifiable
pradel.fit3 <- openCR.fit(
  Ajalapensis.CHxy,
  type = "JSSAfCL",
  model = list(
    phi ~ -1 + Plot * t,
    p ~ -1 + Plot * t,
    f ~ -1 + Plot * t
  ),
  method = "Nelder-Mead",
  details = list(
    control = list(
      maxit = 20000,
      reltol = 1e-6
    )
  ),
  start = pradel.fit0
)
summary(pradel.fit3)
predict(pradel.fit3, all = T)

AIC(pradel.fit0, pradel.fit1, pradel.fit2, pradel.fit3)

# Conclusion: The most complex model is the least parsimonious. Populations vary over time, but not across plots.

# SVL with environment ----------------------------------------------------

# -- Data Preparation --
# Merge capture data with environmental data and scale predictors.
env.data <- readRDS("env_data_SGT.rds")

# Prepare a clean SVL dataset for merging.
svl.data <- data[, c(
  "Fieldtrip",
  "Date",
  "Month",
  "Year",
  "Plot",
  "Trap",
  "VoucherNumber",
  "Recapture",
  "Mass",
  "SVL",
  "TL",
  "TB",
  "Broken",
  "Sex",
  "Eggs",
  "Dead"
)]

names(svl.data) <- c(
  "fieldtrip",
  "date",
  "month",
  "year",
  "plot",
  "trap",
  "ID",
  "recapture",
  "mass",
  "svl",
  "tl",
  "tb",
  "broken",
  "sex",
  "eggs",
  "dead"
)

svl.data$trap <- as.factor(paste(svl.data$plot, svl.data$trap, sep = "_"))
svl.data$fieldtrip <- as.factor(svl.data$fieldtrip)

# Join SVL data with environmental data.
svl.env <- left_join(svl.data, env.data, by = c("fieldtrip", "plot", "trap"))

# Scale all numeric environmental predictors (mean = 0, sd = 1).
# This is important for improving model convergence and comparing effect sizes.
svl.env[, c(18:28)] <- scale(svl.env[, c(18:28)])
head(svl.env)
hist((svl.env$svl))

# -- Full Model and Variable Selection --
# Fit a full model with all potential predictors to investigate their effect on SVL.
mfull.svl <- brm(
  svl ~
    plot +
    TSLF +
    t.med +
    t.max +
    t.max.abs +
    rh.min.abs +
    rh.max.abs +
    VARI.all +
    zentropy +
    tree.density +
    Ajalapensis_perf +
    Ajalapensis_ha90,
  data = svl.env,
  cores = 4
)

summary(mfull.svl)
plot(mfull.svl)
bayes_R2(mfull.svl)

# Perform variable selection using projpred to find the most influential variables.
# First, perform variable selection without cross-validation.
vs <- cv_varsel(mfull.svl, validate_search = F)
summary(vs)
summary(vs, deltas = T)

suggest_size(vs)
ranking(vs)

plot(vs, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.svl)

# Then, perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits <- run_cvfun(
  refm_obj,
  K = 10
)

# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
doParallel::registerDoParallel(5)
# Final cv_varsel() run:
cvvs <- cv_varsel(
  refm_obj,
  cv_method = "kfold",
  cvfits = cv_fits,
  nterms_max = 3,
  parallel = TRUE,
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

summary(cvvs, deltas = T, type = c("mean", "lower", "upper"))


plot(cv_proportions(ranking(cvvs)))

plot(cv_proportions(ranking(cvvs), cumulate = T))

plot(cvvs, deltas = T, stats = "mlpd") + theme_minimal()

# -- Fit Selected Model --
# 'Ajalapensis_perf' is a key predictor.
# We fit a mixed-effects model with random intercepts for the sampling hierarchy
# (traps nested in plots, nested in field trips) to account for non-independence.
msel.svl.re <- brm(
  svl ~ Ajalapensis_perf + (1 | fieldtrip / plot / trap),
  data = svl.env,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.svl.re)

# Calculate Bayesian R-squared.
bayes_R2(msel.svl.re)

# -- Model Diagnostics and Visualization --

# Check trace plots for convergence.
plot(msel.svl.re)

# Posterior predictive check to assess model fit.
pp_check(msel.svl.re)

# Plot the conditional effect of locomotor performance on SVL.
p1 <- conditional_effects(msel.svl.re)

ggplot(p1$Ajalapensis_perf, aes(x = Ajalapensis_perf, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = svl.env,
    aes(x = Ajalapensis_perf, y = svl, color = plot),
    alpha = 0.5,
    size = 2
  ) +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  labs(x = "Locomotor performance", y = "Snout-vent length (mm)") +
  theme_minimal()

# Testing differences among plots
msel.svl.plot.re <- brm(
  svl ~ plot + (1 | fieldtrip / trap),
  data = svl.env,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.svl.plot.re)

# Calculate Bayesian R-squared.
bayes_R2(msel.svl.plot.re)

# -- Model Diagnostics and Visualization --
# Check trace plots for convergence.
plot(msel.svl.plot.re)

# Posterior predictive check to assess model fit.
pp_check(msel.svl.plot.re)

# Conditional effects plot
p1 <- conditional_effects(msel.svl.plot.re)

ggplot(
  svl.env,
  aes(x = plot, y = svl, colour = plot)
) +
  ggdist::stat_halfeye(
    aes(fill = plot),
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA,
    alpha = 0.5
  ) +
  gghalves::geom_half_point(
    aes(colour = plot),
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  geom_pointrange(
    data = p1$plot,
    aes(
      x = plot,
      y = estimate__,
      ymin = lower__,
      ymax = upper__,
      colour = plot
    ),
    size = 1.2,
    linewidth = 1.2
  ) +
  geom_hline(yintercept = 40, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  theme_minimal()

# Pairwise contrasts
emmeans::emmeans(msel.svl.plot.re, pairwise ~ plot)


# Abundance with environment -----------------------------------------------

# -- Data Preparation for Abundance Models --
# This involves creating a dataset with one row per trap per field trip,
# with a column for the total number of captures.
str(svl.env)

svl.env <- svl.env %>%
  group_by(fieldtrip) %>%
  mutate(sampling_day = dense_rank(date)) %>% # Rank the dates within each month
  ungroup()

summary(svl.env$sampling_day)

fire.regimes.traps <- read.csv("fire_regimes_traps_df.csv")

trap_locations <- fire.regimes.traps %>%
  distinct(plot, trap, treatment)

trap_locations$fieldtrip <- substr(trap_locations$trap, 1, 2)
trap_locations$fieldtrip[trap_locations$plot == "A2"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A1"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A3"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A4"] <- NA

trap_locations_heitor <- trap_locations[is.na(trap_locations$fieldtrip), ]
trap_locations_bruna <- trap_locations[!is.na(trap_locations$fieldtrip), ]

trap_locations_heitor <- reshape::expand.grid.df(
  trap_locations_heitor[, -4],
  data.frame(fieldtrip = as.character(c(1:4)))
)

trap_locations <- rbind(trap_locations_bruna, trap_locations_heitor)

# Now, create the full dataset of all possible sampling events
# `crossing` will combine every row from `trap_locations` with every sampling day.
all_sampling_events <- crossing(trap_locations, sampling_day = 1:15)

# Count captures only for existing combinations in your data
capture_counts <- svl.env %>%
  count(fieldtrip, plot, trap, sampling_day, name = "freq")

names(capture_counts) <- c("fieldtrip", "plot", "trap", "sampling_day", "freq")
str(capture_counts)
str(all_sampling_events)

capture_counts$fieldtrip <- as.character(capture_counts$fieldtrip)
capture_counts$plot <- as.character(capture_counts$plot)
capture_counts$trap <- as.character(capture_counts$trap)

# Join the actual captures to the master list of all sampling events
Ajalapensis.captures.day <- all_sampling_events %>%
  left_join(
    capture_counts,
    by = c("fieldtrip", "plot", "trap", "sampling_day")
  ) %>%
  # If a trap on a given day had no captures, its freq is NA. Change that to 0.
  mutate(freq = replace_na(freq, 0))

# View the result
print(Ajalapensis.captures.day)

summary(Ajalapensis.captures.day)
sum(Ajalapensis.captures.day$freq)


Ajalapensis.captures.day <- Ajalapensis.captures.day %>%
  pivot_wider(
    names_from = sampling_day,
    values_from = freq,
    names_prefix = "day_"
  )

summary(Ajalapensis.captures.day)
str(Ajalapensis.captures.day)

env.data$plot <- as.character(env.data$plot)
env.data$trap <- as.character(env.data$trap)
env.data$fieldtrip <- as.character(env.data$fieldtrip)

Ajalapensis.captures.day.env <- left_join(
  env.data,
  Ajalapensis.captures.day,
  by = c("fieldtrip", "plot", "trap")
)
str(env.data)
str(Ajalapensis.captures.day.env)
summary(Ajalapensis.captures.day.env)

## brms------------------------------------------------------------------------
# Create a data frame with total captures per trap per trip.
capts.env <- data.frame(
  Ajalapensis.captures.day.env[, -c(17:31)],
  capts = rowSums(Ajalapensis.captures.day.env[, c(17:31)])
)

# Scale the environmental predictors.
capts.env[, c(4:15)] <- scale(capts.env[, c(4:15)])

# -- Full Model and Variable Selection for Abundance --
mfull.capts <- brm(
  capts ~
    plot +
    TSLF +
    t.med +
    t.max +
    t.max.abs +
    rh.min.abs +
    rh.max.abs +
    VARI.all +
    zentropy +
    tree.density +
    Ajalapensis_perf +
    Ajalapensis_ha90,
  data = capts.env,
  cores = 4,
  # Use Poisson distribution for count data.
  family = "poisson"
)

summary(mfull.capts)
plot(mfull.capts)
bayes_R2(mfull.capts)

# Variable selection for the abundance model.
# First, perform variable selection without cross-validation
vs_capts <- cv_varsel(mfull.capts, validate_search = F)
summary(vs_capts)
summary(vs_capts, deltas = T)

suggest_size(vs_capts)

plot(vs_capts, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.capts)

# Then, perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_capts <- run_cvfun(
  refm_obj,
  K = 10
)

# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
doParallel::registerDoParallel(5)

# Final cv_varsel() run:
cvvs_capts <- cv_varsel(
  refm_obj,
  cv_method = "kfold",
  cvfits = cv_fits_capts,
  nterms_max = 3,
  parallel = TRUE,
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

summary(cvvs_capts, deltas = T, type = c("mean", "lower", "upper"))


plot(cv_proportions(ranking(cvvs_capts)))

plot(cv_proportions(ranking(cvvs_capts), cumulate = T))

plot(cvvs_capts, deltas = T, stats = "mlpd") + theme_minimal()

# -- Fit Selected Model for Abundance --
# 'Ajalapensis_ha90' (hours of activity) is key.
msel.capts.re <- brm(
  capts ~ Ajalapensis_ha90 + (1 | fieldtrip / plot / trap),
  data = capts.env,
  cores = 4,
  family = "poisson",
  iter = 5000,
  control = list(adapt_delta = 0.999)
)

summary(msel.capts.re)
bayes_R2(msel.capts.re)

# -- Model Diagnostics and Visualization --
plot(msel.capts.re)
pp_check(msel.capts.re)

# Plot the conditional effect of hours of activity on captures.
p1 <- conditional_effects(msel.capts.re)

capts.env$plot <- factor(capts.env$plot, levels = c("A1", "A2", "A3", "A4"))

ggplot(p1$Ajalapensis_ha90, aes(x = Ajalapensis_ha90, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = capts.env,
    aes(x = Ajalapensis_ha90, y = capts, color = plot),
    alpha = 0.5,
    size = 2
  ) +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  labs(x = "Hours of activity", y = "Captures") +
  theme_minimal()

# Compare among plots
msel.capts.plot.re <- brm(
  capts ~ plot + (1 | fieldtrip / trap),
  data = capts.env,
  cores = 4,
  family = "poisson",
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.capts.plot.re)
bayes_R2(msel.capts.plot.re)


# -- Model Diagnostics and Visualization --
plot(msel.capts.plot.re)
pp_check(msel.capts.plot.re)

# Plot the conditional effect of plots
p1 <- conditional_effects(msel.capts.plot.re)

ggplot(
  capts.env,
  aes(x = plot, y = capts, colour = plot)
) +
  ggdist::stat_halfeye(
    aes(fill = plot),
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA,
    alpha = 0.5
  ) +
  gghalves::geom_half_point(
    aes(colour = plot),
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  geom_pointrange(
    data = p1$plot,
    aes(
      x = plot,
      y = estimate__,
      ymin = lower__,
      ymax = upper__,
      colour = plot
    ),
    size = 1.2,
    linewidth = 1.2
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Captures") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  theme_minimal()

# Pairwise contrasts
emmeans::emmeans(msel.capts.plot.re, pairwise ~ plot)


## INLA------------------------------------------------------------------------

# -- Data Preparation for INLA --
# N-mixture models require a matrix of counts (sites x occasions) and covariates.

# Create the count matrix.
(y.mat <- as.matrix(Ajalapensis.captures.day.env[, c(17:31)]))

# Create the `inla.mdata` object, which structures the data for the nmix family.
counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  Ajalapensis.captures.day.env$fieldtrip,
  paste(
    Ajalapensis.captures.day.env$fieldtrip,
    Ajalapensis.captures.day.env$plot,
    sep = "_"
  ),
  paste(
    Ajalapensis.captures.day.env$fieldtrip,
    Ajalapensis.captures.day.env$plot,
    Ajalapensis.captures.day.env$trap,
    sep = "_"
  ),
  scale(Ajalapensis.captures.day.env$Ajalapensis_ha90)
)

print(counts.and.count.covs)

# -- Fit N-Mixture Models --
# We fit a series of models to test different hypotheses.

# Model 0: Intercept-only model (null model).
out.inla.env0 <- inla(
  counts.and.count.covs[, c(1:16)] ~ 1,
  data = list(counts.and.count.covs = counts.and.count.covs[, c(1:16)]),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(
        prior = "flat",
        param = numeric()
      )
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.env0, digits = 3)


plot(
  out.inla.env0,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = T
)


# Model 1: Only random effects
out.inla.env1 <- inla(
  counts.and.count.covs[, c(1:16)] ~
    1 +
    f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs[, c(1:16)],
    fieldtrip = counts.and.count.covs$X2,
    plot_id = counts.and.count.covs$X3
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.1,
    prec.intercept = 0.1
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, .1)),
      theta2 = list(param = c(0, .1)),
      theta3 = list(prior = "gamma", param = c(2, 0.5)),
      theta4 = list(prior = "pc.prec", param = c(1, 0.1)),
      theta5 = list(param = c(0, .1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.env1, digits = 3)

# Plot capture and abundance parameters
plot(
  out.inla.env1,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

# Model 2: Hours of activity effect

counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Ajalapensis.captures.day.env$Ajalapensis_ha90),
  paste(
    Ajalapensis.captures.day.env$fieldtrip,
    Ajalapensis.captures.day.env$plot,
    sep = "_"
  )
)

print(counts.and.count.covs)

out.inla.env2 <- inla(
  counts.and.count.covs ~ 1 + ha_90 + f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    plot_id = counts.and.count.covs$X3,
    ha_90 = counts.and.count.covs$X2
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.env2, digits = 3)

# Plot capture and abundance parameters and diagnostics
plot(
  out.inla.env2,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

# Model 3: Plots effects
dummy_vars <- model.matrix(
  ~ factor(
    Ajalapensis.captures.day.env$plot,
    levels = c("A1", "A2", "A3", "A4")
  )
)

counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  dummy_vars[, 2],
  dummy_vars[, 3],
  dummy_vars[, 4]
)

print(counts.and.count.covs)


out.inla.env3 <- inla(
  counts.and.count.covs ~ 1 + plot + f(fieldtrip, model = "iid"),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    plot = factor(
      Ajalapensis.captures.day.env$plot,
      levels = c("A1", "A2", "A3", "A4")
    ), # 'plot' factor for the detection model
    fieldtrip = Ajalapensis.captures.day.env$fieldtrip # 'fieldtrip' for the detection model
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.env3, digits = 3)

# Plot capture and abundance parameters and diagnostics
plot(
  out.inla.env3,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

counts.and.count.covs <- inla.mdata(
  y.mat,
  1
)

print(counts.and.count.covs)

out.inla.env4 <- inla(
  counts.and.count.covs ~ 1 + f(fieldtrip, model = "iid"),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    fieldtrip = Ajalapensis.captures.day.env$fieldtrip # 'fieldtrip' for the detection model
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.env4, digits = 3)

# Plot capture and abundance parameters and diagnostics
plot(
  out.inla.env4,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.env4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

# -- Model Comparison using WAIC --
# Compare models to find the best explanation for abundance patterns.
out.inla.env0$waic$waic
out.inla.env1$waic$waic
out.inla.env2$waic$waic
out.inla.env3$waic$waic
out.inla.env4$waic$waic

inla_env_models <- list(
  out.inla.env0,
  out.inla.env1,
  out.inla.env2,
  out.inla.env3,
  out.inla.env4
)

waic_values <- sapply(inla_env_models, function(model) model$waic$waic)
peff_values <- sapply(inla_env_models, function(model) model$waic$p.eff)
mlik_values <- sapply(inla_env_models, function(model) model$mlik[1])


delta_waic <- waic_values - min(waic_values)

w_akaike <- exp(-0.5 * delta_waic) / sum(exp(-0.5 * delta_waic))

waic_table <- data.frame(
  Model = paste0("Model ", 0:4),
  mlik = mlik_values,
  WAIC = waic_values,
  deltaWAIC = delta_waic,
  pWAIC = peff_values,
  Weight = round(w_akaike, 2)
)

waic_table$formula[1] <- "(lambda ~ 1), (p ~ 1)"
waic_table$formula[2] <- "(lambda ~ 1), (p ~ (1|fieldtrip:plot))"
waic_table$formula[3] <- "(lambda ~ ha90), (p ~ ha90 + (1|fieldtrip:plot))"
waic_table$formula[4] <- "(lambda ~ plot), (p ~ plot + (1|fieldtrip))"
waic_table$formula[5] <- "(lambda ~ 1), (p ~ 1 + (1|fieldtrip))"

# Sort by WAIC (best model first)
waic_table <- waic_table[order(waic_table$WAIC), ]

# Print table
print(waic_table)

# Best models' summaries
summary(out.inla.env2)
summary(out.inla.env3)

summary(plogis(as.matrix(out.inla.env2$summary.linear.predictor[, c(
  "mean",
  "sd",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "mode"
)])))
summary(out.inla.env2$summary.fitted.values)


# Hours of activity effects

ha_xx <- seq(
  min(scale(Ajalapensis.captures.day.env$Ajalapensis_ha90)),
  max(scale(Ajalapensis.captures.day.env$Ajalapensis_ha90)),
  by = 0.1
)

# Draw 1000 samples from the model's posterior
posterior_samples <- inla.posterior.sample(10000, out.inla.env2)

# The 'contents' table describes the structure of the model's output
contents <- out.inla.env2$misc$configs$contents

# Find the numeric index for the intercept
intercept_idx <- contents$start[contents$tag == "(Intercept)"]

# Find the numeric index for the slope
slope_idx <- contents$start[contents$tag == "ha_90"]

# Run the sapply loop
predicted_lines <- sapply(posterior_samples, function(s) {
  # Access the column by its position [ ,1]
  intercept <- s$latent[intercept_idx, 1]
  slope <- s$latent[slope_idx, 1]
  linear_predictor <- intercept + slope * ha_xx
  return(plogis(linear_predictor))
})

rel_p_ha_90 <- data.frame(
  ha_90 = ha_xx,
  median = apply(predicted_lines, 1, quantile, 0.5),
  lower = apply(predicted_lines, 1, quantile, 0.025),
  upper = apply(predicted_lines, 1, quantile, 0.975)
)

pred_df <- data.frame(
  ha_90 = scale(Ajalapensis.captures.day.env$Ajalapensis_ha90),
  trap = Ajalapensis.captures.day.env$trap,
  plot = factor(
    Ajalapensis.captures.day.env$plot,
    levels = c("A1", "A2", "A3", "A4")
  ),
  p = plogis(out.inla.env2$summary.linear.predictor$`0.5quant`),
  ncaps = rowSums(y.mat)
)

ggplot() +
  geom_point(
    data = pred_df,
    aes(x = ha_90, y = p, colour = plot),
    size = 3,
    alpha = 0.3
  ) +
  scale_color_manual(values = turbo(4)) +
  geom_line(
    data = rel_p_ha_90,
    aes(x = ha_90, y = median),
    color = "blue",
    linewidth = 2
  ) +
  geom_ribbon(
    data = rel_p_ha_90,
    aes(x = ha_90, ymin = lower, ymax = upper),
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  labs(x = "Hours of activity", y = "Capture") +
  theme_minimal()

ggplot(pred_df, aes(x = ncaps, y = p)) +
  geom_point(size = 3, alpha = 0.3)

ggplot(pred_df, aes(x = ha_90, y = ncaps)) +
  geom_point(size = 3, alpha = 0.3)

# Plots effects

# Capture

# Get posterior samples for the fixed effects ---
# We'll use the robust index.
post_samples_fixed <- inla.posterior.sample(10000, out.inla.env3)
contents <- out.inla.env3$misc$configs$contents

# Get the numeric index for each fixed effect ---
idx_intercept <- contents$start[contents$tag == "(Intercept)"]
idx_A2 <- contents$start[contents$tag == "plotA2"]
idx_A3 <- contents$start[contents$tag == "plotA3"]
idx_A4 <- contents$start[contents$tag == "plotA4"]

# Extract the posterior samples for each coefficient ---
# This loop gets the values for each sample and stores them in vectors
b_intercept <- sapply(post_samples_fixed, function(s) {
  s$latent[idx_intercept, 1]
})
b_A2 <- sapply(post_samples_fixed, function(s) s$latent[idx_A2, 1])
b_A3 <- sapply(post_samples_fixed, function(s) s$latent[idx_A3, 1])
b_A4 <- sapply(post_samples_fixed, function(s) s$latent[idx_A4, 1])

# Calculate all 6 pairwise contrasts on the logit scale ---
# The coefficients themselves are already differences relative to the intercept.
# To compare other levels, we subtract their respective coefficients.
contrasts_p <- list(
  `A2 - A1` = b_A2,
  `A3 - A1` = b_A3,
  `A4 - A1` = b_A4,
  `A2 - A3` = b_A2 - b_A3,
  `A2 - A4` = b_A2 - b_A4,
  `A3 - A4` = b_A3 - b_A4
)

# Summarize the results into a data frame ---
contrast_summary_p <- t(sapply(
  contrasts_p,
  quantile,
  probs = c(0.025, 0.5, 0.975)
))

contrast_summary_p_df <- as.data.frame(contrast_summary_p)
colnames(contrast_summary_p_df) <- c("Lower.HPD", "Median", "Upper.HPD")

print(contrast_summary_p_df, digits = 3)

# Calculate detection probability for each plot type ---
# The model is: logit(p) = Intercept + effect_A2 + effect_A3 + effect_A4
# We use plogis() to transform back to the probability scale

# Loop through samples to get posteriors for p
p_A1_post <- sapply(post_samples_fixed, function(s) {
  plogis(s$latent[idx_intercept, 1])
})
p_A2_post <- sapply(post_samples_fixed, function(s) {
  plogis(s$latent[idx_intercept, 1] + s$latent[idx_A2, 1])
})
p_A3_post <- sapply(post_samples_fixed, function(s) {
  plogis(s$latent[idx_intercept, 1] + s$latent[idx_A3, 1])
})
p_A4_post <- sapply(post_samples_fixed, function(s) {
  plogis(s$latent[idx_intercept, 1] + s$latent[idx_A4, 1])
})

# Summarize the posteriors for plotting ---

p_post <- data.frame(
  A1 = p_A1_post,
  A2 = p_A2_post,
  A3 = p_A3_post,
  A4 = p_A4_post
)


p_post_long <- p_post %>%
  pivot_longer(
    cols = everything(),
    names_to = "Plot",
    values_to = "p"
  ) %>%
  mutate(Plot = factor(Plot, levels = c("A1", "A2", "A3", "A4")))

ggplot(p_post_long, aes(x = Plot, y = p, fill = Plot, color = Plot)) +
  ggdist::stat_halfeye(
    adjust = 1,
    width = .8,
    justification = -.1,
    alpha = 0.5,
    # Sets the credible interval to 95%
    .width = 0.95,
    # Display the median and Highest Density Interval
    point_interval = "median_hdi"
  ) +

  # Wrap the labels argument in a function
  scale_y_log10() +

  # Add other labels and scales
  labs(
    y = "Capture probability",
    x = "Fire regime"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  # Hide the redundant legend
  guides(fill = "none")

# Abundance

# Get samples from the posterior hyperparameters ---
# Increase sample size for more stable estimates
post_samples <- inla.hyperpar.sample(10000, out.inla.env3)

# Extract the posterior samples for each beta coefficient for abundance.
# Use the exact names from our model's summary output.
beta1 <- post_samples[, "beta[1] for NMixNB observations"] # Intercept (A1)
beta2 <- post_samples[, "beta[2] for NMixNB observations"] # A2 vs A1
beta3 <- post_samples[, "beta[3] for NMixNB observations"] # A3 vs A1
beta4 <- post_samples[, "beta[4] for NMixNB observations"] # A4 vs A1

# Calculate all 6 pairwise contrasts ---
contrasts <- list(
  `A2 - A1` = beta2,
  `A3 - A1` = beta3,
  `A4 - A1` = beta4,
  `A2 - A3` = beta2 - beta3,
  `A2 - A4` = beta2 - beta4,
  `A3 - A4` = beta3 - beta4
)

# Summarize the results into a data frame
contrast_summary <- t(sapply(contrasts, quantile, probs = c(0.025, 0.5, 0.975)))

# Make it a clean data frame and print it
contrast_summary_df <- as.data.frame(contrast_summary)
colnames(contrast_summary_df) <- c("Lower.HPD", "Median", "Upper.HPD")

print(contrast_summary_df, digits = 3)

# Calculate lambda for each plot type using the posterior samples ---
# These are on the log scale, so we use exp() to transform them back
lambda_post <- data.frame(
  A1 = exp(beta1),
  A2 = exp(beta1 + beta2),
  A3 = exp(beta1 + beta3),
  A4 = exp(beta1 + beta4)
)


lambda_post_long <- lambda_post %>%
  pivot_longer(
    cols = everything(),
    names_to = "Plot",
    values_to = "lambda"
  ) %>%
  mutate(Plot = factor(Plot, levels = c("A1", "A2", "A3", "A4")))

y_limits <- quantile(lambda_post_long$lambda, c(0.025, 0.975))
y_min <- y_limits[1]
y_max <- y_limits[2]

ggplot(lambda_post_long, aes(x = Plot, y = lambda, fill = Plot, color = Plot)) +

  ggdist::stat_halfeye(
    adjust = 1,
    width = .8,
    justification = -.1,
    alpha = 0.5,
    # Sets the credible interval to 95%
    .width = 0.95,
    # Display the median and Highest Density Interval
    point_interval = "median_hdi"
  ) +

  scale_y_log10(
    labels = function(x) scales::number(x, accuracy = 1)
  ) +

  # Add other labels and scales
  labs(
    y = "Abundance",
    x = "Plot Type"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_color_viridis_d(option = "turbo") +
  coord_cartesian(ylim = c(y_min, y_max)) +
  theme_minimal() +
  guides(fill = "none")


# Sex with environment -----------------------------------------------
head(svl.env)

# Set immature individuals as NA
svl.env$sex[svl.env$sex == "I"] <- NA

# Count captures by sex
table(svl.env$sex)

# Fit full model
mfull.sex <- brm(
  sex ~
    plot +
    TSLF +
    t.med +
    t.max +
    t.max.abs +
    rh.min.abs +
    rh.max.abs +
    VARI.all +
    zentropy +
    tree.density +
    Ajalapensis_perf +
    Ajalapensis_ha90,
  data = svl.env,
  cores = 4,
  family = "bernoulli"
)

summary(mfull.sex)
plot(mfull.sex)

# Perform variable selection using projpred
# First, perform variable selection without cross-validation
vs_sex <- cv_varsel(mfull.sex, validate_search = F)
summary(vs_sex)
summary(vs_sex, deltas = T)

suggest_size(vs_sex)

plot(vs_sex, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.sex)

# Then, perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_sex <- run_cvfun(
  refm_obj,

  K = 10
  ###
)
# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
doParallel::registerDoParallel(5)
# Final cv_varsel() run:
cvvs_sex <- cv_varsel(
  refm_obj,
  cv_method = "kfold",
  cvfits = cv_fits_sex,
  nterms_max = 4,
  parallel = TRUE,
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

summary(cvvs_sex, deltas = T, type = c("mean", "lower", "upper"))


plot(cv_proportions(ranking(cvvs_sex)))

plot(cv_proportions(ranking(cvvs_sex), cumulate = T))

plot(cvvs_sex, deltas = T, stats = "mlpd") + theme_minimal()

# Fit selected model
# 'Ajalapensis_perf' is key
msel.sex.re <- brm(
  sex ~ Ajalapensis_perf + (1 | fieldtrip / plot / trap),
  data = svl.env,
  cores = 4,
  family = "bernoulli",
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.sex.re)
bayes_R2(msel.sex.re)

# Fit null model
msel.null.re <- brm(
  sex ~ 1 + (1 | fieldtrip / plot / trap),
  data = svl.env,
  cores = 4,
  family = "bernoulli",
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.null.re)
bayes_R2(msel.null.re)

# Fit plot effects model
msel.sex.plot.re <- brm(
  sex ~ plot + (1 | fieldtrip / trap),
  data = svl.env,
  cores = 4,
  family = "bernoulli",
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

summary(msel.sex.plot.re)
plot(msel.sex.plot.re)

pp_check(msel.sex.plot.re)
bayes_R2(msel.sex.plot.re)

# Compare models
loo(msel.sex.re, msel.null.re, msel.sex.plot.re)

# Ajalapensis_perf is not significant
conditional_effects(msel.sex.re)
par(mfrow = c(1, 1))
cdplot(as.factor(sex) ~ Ajalapensis_perf, data = svl.env)

svl.env$sex <- as.integer(as.factor(svl.env$sex)) - 1

# Conditional effects plot to visualize differences among plots
p1 <- conditional_effects(msel.sex.plot.re)

ggplot(
  svl.env,
  aes(x = plot, y = sex, colour = plot)
) +

  ggdist::geom_dotsinterval(aes(fill = plot)) +
  geom_pointrange(
    data = p1$plot,
    aes(
      x = plot,
      y = estimate__,
      ymin = lower__,
      ymax = upper__,
      colour = plot
    ),
    size = 1.2,
    linewidth = 1.2
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Proportion of males") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity")

# Pairwise comparisons
emmeans::emmeans(msel.sex.plot.re, pairwise ~ plot)

# Count sex per plot
table(svl.env$sex, svl.env$plot)


# Body condition with environment -----------------------------------------------

# Remove cases where mass = NA
completos <- complete.cases(svl.env[, c("mass")])
condition.data <- droplevels(svl.env[completos, ])

# View data (svl, mass)
plot(condition.data$svl, condition.data$mass, las = 1, bty = "n")
plot(log(condition.data$svl), log(condition.data$mass), las = 1, bty = "n")

# -- Calculate Scaled Mass Index (SMI) --
# SMI is a measure of body condition that corrects for body size.

# Calculate body condition as "scaled mass index" (Peig & Green, 2009)
# install.packages("lmodel2", dependencies=T)
library(lmodel2)

# Perform a Standardized Major Axis (SMA) regression of log(mass) on log(svl).
sma.jalapensis <- lmodel2(
  log(mass + 1) ~ log(svl),
  data = condition.data,
  "relative",
  "relative",
  99
)
sma.jalapensis
par(mfrow = c(1, 2))

plot(sma.jalapensis, "OLS")
plot(sma.jalapensis, "SMA")
plot(sma.jalapensis, "MA")
plot(sma.jalapensis, "RMA")

par(mfrow = c(1, 1))

attach(condition.data)

# Calculate SMI for each individual.
smi <- mass * ((mean(svl) / svl)^2.027180) # 2.027180 = slope of SMA regression
boxplot(smi)
plot(svl, smi, las = 1, bty = "n")
detach(condition.data)
condition.data$smi <- smi
boxplot(smi ~ sex, data = condition.data)

rm(smi)
detach("package:lmodel2", unload = TRUE)

# -- Model Body Condition --
# Test the effect of environmental variables on SMI.

# Fit full model
mfull.smi <- brm(
  smi ~
    plot +
    TSLF +
    t.med +
    t.max +
    t.max.abs +
    rh.min.abs +
    rh.max.abs +
    VARI.all +
    zentropy +
    tree.density +
    Ajalapensis_perf +
    Ajalapensis_ha90,
  data = condition.data,
  cores = 4
)

summary(mfull.smi)
plot(mfull.smi)
bayes_R2(mfull.smi)

# Perform variable selection without cross-validation
vs_smi <- cv_varsel(mfull.smi, validate_search = F)
summary(vs_smi)
summary(vs_smi, deltas = T)

suggest_size(vs_smi)

plot(vs_smi, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.smi)

# Perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_smi <- run_cvfun(
  refm_obj,

  K = 10
  ###
)
# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
doParallel::registerDoParallel(5)
# Final cv_varsel() run:
cvvs_smi <- cv_varsel(
  refm_obj,
  cv_method = "kfold",
  cvfits = cv_fits_smi,
  nterms_max = 3,
  parallel = TRUE,
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

summary(cvvs_smi, deltas = T, type = c("mean", "lower", "upper"))


plot(cv_proportions(ranking(cvvs_smi)))

plot(cv_proportions(ranking(cvvs_smi), cumulate = T))

plot(cvvs_smi, deltas = T, stats = "mlpd") + theme_minimal()

# Fit selected model.
# 'Ajalapensis_perf' is key
msel.smi.re <- brm(
  smi ~ Ajalapensis_perf + (1 | fieldtrip / plot / trap),
  data = condition.data,
  cores = 4,
  iter = 5000,
  control = list(adapt_delta = 0.999),
  silent = 0
)

summary(msel.smi.re)
plot(msel.smi.re)

pp_check(msel.smi.re)
bayes_R2(msel.smi.re)

# Plot the conditional effect of performance on body condition.
p1 <- conditional_effects(msel.smi.re)

ggplot(p1$Ajalapensis_perf, aes(x = Ajalapensis_perf, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = condition.data,
    aes(x = Ajalapensis_perf, y = smi, color = plot),
    alpha = 0.5,
    size = 2
  ) +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  labs(x = "Locomotor performance", y = "Scaled mass index") +
  theme_minimal()

# Compare among plots
msel.smi.plot.re <- brm(
  smi ~ plot + (1 | fieldtrip / trap),
  data = condition.data,
  cores = 4,
  iter = 5000,
  control = list(adapt_delta = 0.999)
)

summary(msel.smi.plot.re)
plot(msel.smi.plot.re)

pp_check(msel.smi.plot.re)
bayes_R2(msel.smi.plot.re)

# Conditional effects plot
p1 <- conditional_effects(msel.smi.plot.re)

ggplot(
  condition.data,
  aes(x = plot, y = smi, colour = plot)
) +
  ggdist::stat_halfeye(
    aes(fill = plot),
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA,
    alpha = 0.5
  ) +
  gghalves::geom_half_point(
    aes(colour = plot),
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  geom_pointrange(
    data = p1$plot,
    aes(
      x = plot,
      y = estimate__,
      ymin = lower__,
      ymax = upper__,
      colour = plot
    ),
    size = 1.2,
    linewidth = 1.2
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Scaled mass index") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  theme_minimal()

# Pairwise comparisons
emmeans::emmeans(msel.smi.plot.re, pairwise ~ plot)


# Fire regime components (fixed and rotating plots) ---------------------------

# -- Data Preparation --
# This section merges the fixed plot data with data from rotating plots
# (from a different source) and links them to fire regime variables.

# Load data from rotating plots.
brunadata <- readxl::read_excel(
  "Ameivula_jalapensis_EESGT_BrunaGomes.xlsx",
  na = "NA"
)
glimpse(brunadata)

brunadata <- dplyr::rename(brunadata, "plot" = "sampling_point")
brunadata$trap <- paste0(brunadata$fieldtrip, brunadata$Y)
glimpse(brunadata)
summary(brunadata)
brunadata$plot <- as.character(brunadata$plot)

# Load fire regime data for all traps.
fire.regimes.traps <- read.csv("fire_regimes_traps_df.csv")

# Merge the datasets.
brunadata <- left_join(
  brunadata,
  fire.regimes.traps[, -1],
  by = c("plot", "trap", "treatment")
)
summary(brunadata)

## SVL -------------------------------------------------------------------

# Simple visualization
ggplot(brunadata, aes(y = svl_mm, x = jitter(severity))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(brunadata, aes(y = svl_mm, x = jitter(freq))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(brunadata, aes(y = svl_mm, x = MeanTSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(brunadata, aes(y = svl_mm, x = jitter(TUF))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(brunadata, aes(y = svl_mm, x = treatment)) +
  geom_boxplot()

ggplot(brunadata, aes(y = svl_mm, x = treat_code)) +
  geom_boxplot()

ggplot(brunadata, aes(y = svl_mm, x = as.factor(TUF))) +
  geom_boxplot()

ggplot(brunadata, aes(y = svl_mm, x = date)) +
  geom_point(alpha = 0.5)

ggplot(brunadata, aes(x = svl_mm, y = weight_animal)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(brunadata, aes(x = log(svl_mm), y = log(weight_animal))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")

plot(svl_mm ~ freq, data = brunadata)
plot(svl_mm ~ MeanTSLF, data = brunadata)
plot(svl_mm ~ TUF, data = brunadata)
boxplot(svl_mm ~ treatment, data = brunadata)

# Filter variables
svlbruna <- brunadata[, c(
  "fieldtrip",
  "date",
  "plot",
  "trap_code",
  "TUF",
  "treatment",
  "treat_code",
  "recapture",
  "svl_mm",
  "weight_animal",
  "cicatriz_umbilical",
  "sex",
  "ovada",
  "trap",
  "freq",
  "MeanTSLF",
  "MFRI",
  "severity",
  "TSLF"
)]

svlbruna <- filter(svlbruna, svl_mm < 70)

ggplot(svlbruna, aes(y = svl_mm, x = jitter(severity))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(y = svl_mm, x = jitter(MFRI))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(y = svl_mm, x = jitter(freq))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(y = svl_mm, x = MeanTSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(y = svl_mm, x = jitter(TUF))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(y = svl_mm, x = treatment)) +
  geom_boxplot()

ggplot(svlbruna, aes(y = svl_mm, x = treat_code)) +
  geom_boxplot()

ggplot(svlbruna, aes(y = svl_mm, x = as.factor(TUF))) +
  geom_boxplot()

ggplot(svlbruna, aes(y = svl_mm, x = date)) +
  geom_point(alpha = 0.5)

ggplot(svlbruna, aes(x = svl_mm, y = weight_animal)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svlbruna, aes(x = log(svl_mm), y = log(weight_animal))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")

names(svlbruna) <- c(
  "fieldtrip",
  "date",
  "plot",
  "trap_code",
  "TUF",
  "treatment",
  "treat_code",
  "recapture",
  "svl",
  "mass",
  "umb_scar",
  "sex",
  "eggs",
  "trap",
  "freq",
  "MeanTSLF",
  "MFRI",
  "severity",
  "TSLF"
)

# Number of recaptures
table(svlbruna$recapture)

# Number of captures (individuals)
nrow(svlbruna) - table(svlbruna$recapture)

# Mean recapture rate
table(svlbruna$recapture) / nrow(svlbruna)

# Sex
table(svlbruna$sex)

# Filter data
svl.data <- data[, c(
  "Fieldtrip",
  "Date",
  "Month",
  "Year",
  "Plot",
  "Trap",
  "VoucherNumber",
  "Recapture",
  "SVL",
  "Mass",
  "Sex",
  "Eggs"
)]

# Rename data
names(svl.data) <- c(
  "fieldtrip",
  "date",
  "month",
  "year",
  "plot",
  "trap",
  "ID",
  "recapture",
  "svl",
  "mass",
  "sex",
  "eggs"
)

svl.data$trap <- as.factor(paste(svl.data$plot, svl.data$trap, sep = "_"))
svl.data$fieldtrip <- as.factor(svl.data$fieldtrip)

# Merge datasets
svl.fire <- left_join(
  svl.data,
  fire.regimes.traps[, -c(1:2, 10)],
  by = c("plot", "trap")
)

svl.fire <- left_join(
  svl.fire,
  env.data[, c(1:3, 5)],
  by = c("plot", "trap", "fieldtrip")
)

svl.fire$date <- as.POSIXct(svl.fire$date, format = "%d/%m/%Y", tz = "UTC")
summary(svl.fire)
unique(svl.fire$severity)
plot(svl ~ severity, data = svl.fire)

svl.merged <- dplyr::full_join(svl.fire, svlbruna)
svl.merged$month <- month(svl.merged$date)
svl.merged$year <- year(svl.merged$date)

summary(svl.merged)
table(svl.merged$TUF == svl.merged$TSLF)

# Transform TSLF to months and set to zero, when recently burned ("Q")
svl.merged$TSLF[!is.na(svl.merged$TUF)] <- svl.merged$TSLF[
  !is.na(svl.merged$TUF)
] *
  12
svl.merged$TSLF[svl.merged$treatment == "Q"] <- 0

# Simple plots
ggplot(svl.merged, aes(y = svl, x = jitter(severity))) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = jitter(freq))) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = MeanTSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = TSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = jitter(TUF))) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = treat_code)) +
  geom_boxplot()

ggplot(svl.merged, aes(y = svl, x = as.factor(TUF))) +
  geom_boxplot()

ggplot(svl.merged, aes(y = svl, x = as.Date(date))) +
  labs(x = "Date", y = "Snout-vent length (mm)") +
  scale_x_date(date_minor_breaks = "1 month") +
  theme_minimal() +
  geom_hline(aes(yintercept = 40), linetype = 2, alpha = 0.5) +
  geom_point(alpha = 0.5)

svl.merged %>%
  filter(!is.na(month)) %>%
  # In aes(), use month.abb[month] to get the abbreviation.
  # Set the levels to month.abb to ensure correct chronological order.
  ggplot(aes(y = svl, x = factor(month.abb[month], levels = month.abb))) +
  labs(x = "", y = "Snout-vent length (mm)") +
  theme_minimal() +
  geom_hline(aes(yintercept = 40), linetype = 2, alpha = 0.5) +
  geom_boxplot()

table(svl.merged$month)

ggplot(svl.merged, aes(x = svl, y = mass)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(svl.merged, aes(x = svl, y = mass)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")

ggplot(svl.merged, aes(x = log(svl), y = log(mass))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")

# Fire regime variables pairwise plots
ggpairs(svl.merged[, c(13:16, 18)])

# Remove collinear variables
usdm::vifstep(svl.merged[, c(13:16, 18)], keep = "severity", th = 2)
usdm::vif(svl.merged[, c(13:16, 18)])

# Scale the fire regime predictors.
svl.merged[, c(13:16, 18)] <- scale(svl.merged[, c(13:16, 18)])

# Testing the effects of fire components
# Model with Mean Time Since Last Fire (MeanTSLF).
msel.svl.MeanTSLF.re <- brm(
  svl ~ MeanTSLF + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  cores = 4
)

summary(msel.svl.MeanTSLF.re)
bayes_R2(msel.svl.MeanTSLF.re)

# Diagnostic plots
plot(msel.svl.MeanTSLF.re)
pp_check(msel.svl.MeanTSLF.re)

# Model with Time Since Last Fire (TSLF).
msel.svl.TSLF.re <- brm(
  svl ~ TSLF + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  cores = 4
)

summary(msel.svl.TSLF.re)
bayes_R2(msel.svl.TSLF.re)

# Diagnostic plots
plot(msel.svl.TSLF.re)
pp_check(msel.svl.TSLF.re)

# Model with fire severity.
msel.svl.severity.re <- brm(
  svl ~ severity + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  cores = 4
)

summary(msel.svl.severity.re)
bayes_R2(msel.svl.severity.re)

# Diagnostic plots
plot(msel.svl.severity.re)
pp_check(msel.svl.severity.re)

# Null model for comparison.
msel.svl.firenull.re <- brm(
  svl ~ 1 + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  cores = 4
)

summary(msel.svl.firenull.re)
bayes_R2(msel.svl.firenull.re)

# Diagnostic plots
plot(msel.svl.firenull.re)
pp_check(msel.svl.firenull.re)

# Compare models using Leave-One-Out Cross-Validation (LOO).
loo(
  msel.svl.MeanTSLF.re,
  msel.svl.TSLF.re,
  msel.svl.severity.re,
  msel.svl.firenull.re
)


# Conditional effects plot
p1 <- conditional_effects(msel.svl.MeanTSLF.re)

# Plot the conditional effect of the best predictor (e.g., MeanTSLF).
ggplot(p1$MeanTSLF, aes(x = MeanTSLF, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = 40), linetype = 2) +
  geom_point(
    data = svl.merged,
    aes(x = MeanTSLF, y = svl),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Mean fire interval", y = "Snout-vent length (mm)") +
  theme_minimal()

# Conditional effects plot
p1 <- conditional_effects(msel.svl.severity.re)

ggplot(p1$severity, aes(x = severity, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = 40), linetype = 2) +
  geom_point(
    data = svl.merged,
    aes(x = severity, y = svl),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Mean fire severity", y = "Snout-vent length (mm)") +
  theme_minimal()


## Body condition -----------------------------------------------

# Remove cases where mass = NA
complete <- complete.cases(svl.merged[, c("mass")])
condition.data <- droplevels(svl.merged[complete, ])

# View data (svl, mass)
plot(condition.data$svl, condition.data$mass, las = 1, bty = "n")
plot(log(condition.data$svl), log(condition.data$mass), las = 1, bty = "n")


mod1 <- lm(condition.data$mass ~ condition.data$svl)
summary(mod1)

mod2 <- lm(sqrt(condition.data$mass) ~ condition.data$svl)
summary(mod2)

mod3 <- lm((condition.data$mass^(1 / 3)) ~ condition.data$svl)
summary(mod3)

mod4 <- lm(log(condition.data$mass^(1 / 3)) ~ log(condition.data$svl))
summary(mod4)

AIC(mod1, mod2, mod3, mod4)

plot(
  condition.data$svl,
  log(condition.data$mass^(1 / 3)),
  pch = 21,
  cex = 1,
  col = rgb(217, 72, 1, 100, maxColorValue = 255),
  bg = rgb(253, 174, 107, 100, maxColorValue = 255),
  las = 1,
  bty = "n",
  ylab = "log(Mass ^ 1/3 (g))",
  xlab = "Comprimento Rostro-Cloacal (mm)"
) +
  abline(
    a = lm(log(condition.data$mass^(1 / 3)) ~ condition.data$svl),
    col = "darkred"
  )

zcrit <- qnorm(0.9999)
summary(rstandard(mod4))
st_residuals <- condition.data %>%
  filter(abs(rstandard(mod4)) > zcrit)
st_residuals

# Remove outliers
condition.data <- condition.data %>%
  filter(abs(rstandard(mod4)) < zcrit)

mod4 <- lm(log(condition.data$mass^(1 / 3)) ~ log(condition.data$svl))
summary(mod4)

plot(
  condition.data$svl,
  log(condition.data$mass^(1 / 3)),
  pch = 21,
  cex = 1,
  col = rgb(217, 72, 1, 100, maxColorValue = 255),
  bg = rgb(253, 174, 107, 100, maxColorValue = 255),
  las = 1,
  bty = "n",
  ylab = "log(Mass ^ 1/3 (g))",
  xlab = "Comprimento Rostro-Cloacal (mm)"
) +
  abline(
    a = lm(log(condition.data$mass^(1 / 3)) ~ condition.data$svl),
    col = "darkred"
  )

st_residuals <- condition.data %>%
  filter(abs(rstandard(mod4)) > zcrit)
st_residuals

p <- length(coefficients(mod4))
n <- length(fitted(mod4))
ratio <- p / n
leverage <- condition.data %>%
  filter(hatvalues(mod4) > (3 * ratio))
leverage

cutoff <- 4 / (dim(condition.data)[1])
cooks_dis <- condition.data %>%
  filter(cooks.distance(mod4) > cutoff)
cooks_dis

bonf <- car::outlierTest(mod4, n.max = 9999)

condition.data[names(bonf$rstudent), ]

# Calculate body condition as "scaled mass index" (Peig & Green, 2009)
# install.packages("lmodel2", dependencies=T)
library(lmodel2)

sma.jalapensis <- lmodel2(
  log(mass + 1) ~ log(svl),
  data = condition.data,
  "relative",
  "relative",
  99
)
sma.jalapensis
par(mfrow = c(1, 2))

plot(sma.jalapensis, "OLS")
plot(sma.jalapensis, "SMA")
plot(sma.jalapensis, "MA")
plot(sma.jalapensis, "RMA")

par(mfrow = c(1, 1))

attach(condition.data)
smi <- mass * ((mean(svl) / svl)^1.929048) # 1.929048= slope of SMA regression
boxplot(smi)
plot(svl, smi, las = 1, bty = "n")

detach(condition.data)
condition.data$smi <- smi
rm(smi)
detach("package:lmodel2", unload = TRUE)

# Testing the effects of fire components
# Mean time since last fire
msel.smi.MeanTSLF.re <- brm(
  smi ~ MeanTSLF + (1 | fieldtrip / plot / trap),
  data = condition.data,
  cores = 4
)

summary(msel.smi.MeanTSLF.re)
bayes_R2(msel.smi.MeanTSLF.re)

# Diagnostic plots
plot(msel.smi.MeanTSLF.re)
pp_check(msel.smi.MeanTSLF.re)

# Time since last fire
msel.smi.TSLF.re <- brm(
  smi ~ TSLF + (1 | fieldtrip / plot / trap),
  data = condition.data,
  cores = 4
)

summary(msel.smi.TSLF.re)
bayes_R2(msel.smi.TSLF.re)

# Diagnostic plots
plot(msel.smi.TSLF.re)
pp_check(msel.smi.TSLF.re)

# Fire severity
msel.smi.severity.re <- brm(
  smi ~ severity + (1 | fieldtrip / plot / trap),
  data = condition.data,
  cores = 4
)

summary(msel.smi.severity.re)
bayes_R2(msel.smi.severity.re)

# Diagnostic plots
plot(msel.smi.severity.re)
pp_check(msel.smi.severity.re)

msel.smi.firenull.re <- brm(
  smi ~ 1 + (1 | fieldtrip / plot / trap),
  data = condition.data,
  cores = 4
)

summary(msel.smi.firenull.re)
bayes_R2(msel.smi.firenull.re)

# Diagnostic plots
plot(msel.smi.firenull.re)
pp_check(msel.smi.firenull.re)

# Model comparisons
loo(
  msel.smi.MeanTSLF.re,
  msel.smi.TSLF.re,
  msel.smi.severity.re,
  msel.smi.firenull.re
)


# Conditional effects plot
p1 <- conditional_effects(msel.smi.TSLF.re)

ggplot(p1$TSLF, aes(x = TSLF, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = condition.data,
    aes(x = TSLF, y = smi),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Time since last fire", y = "Scaled mass index") +
  theme_minimal()

# Conditional effects plot
p1 <- conditional_effects(msel.smi.severity.re)

ggplot(p1$severity, aes(x = severity, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = condition.data,
    aes(x = severity, y = smi),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Mean fire severity", y = "Scaled mass index") +
  theme_minimal()

## Sex ratio -----------------------------------------------
head(svl.merged)
table(svl.merged$sex)

# Standardize data
svl.merged$sex[svl.merged$sex == "I"] <- NA
svl.merged$sex[svl.merged$sex == "j"] <- NA
svl.merged$sex[svl.merged$sex == "J"] <- NA
svl.merged$sex[svl.merged$sex == "f"] <- "F"
svl.merged$sex[svl.merged$sex == "m"] <- "M"

table(svl.merged$sex)

# Simple plots
cdplot(as.factor(sex) ~ TSLF, data = svl.merged)
cdplot(as.factor(sex) ~ MeanTSLF, data = svl.merged)
cdplot(as.factor(sex) ~ severity, data = svl.merged)
cdplot(as.factor(sex) ~ freq, data = svl.merged)

# Testing the effects of fire components
# Mean time since last fire
msel.sex.MeanTSLF.re <- brm(
  sex ~ MeanTSLF + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  family = "bernoulli",
  cores = 4
)

summary(msel.sex.MeanTSLF.re)
bayes_R2(msel.sex.MeanTSLF.re)

# Diagnostic plots
plot(msel.sex.MeanTSLF.re)
pp_check(msel.sex.MeanTSLF.re)

# Time since last fire
msel.sex.TSLF.re <- brm(
  sex ~ TSLF + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  family = "bernoulli",
  cores = 4
)

summary(msel.sex.TSLF.re)
bayes_R2(msel.sex.TSLF.re)

# Diagnostic plots
plot(msel.sex.TSLF.re)
pp_check(msel.sex.TSLF.re)

# Fire severity
msel.sex.severity.re <- brm(
  sex ~ severity + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  family = "bernoulli",
  cores = 4
)

summary(msel.sex.severity.re)
bayes_R2(msel.sex.severity.re)

# Diagnostic plots
plot(msel.sex.severity.re)
pp_check(msel.sex.severity.re)

# Null model
msel.sex.firenull.re <- brm(
  sex ~ 1 + (1 | fieldtrip / plot / trap),
  data = svl.merged,
  family = "bernoulli",
  cores = 4
)

summary(msel.sex.firenull.re)
bayes_R2(msel.sex.firenull.re)

# Diagnostic plots
plot(msel.sex.firenull.re)
pp_check(msel.sex.firenull.re)

# Model comparisons
loo(
  msel.sex.MeanTSLF.re,
  msel.sex.TSLF.re,
  msel.sex.severity.re,
  msel.sex.firenull.re
)


# Conditional effects plot
p1 <- conditional_effects(msel.sex.TSLF.re)

ggplot(p1$TSLF, aes(x = TSLF, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = svl.merged,
    aes(x = TSLF, y = as.integer(sex) - 1),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Time since last fire", y = "Proportion of males") +
  theme_minimal()

## Abundance------------------------------------------------------------------------
# Data preparation
str(svl.merged)


svl.merged <- svl.merged %>%
  group_by(fieldtrip) %>%
  mutate(sampling_day = dense_rank(date)) %>% # Rank the dates within each month
  ungroup()

summary(svl.merged$sampling_day)

trap_locations <- fire.regimes.traps %>%
  distinct(plot, trap, treatment)

trap_locations$fieldtrip <- substr(trap_locations$trap, 1, 2)
trap_locations$fieldtrip[trap_locations$plot == "A2"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A1"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A3"] <- NA
trap_locations$fieldtrip[trap_locations$plot == "A4"] <- NA

trap_locations_heitor <- trap_locations[is.na(trap_locations$fieldtrip), ]
trap_locations_bruna <- trap_locations[!is.na(trap_locations$fieldtrip), ]

trap_locations_heitor <- reshape::expand.grid.df(
  trap_locations_heitor[, -4],
  data.frame(fieldtrip = as.character(c(1:4)))
)

trap_locations <- rbind(trap_locations_bruna, trap_locations_heitor)

# Now, create the full dataset of all possible sampling events
# `crossing` will combine every row from `trap_locations` with every sampling day.
all_sampling_events <- crossing(trap_locations, sampling_day = 1:15)

# Count captures only for existing combinations in your data
capture_counts <- svl.merged %>%
  count(fieldtrip, plot, trap, sampling_day, name = "freq")

# Join the actual captures to the master list of all sampling events
Ajalapensis.captures.day <- all_sampling_events %>%
  left_join(
    capture_counts,
    by = c("fieldtrip", "plot", "trap", "sampling_day")
  ) %>%
  # If a trap on a given day had no captures, its freq is NA. Change that to 0.
  mutate(freq = replace_na(freq, 0))

# View the result
print(Ajalapensis.captures.day)


summary(Ajalapensis.captures.day)
sum(Ajalapensis.captures.day$freq)


Ajalapensis.captures.day <- Ajalapensis.captures.day %>%
  pivot_wider(
    names_from = sampling_day,
    values_from = freq,
    names_prefix = "day_"
  )

summary(Ajalapensis.captures.day)
str(Ajalapensis.captures.day)

Ajalapensis.captures.day.fire <- left_join(
  Ajalapensis.captures.day,
  fire.regimes.traps[, -c(1, 2)],
  by = c("plot", "trap", "treatment")
)

Ajalapensis.captures.day.fire <- left_join(
  Ajalapensis.captures.day.fire,
  env.data[, c(1:3, 5)],
  by = c("plot", "trap", "fieldtrip")
)
summary(Ajalapensis.captures.day.fire)

Ajalapensis.captures.day.fire$TSLF.y[is.na(
  Ajalapensis.captures.day.fire$TSLF.y
)] <- Ajalapensis.captures.day.fire$TSLF.x[is.na(
  Ajalapensis.captures.day.fire$TSLF.y
)] *
  12

Ajalapensis.captures.day.fire$TSLF.y[
  Ajalapensis.captures.day.fire$treatment == "Q"
] <- 0
Ajalapensis.captures.day.fire <- Ajalapensis.captures.day.fire[, -24]
Ajalapensis.captures.day.fire <- dplyr::rename(
  Ajalapensis.captures.day.fire,
  "TSLF" = "TSLF.y"
)

summary(Ajalapensis.captures.day.fire)

## brms------------------------------------------------------------------------
Ajalapensis.captures.fire <- Ajalapensis.captures.day.fire[, c(1:4, 20:24)]

Ajalapensis.captures.fire$capts <- rowSums(Ajalapensis.captures.day.fire[, c(
  5:19
)])
summary(Ajalapensis.captures.fire)
hist(Ajalapensis.captures.fire$capts)

# Simple plots
ggplot(data = Ajalapensis.captures.fire, aes(y = capts, x = TSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(data = Ajalapensis.captures.fire, aes(y = capts, x = MeanTSLF)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(data = Ajalapensis.captures.fire, aes(y = capts, x = severity)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

ggplot(data = Ajalapensis.captures.fire, aes(y = capts, x = freq)) +
  geom_point(alpha = 0.5) +
  geom_smooth()

# Testing the effects of fire components
# Mean time since last fire
msel.capts.MeanTSLF.re <- brm(
  capts ~ MeanTSLF + (1 | fieldtrip / plot / trap),
  data = Ajalapensis.captures.fire,
  family = "poisson",
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(msel.capts.MeanTSLF.re)
bayes_R2(msel.capts.MeanTSLF.re)

# Diagnostic plots
plot(msel.capts.MeanTSLF.re)
pp_check(msel.capts.MeanTSLF.re)
plot(conditional_effects(msel.capts.MeanTSLF.re), points = T)

# Time since last fire
msel.capts.TSLF.re <- brm(
  capts ~ TSLF + (1 | fieldtrip / plot / trap),
  data = Ajalapensis.captures.fire,
  family = "poisson",
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(msel.capts.TSLF.re)
bayes_R2(msel.capts.TSLF.re)
plot(conditional_effects(msel.capts.TSLF.re), points = T)


# Diagnostic plots
plot(msel.capts.TSLF.re)
pp_check(msel.capts.TSLF.re)

# Fire severity
msel.capts.severity.re <- brm(
  capts ~ severity + (1 | fieldtrip / plot / trap),
  data = Ajalapensis.captures.fire,
  family = "poisson",
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(msel.capts.severity.re)
bayes_R2(msel.capts.severity.re)

plot(conditional_effects(msel.capts.severity.re), points = T)

# Diagnostic plots
plot(msel.capts.severity.re)
pp_check(msel.capts.severity.re)

# Null model
msel.capts.firenull.re <- brm(
  capts ~ 1 + (1 | fieldtrip / plot / trap),
  data = Ajalapensis.captures.fire,
  family = "poisson",
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(msel.capts.firenull.re)
bayes_R2(msel.capts.firenull.re)

# Diagnostic plots
plot(msel.capts.firenull.re)
pp_check(msel.capts.firenull.re)

# Model comparisons
loo(
  msel.capts.MeanTSLF.re,
  msel.capts.TSLF.re,
  msel.capts.severity.re,
  msel.capts.firenull.re
)


# Conditional effects plot
p1 <- conditional_effects(msel.capts.severity.re)

ggplot(p1$severity, aes(x = severity, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = Ajalapensis.captures.fire,
    aes(x = jitter(severity), y = capts),
    alpha = 0.4,
    size = 2
  ) +
  labs(x = "Mean fire severity", y = "Number of captures") +
  theme_minimal()

p1 <- conditional_effects(msel.capts.MeanTSLF.re)

ggplot(p1$MeanTSLF, aes(x = MeanTSLF, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = Ajalapensis.captures.fire,
    aes(x = MeanTSLF, y = capts),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Mean fire interval", y = "Number of captures") +
  theme_minimal()

p1 <- conditional_effects(msel.capts.TSLF.re)

ggplot(p1$TSLF, aes(x = TSLF, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = Ajalapensis.captures.fire,
    aes(x = TSLF, y = capts),
    alpha = 0.5,
    size = 2
  ) +
  labs(x = "Time since last fire", y = "Number of captures")

### INLA------------------------------------------------------------------------
# -- Data Preparation for INLA --
# N-mixture models require a matrix of counts (sites x occasions) and covariates.

# Create the count matrix.
(y.mat <- as.matrix(Ajalapensis.captures.day.fire[, c(5:19)]))

# Create the `inla.mdata` object, which structures the data for the nmix family.
counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  Ajalapensis.captures.day.fire$fieldtrip,
  paste(
    Ajalapensis.captures.day.fire$fieldtrip,
    Ajalapensis.captures.day.fire$plot,
    sep = "_"
  ),
  paste(
    Ajalapensis.captures.day.fire$fieldtrip,
    Ajalapensis.captures.day.fire$plot,
    Ajalapensis.captures.day.fire$trap,
    sep = "_"
  ),
  scale(Ajalapensis.captures.day.fire$freq),
  scale(Ajalapensis.captures.day.fire$MeanTSLF),
  scale(Ajalapensis.captures.day.fire$severity),
  scale(Ajalapensis.captures.day.fire$TSLF)
)

print(counts.and.count.covs)

# Only intercept
out.inla.0 <- inla(
  counts.and.count.covs[, c(1:16)] ~ 1,
  data = list(counts.and.count.covs = counts.and.count.covs[, c(1:16)]),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(
        prior = "flat",
        param = numeric()
      )
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.0, digits = 3)


plot(
  out.inla.0,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = T
)


# Only random effects
out.inla.1 <- inla(
  counts.and.count.covs[, c(1:18)] ~ 1 + f(fieldtrip) + f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs[, c(1:18)],
    fieldtrip = counts.and.count.covs$X2,
    plot_id = counts.and.count.covs$X3
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.1, digits = 3)


plot(
  out.inla.1,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)


# Mean time since last fire effect

counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Ajalapensis.captures.day.fire$MeanTSLF),
  Ajalapensis.captures.day.fire$fieldtrip,
  paste(
    Ajalapensis.captures.day.fire$fieldtrip,
    Ajalapensis.captures.day.fire$plot,
    sep = "_"
  )
)

print(counts.and.count.covs)

out.inla.2 <- inla(
  counts.and.count.covs ~ 1 + MeanTSLF + f(fieldtrip) + f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    fieldtrip = counts.and.count.covs$X3,
    plot_id = counts.and.count.covs$X4,
    MeanTSLF = counts.and.count.covs$X2
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.2, digits = 3)


plot(
  out.inla.2,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

# Mean fire severity effect

counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Ajalapensis.captures.day.fire$severity),
  Ajalapensis.captures.day.fire$fieldtrip,
  paste(
    Ajalapensis.captures.day.fire$fieldtrip,
    Ajalapensis.captures.day.fire$plot,
    sep = "_"
  )
)

print(counts.and.count.covs)

out.inla.3 <- inla(
  counts.and.count.covs ~ 1 + severity + f(fieldtrip) + f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    fieldtrip = counts.and.count.covs$X3,
    plot_id = counts.and.count.covs$X4,
    severity = counts.and.count.covs$X2
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.3, digits = 3)


plot(
  out.inla.3,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

# Time since last fire effect

counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Ajalapensis.captures.day.fire$TSLF),
  Ajalapensis.captures.day.fire$fieldtrip,
  paste(
    Ajalapensis.captures.day.fire$fieldtrip,
    Ajalapensis.captures.day.fire$plot,
    sep = "_"
  )
)

print(counts.and.count.covs)

out.inla.4 <- inla(
  counts.and.count.covs ~ 1 + TSLF + f(fieldtrip) + f(plot_id),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    fieldtrip = counts.and.count.covs$X3,
    plot_id = counts.and.count.covs$X4,
    TSLF = counts.and.count.covs$X2
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.1)),
      theta2 = list(param = c(0, 0.1)),
      theta3 = list(param = c(0, 0.1)),
      theta4 = list(param = c(0, 0.1)),
      theta5 = list(param = c(0, 0.1)),
      theta6 = list(param = c(0, 0.1))
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)
summary(out.inla.4, digits = 3)


plot(
  out.inla.4,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = TRUE,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)

# Compare models
out.inla.0$waic$waic
out.inla.1$waic$waic
out.inla.2$waic$waic
out.inla.3$waic$waic
out.inla.4$waic$waic

inla_models <- list(
  out.inla.0,
  out.inla.1,
  out.inla.2,
  out.inla.3,
  out.inla.4
)

waic_values <- sapply(inla_models, function(model) model$waic$waic)
peff_values <- sapply(inla_models, function(model) model$waic$p.eff)
mlik_values <- sapply(inla_models, function(model) model$mlik[1])


delta_waic <- waic_values - min(waic_values)

w_akaike <- exp(-0.5 * delta_waic) / sum(exp(-0.5 * delta_waic))

waic_table <- data.frame(
  Model = paste0("Model ", 0:4),
  mlik = mlik_values,
  WAIC = waic_values,
  deltaWAIC = delta_waic,
  pWAIC = peff_values,
  Weight = round(w_akaike, 2)
)

waic_table$formula[1] <- "(lambda ~ 1), (p ~ 1)"
waic_table$formula[2] <- "(lambda ~ 1), (p ~ (1|fieldtrip/plot))"
waic_table$formula[
  3
] <- "(lambda ~ MeanTSLF), (p ~ MeanTSLF + (1|fieldtrip/plot))"
waic_table$formula[
  4
] <- "(lambda ~ FireSev), (p ~ FireSev + (1|fieldtrip/plot))"
waic_table$formula[5] <- "(lambda ~ TSLF), (p ~ TSLF + (1|fieldtrip/plot))"


# Sort by WAIC (best model first)
waic_table <- waic_table[order(waic_table$WAIC), ]

# Print table
print(waic_table)

summary(out.inla.3)
summary(out.inla.4)

summary(plogis(as.matrix(out.inla.3$summary.linear.predictor[, c(
  "mean",
  "sd",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "mode"
)])))
summary(out.inla.3$summary.fitted.values)

inla.nmix.lambda.fitted <- function(
  result,
  sample.size = 1000,
  return.posteriors = FALSE,
  scale = "exp"
) {
  fam <- result$.args$family
  if (length(grep(pattern = "nmix", x = fam)) == 0) {
    stop("This function is only for models with 'nmix' or 'nmixnb' likelihoods")
  }
  if (missing(result)) {
    stop("Please specify a model result")
  }
  s.check <- as.numeric(scale == "exp") + as.numeric(scale == "log")
  if (s.check == 0) {
    stop("Scale must be set to 'exp' or 'log'")
  }
  if (sample.size < 500) {
    warning("Please increase the sample size")
  }
  mdata.obj <- result$.args$data[[1]]
  counts <- as.data.frame(mdata.obj[grep(pattern = "Y", names(mdata.obj))])
  lambda.covs <- as.data.frame(mdata.obj[grep(
    pattern = "X",
    names(mdata.obj)
  )])
  lambda.covs <- as.matrix(lambda.covs[, -c(3, 4)])
  n.counts <- ncol(counts)
  n.lambda.covs <- ncol(lambda.covs)
  n.data <- nrow(counts)
  hyperpar.samples <- inla.hyperpar.sample(sample.size, result)
  s.names <- rownames(hyperpar.samples)
  if (fam == "nmixnb") {
    hyperpar.samples <- hyperpar.samples[, -c(3:ncol(hyperpar.samples))]
  }
  n.samp.covs <- ncol(hyperpar.samples)
  fitted.posteriors <- matrix(-1, nrow = n.data, ncol = sample.size)
  for (i in 1:n.data) {
    obs <- lambda.covs[i, ]
    for (j in 1:sample.size) {
      post <- hyperpar.samples[j, ]
      fitted <- sum(obs * post)
      fitted.posteriors[i, j] <- fitted
    }
  }
  index <- 1:n.data
  fitted.posteriors <- as.data.frame(fitted.posteriors)
  row.names(fitted.posteriors) <- NULL
  names(fitted.posteriors) <- s.names
  if (scale == "exp") {
    fitted.posteriors <- exp(fitted.posteriors)
  }
  fitted.meds <- round(
    apply(fitted.posteriors, 1, median),
    4
  )
  fitted.means <- round(
    apply(fitted.posteriors, 1, mean),
    4
  )
  fitted.sds <- round(apply(fitted.posteriors, 1, sd), 4)
  fitted.q025 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.025), 4)
  fitted.q500 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.5), 4)
  fitted.q975 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.975), 4)
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  fitted.modes <- round(
    apply(fitted.posteriors, 1, Mode),
    4
  )
  fitted.summary <- data.frame(
    mean.lambda = fitted.means,
    sd.lambda = fitted.sds,
    quant025.lambda = fitted.q025,
    median.lambda = fitted.q500,
    quant975.lambda = fitted.q975,
    mode.lambda = fitted.modes
  )
  fitted.summary <- cbind(index, fitted.summary)
  fitted.posteriors <- cbind(index, fitted.posteriors)
  if (return.posteriors == TRUE) {
    out <- list(
      fitted.summary = fitted.summary,
      fitted.posteriors = fitted.posteriors
    )
  } else {
    out <- list(fitted.summary = fitted.summary)
  }
  return(out)
}


# Severity effect
out.inla.3.lambda.fits <- inla.nmix.lambda.fitted(
  result = out.inla.3,
  sample.size = 10000,
  return.posteriors = FALSE
)$fitted.summary

ggplot(
  stack(out.inla.3.lambda.fits[, c(
    "mean.lambda",
    "median.lambda",
    "mode.lambda"
  )]),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "Abundance") +
  theme_classic()

ggplot(out.inla.3.lambda.fits, aes(x = median.lambda)) +
  geom_density(fill = "lightblue", alpha = 0.7, color = "black") +
  labs(
    y = "Density",
    x = "Abundance"
  ) +
  theme_classic()

ggplot(
  stack(log(out.inla.3.lambda.fits[, c(
    "mean.lambda",
    "median.lambda",
    "mode.lambda"
  )])),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "log(Abundance)") +
  theme_classic()

severity_xx <- seq(
  min(scale(Ajalapensis.captures.day.fire$severity)),
  max(scale(Ajalapensis.captures.day.fire$severity)),
  by = 0.1
)

# Draw 1000 samples from the model's posterior
posterior_samples <- inla.posterior.sample(1000, out.inla.3)

# The 'contents' table describes the structure of the model's output
contents <- out.inla.3$misc$configs$contents

# Find the numeric index for the intercept
intercept_idx <- contents$start[contents$tag == "(Intercept)"]

# Find the numeric index for the slope
slope_idx <- contents$start[contents$tag == "severity"]
# Initialize an empty list to store the results
predicted_lines_list <- list()

# Run the corrected sapply loop
predicted_lines <- sapply(posterior_samples, function(s) {
  # Access the column by its position [ ,1] instead of by name
  intercept <- s$latent[intercept_idx, 1]
  slope <- s$latent[slope_idx, 1]

  # The rest of the code is the same
  linear_predictor <- intercept + slope * severity_xx
  return(plogis(linear_predictor))
})

rel_p_severity <- data.frame(
  severity = severity_xx,
  median = apply(predicted_lines, 1, quantile, 0.5),
  lower = apply(predicted_lines, 1, quantile, 0.025),
  upper = apply(predicted_lines, 1, quantile, 0.975)
)

pred_df <- data.frame(
  severity = scale(Ajalapensis.captures.day.fire$severity),
  trap = Ajalapensis.captures.day.fire$trap,
  N = out.inla.3.lambda.fits$median.lambda,
  Nlow = out.inla.3.lambda.fits$quant025.lambda,
  Nhigh = out.inla.3.lambda.fits$quant975.lambda,
  p = plogis(out.inla.3$summary.linear.predictor$`0.5quant`),
  ncaps = rowSums(y.mat)
)

ggplot(pred_df, aes(x = severity, y = p)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_line(
    data = rel_p_severity,
    aes(x = severity, y = median),
    col = "blue",
    linewidth = 2
  ) +
  geom_ribbon(
    data = rel_p_severity,
    aes(ymin = lower, ymax = upper),
    alpha = 0.3
  ) +
  labs(x = "Mean fire severity", y = "Capture probability") +
  scale_y_log10(labels = number) +
  theme_minimal()

# Other simple plots
ggplot(pred_df, aes(x = ncaps, y = p)) +
  geom_point(size = 3, alpha = 0.3)

ggplot(pred_df, aes(x = ncaps, y = N)) +
  geom_point(size = 3, alpha = 0.3)

ggplot(pred_df, aes(x = p, y = N)) +
  geom_point(size = 3, alpha = 0.3)

ggplot(pred_df, aes(x = severity, y = N)) +
  geom_line(col = "blue") +
  geom_ribbon(data = pred_df, aes(ymin = Nlow, ymax = Nhigh), alpha = 0.3) +
  scale_y_log10()

ggplot(pred_df, aes(x = severity, y = ncaps)) +
  geom_point(size = 3, alpha = 0.3)

# Save workspace
save.image(file = "Ajalapensis_demography.RData")

# Map location ------------------------------------------------------------

library(tidyverse)
library(readxl)
library(geobr)
library(grid)
library(gt)
library(sp)
library(sf)
library(maptiles)
library(terra)
library(tidyterra)
library(spData)

# Load the multi-layer raster of burn dates
SGT_fire_r <- rast("SGT_Fire.tif")

# First, create a binary raster stack (1 if a fire occurred, 0 if not)
fire_binary_stack <- SGT_fire_r > 0

# Get the month for each layer from its name
raster_months <- as.integer(substr(names(SGT_fire_r), 6, 7))

# Use tapp() to group layers by month and sum them.
# This creates a 12-layer stack where each layer is the total fire
# frequency for that month across all years.
monthly_freq_stack <- tapp(
  fire_binary_stack,
  index = raster_months,
  fun = "sum",
  na.rm = TRUE
)

plot(monthly_freq_stack)

# Apply Seasonal Weights
# Define the weights for each of the 12 months
seasonal_weights <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3)

# Multiply the 12-layer frequency stack by the 12 weights
monthly_severity_stack <- monthly_freq_stack * seasonal_weights

plot(monthly_severity_stack)

# Calculate the Final Mean Severity ---
# As in your script, calculate the mean of the 12 monthly severity scores for each pixel
severity_raster <- app(monthly_severity_stack, "mean", na.rm = TRUE)
names(severity_raster) <- "severity"

# Optional: Save your new raster
writeRaster(severity_raster, "SGT_Severity.tif", overwrite = TRUE)

# Plot to check
plot(severity_raster)

# Read and plot points coordinates
traps.pts <- read.table("Points_Traps_Bruna_Heitor.txt", h = T)
summary(traps.pts)
traps.pts$sampdes <- c(
  rep("Rotating", 172),
  rep("Fixed", 48)
)

coordinates(traps.pts) <- c("Longitude", "Latitude")
proj4string(traps.pts) <- CRS("+proj=longlat +datum=WGS84")
traps.pts
traps.pts <- st_as_sf(traps.pts)

plot(severity_raster)
plot(traps.pts, add = T)

# Read biomes
biomes <- read_biomes()
biomes

Cerrado <- biomes[biomes$name_biome == "Cerrado", ]

# Read Conservation Units
CUs <- read_conservation_units()

SGT <- CUs[grep("SERRA GERAL DO TOCANTINS", CUs$name_conservation_unit), ]

# Read states
states <- read_state()

# Read South America
world <- spData::world
americas <- filter(world, region_un == "Americas")

ggplot() +
  geom_sf(data = biomes, color = NA, aes(fill = name_biome)) +
  geom_sf(data = states, color = "black", fill = NA) +
  scale_fill_manual(values = viridis::turbo(7)) +
  theme_void()

ggplot() +
  geom_sf(data = Cerrado, color = NA) +
  geom_sf(data = states, color = "black", fill = NA) +
  theme_void()

inset.cerrado <- ggplot() +
  geom_sf(data = americas) +
  geom_sf(data = Cerrado, color = NA, fill = "purple", alpha = 0.5) +
  geom_sf(data = states, fill = NA) +
  geom_sf(data = SGT, color = "forestgreen", fill = NA) +
  geom_rect(
    aes(xmin = -47.3, xmax = -45.8, ymin = -11.5, ymax = -10.4),
    color = "red",
    fill = NA
  ) +
  theme_test() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.background = element_rect(fill = "lightblue")
  ) +
  coord_sf(xlim = c(-74, -35), ylim = c(-34, 6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

inset.sgt <- ggplot() +
  geom_sf(data = SGT, color = "forestgreen", fill = NA, linewidth = 1.25) +
  # Add your trap locations
  geom_sf(
    data = traps.pts,
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c("Fixed" = 22, "Rotating" = 21),
    name = "Sampling Design"
  ) +
  geom_rect(
    aes(xmin = -47.3, xmax = -45.8, ymin = -11.5, ymax = -10.4),
    color = "red",
    fill = NA
  ) +
  theme_test() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.background = element_rect(fill = "lightgray")
  ) +
  coord_sf(xlim = c(-47.3, -45.8), ylim = c(-11.5, -10.4), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

# dowload tiles and compose raster (SpatRaster)
bbox.traps <- st_bbox(traps.pts)
bbox.traps[1:4] <- c(
  bbox.traps[1] - 0.05,
  bbox.traps[2] - 0.01,
  bbox.traps[3] + 0.05,
  bbox.traps[4] + 0.01
)

imagery <- get_tiles(
  bbox.traps,
  crop = TRUE,
  provider = "Esri.WorldImagery",
  zoom = 15
)
plot(imagery)
plot(st_geometry(traps.pts), col = "red", pch = 21, add = T)

map.traps.sgt <- ggplot() +
  geom_spatraster_rgb(data = imagery, interpolate = F) +
  # Add your trap locations
  geom_sf(
    data = traps.pts,
    size = 3,
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3
  ) +
  scale_shape_manual(
    values = c("Fixed" = 22, "Rotating" = 21),
    name = "Sampling Design"
  ) +
  xlim(bbox.traps[1], bbox.traps[3]) +
  ylim(bbox.traps[2], bbox.traps[4]) +
  theme_test() +
  guides(size = "none") +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white"),
    pad_y = unit(0.2, "in")
  ) +
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.53, "in"),
    pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +
  coord_sf(
    xlim = c(bbox.traps[1], bbox.traps[3]),
    ylim = c(bbox.traps[2], bbox.traps[4]),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude")


map.traps.sgt

# Combining both maps
print(inset.cerrado, vp = viewport(0.7, 0.8, width = 0.2, height = 0.2))
print(inset.sgt, vp = viewport(0.68, 0.6, width = 0.2, height = 0.2))

# Map with fire severity

#--- Create the Final Map ---
map.traps.severity <- ggplot() +
  # Use geom_spatraster to plot your severity raster
  geom_spatraster(data = severity_raster) +

  # Add a color scale for severity (e.g., "inferno" or "viridis")
  scale_fill_viridis_c(
    option = "inferno",
    name = "Fire Severity",
    na.value = "transparent"
  ) +

  # Add your trap locations
  geom_sf(
    data = traps.pts,
    size = 3,
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3
  ) +
  scale_shape_manual(
    values = c("Fixed" = 22, "Rotating" = 21),
    name = "Sampling Design"
  ) +

  # Add scale bar and north arrow
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white")
  ) +
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.53, "in"),
    pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +

  # Set coordinates and labels
  coord_sf(
    xlim = c(bbox.traps[1], bbox.traps[3]),
    ylim = c(bbox.traps[2], bbox.traps[4]),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()

# Display the map
print(map.traps.severity)

print(map.traps.severity)

# Combining both maps
print(inset.cerrado, vp = viewport(0.7, 0.8, width = 0.2, height = 0.2))
print(inset.sgt, vp = viewport(0.68, 0.6, width = 0.2, height = 0.2))

# Power analysis (PJS) -----------------------------------------------
## 0. Setup -------------------------------------------------------
library(openCR)
library(parallel)
set.seed(123)

## 1. Define empirical parameters (from your fitted model) ------------------------
plots <- c("A1", "A2", "A3", "A4")

N_true <- c(25, 25, 25, 25)
names(N_true) <- plots

p_true <- c(0.007, 0.011, 0.027, 0.016)
names(p_true) <- plots

phi_true <- c(0.75, 0.84, 0.72, 0.91)

## 2. Match your sampling design ---------------------------------------------------
nsessions <- 4 # field trips
noccasions <- 15 # days per field trip

## 3. Function to simulate one plot ------------------------------------------------

sim_one_plot <- function(N, phi, p) {
  sim.nonspatial(
    N = N,
    turnover = list(phi = phi, lambda = 1),
    p = p,
    nsessions = nsessions,
    noccasions = noccasions
  )
}

## 4. Simulate one full dataset (4 plots combined) -------------------------------
sim_one_dataset <- function(phi, p, N = 400) {
  sim.nonspatial(
    N = N,
    turnover = list(phi = phi, lambda = 1),
    p = p,
    nsessions = 4,
    noccasions = 15
  )
}


## 5. Fit the same model you used (Pradel / JSSAfCL) ------------------------------
fit_pradel <- function(ch) {
  openCR.fit(
    ch,
    type = "JSSAfCL",
    model = list(
      p ~ 1,
      phi ~ 1,
      f ~ 1
    )
  )
}

## 6. Define a power criterion --------------------------------------
detect_p_difference <- function(fit) {
  est <- predict(fit)$p
  est <- est[!is.na(est$estimate), ]

  combs <- combn(plots, 2, simplify = FALSE)

  for (cc in combs) {
    a <- est[est$Plot == cc[1], ]
    b <- est[est$Plot == cc[2], ]

    overlap <- !(a$ucl < b$lcl | b$ucl < a$lcl)
    if (!overlap) return(TRUE)
  }

  FALSE
}

phi_estimable <- function(fit) {
  est <- predict(fit)$phi
  est <- est[!is.na(est$estimate), ]

  if (nrow(est) == 0) {
    return(FALSE)
  }

  (est$ucl - est$lcl) < 0.2
}

## 7. Run the simulation (power analysis) ------------------------------
phi_vals <- c(0.75, 0.85, 0.90)
p_vals <- c(0.007, 0.015, 0.03)

nrepl <- 1000

out <- expand.grid(phi = phi_vals, p = p_vals)
out$power <- NA

for (i in seq_len(nrow(out))) {
  res <- replicate(nrepl, {
    ch <- sim_one_dataset(phi = out$phi[i], p = out$p[i])
    fit <- try(fit_pradel(ch), silent = TRUE)

    if (inherits(fit, "try-error")) {
      return(FALSE)
    }
    phi_estimable(fit)
  })

  out$power[i] <- mean(res)
}


## 8. Power estimate -----------------------------------------------

out
