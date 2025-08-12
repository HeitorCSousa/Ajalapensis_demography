rm(list = ls())

library(tidyverse)
#library(Amelia)
library(missForest)
library(Boruta)
library(GGally)
library(viridis)
library(ggfortify)


setwd("~/Documents/UFT/Geisa")
#Environmental variables
env.vars <- readRDS("env_vars.rds")
env.vars$plot <- factor(env.vars$plot, levels = c("C2",
                                                  "C1",
                                                  "QP",
                                                  "QT"))
env.vars <- as.data.frame(env.vars)

env.vars <- do.call(data.frame,                      # Replace Inf in data by NA
                   lapply(env.vars,
                          function(x) replace(x, is.infinite(x), NA)))

summary(env.vars)


#Ecophysiology#
###############

#Temperaturas criticas#
#######################
#setwd("/Volumes/Extreme SSD/Heitor/Documentos/UFT/Discentes/Raquel_Acacio/Analises")

#Le tabela
ct.sgt <- read.table("Dados_CTs.txt", h=T)

ct.Ajalapensis <- ct.sgt[ct.sgt$Especie=="Ameivula_jalapensis",]
ct.Ajalapensis

ct.Toreadicus <- ct.sgt[ct.sgt$Especie=="Tropidurus_oreadicus",]
ct.Toreadicus


table(ct.Ajalapensis$Sexo)

#Desempenho locomotor#
######################
des.loc <- read.table("desemp_loc_SGT.txt",h=T)
str(des.loc)
head(des.loc)

table(des.loc$sp,des.loc$SGT)

summary(abs(des.loc$v))
summary(abs(des.loc$a))
#windows(10,10)
boxplot(abs(des.loc$v))

dados <- aggregate(abs(des.loc$v),by=list(des.loc$sp,
                                          des.loc$SGT,
                                          des.loc$Temp,
                                          des.loc$Corrida,
                                          des.loc$Sexo,
                                          des.loc$CRC),
                   FUN=max, na.rm=T)

head(dados)
names(dados)<-c("especie","SGT","temp",
                "corrida","sexo","CRC",
                "Veloc")

head(dados)
boxplot(Veloc ~ especie, data = dados)

dados.ctmin <- data.frame("especie"=ct.sgt$Especie,
                          "SGT"=ct.sgt$ID,
                          "temp"=ct.sgt$Ctmin,
                          "corrida"=c(rep(NA,length(ct.sgt$Ctmin))),
                          "sexo"=ct.sgt$Sexo,
                          "CRC"=ct.sgt$CRC,
                          "Veloc"=c(rep(0,length(ct.sgt$Ctmin))))

dados.ctmax <- data.frame("especie"=ct.sgt$Especie,
                          "SGT"=ct.sgt$ID,
                          "temp"=ct.sgt$Ctmax,
                          "corrida"=c(rep(NA,length(ct.sgt$Ctmin))),
                          "sexo"=ct.sgt$Sexo,
                          "CRC"=ct.sgt$CRC,
                          "Veloc"=c(rep(0,length(ct.sgt$Ctmin))))

dados.completos <- rbind(dados,
                         dados.ctmin,
                         dados.ctmax)

head(dados.completos)
tail(dados.completos)
table(dados.completos$especie)

#Standardize species' names
dados.completos$especie[dados.completos$especie=="A_ameiva"] <- "Ameiva_ameiva"
dados.completos$especie[dados.completos$especie=="A_jalapensis"] <- "Ameivula_jalapensis"
dados.completos$especie[dados.completos$especie=="B_heathi"] <- "Brasiliscincus_heathi"
dados.completos$especie[dados.completos$especie=="B_oxyrhina"] <- "Bachia_oxyrhina"
dados.completos$especie[dados.completos$especie=="B_oyrhina"] <- "Bachia_oxyrhina"
dados.completos$especie[dados.completos$especie=="C_nigropunctatum"] <- "Copeoglossum_nigropunctatum"
dados.completos$especie[dados.completos$especie=="H_brasilianus"] <- "Hemidactylus_brasilianus"
dados.completos$especie[dados.completos$especie=="T_oreadicus"] <- "Tropidurus_oreadicus"
dados.completos$especie[dados.completos$especie=="V_savanicola"] <- "Vanzosaura_savanicola"

table(dados.completos$especie)


#Curvas de desempenho locomotor#
################################
#install.packages("gamm4",dep=T)
library(gamm4)
library(brms)
library(ggplot2)

# options(brms.backend = "cmdstanr")


# detach(des.loc)

#Ameivula jalapensis
dados.Ajalapensis <- dados.completos[dados.completos$especie=="Ameivula_jalapensis",]

dados.Ajalapensis$sexo <- as.factor(dados.Ajalapensis$sexo)
summary(dados.Ajalapensis$sexo)
#Modelo generalizado aditivo de efeitos mistos
m.Ajalapensis.sprint <- gamm4(Veloc~t2(temp),
                          random=~(1|SGT),
                          data=dados.completos[dados.completos$especie=="Ameivula_jalapensis",])
summary(m.Ajalapensis.sprint$gam)

m.Ajalapensis.sex.sprint <- gamm4(Veloc~t2(temp, by = sexo),
                              random=~(1|SGT),
                              data=dados.Ajalapensis)

summary(m.Ajalapensis.sex.sprint$gam)
plot(m.Ajalapensis.sex.sprint$gam)

AIC(m.Ajalapensis.sprint$mer)
AIC(m.Ajalapensis.sex.sprint$mer)

m.Ajalapensis.sprint.refit <- gamm4(Veloc~t2(temp),
                              random=~(1|SGT),
                              data=dados.Ajalapensis[!is.na(dados.Ajalapensis$sexo),])
anova(m.Ajalapensis.sprint.refit$mer, m.Ajalapensis.sex.sprint$mer)

m.Ajalapensis.sprint.brms <- brm(Veloc~t2(temp) + (1|SGT),
                              data=dados.completos[dados.completos$especie=="Ameivula_jalapensis",])
summary(m.Ajalapensis.sprint.brms)
plot(m.Ajalapensis.sprint.brms)
plot(conditional_effects(m.Ajalapensis.sprint.brms), points = T, 
     xlab = "Temperature (°C)", ylab="Predicted Speed (m/s)")

p1 <- conditional_effects(m.Ajalapensis.sprint.brms)
#windows(10,10)
quartz(8,8)
ggplot(data = dados.completos[dados.completos$especie=="Ameivula_jalapensis",],
       aes(x = temp, y = Veloc)) +
  geom_point(alpha = 0.5) +
  geom_line(data = p1$temp, aes(x = temp, y = estimate__), colour = "blue", linewidth = 1.2) +
  geom_ribbon(data = p1$temp, aes(ymin=lower__, ymax=upper__), linetype=3, alpha=.2) +
  labs(y='Predicted Speed (m/s)',x='Temperature (°C)') +
  theme(plot.title = element_text(hjust = 0.5))
     
#windows(10,10)
quartz(8,8)
plot(m.Ajalapensis.sprint$gam,
     main=expression(italic("Ameivula jalapensis")),
     ylab="Desempenho locomotor",xlab="Temperatura (°C)")


preddata_tpc <- data.frame(temp=seq(10,50,0.1))

## Faz a predicao 
pred_tpc_Ajalapensis <- predict(m.Ajalapensis.sprint$gam,
                                preddata_tpc,
                                se.fit=T)

## Anexa as predicoes e erros na tabela de dados
preddata_tpc <- rbind(preddata_tpc)
preddata_tpc$predicted <- c(pred_tpc_Ajalapensis$fit)
preddata_tpc$se <- c(pred_tpc_Ajalapensis$se.fit)

quartz(h = 8, w = 8)
ggplot(preddata_tpc, aes(x=temp, y=predicted)) +  geom_line(color = "blue", linewidth = 1.5) + 
  geom_ribbon(aes(ymin=predicted - se, ymax=predicted+se), linetype=3, alpha=.2)+
  geom_point(data = dados.completos[dados.completos$especie=="Ameivula_jalapensis",],
         aes(x = temp, y = Veloc), alpha = 0.5)+
  labs(y='Predicted Speed (m/s)',x='Temperature (°C)') +
  lims(x=c(15,45),y=c(0,1.8)) + theme(plot.title = element_text(hjust = 0.5))

#Temperatura ótima
preddata_tpc[preddata_tpc$predicted==max(pred_tpc_Ajalapensis$fit),]
p1$temp[p1$temp$estimate__==max(p1$temp$estimate__),]
p1$temp[p1$temp$upper__==max(p1$temp$upper__),]
p1$temp[p1$temp$lower__==max(p1$temp$lower__),]


#31.7 for A. jalapensis

#Prediz com os dados microclimaticos
microclim.SGT <- readRDS("microclimate_EESGT.rds")

env.var.hour <- microclim.SGT %>% 
  group_by(plot, arm.int, campanha, hour, day, month, year) %>% 
  summarise(tmed = mean(temp, na.rm = T),
            rhmed = mean(rh, na.rm = T))

env.var.hour$Ajalapensis_perf <- predict.gam(m.Ajalapensis.sprint$gam,
                                newdata = data.frame(temp = env.var.hour$tmed))

#Tropidurus oreadicus
#Modelo generalizado aditivo de efeitos mistos
m.Toreadicus.sprint <- gamm4(Veloc~t2(temp),
                              random=~(1|SGT),
                              data=dados.completos[dados.completos$especie=="Tropidurus_oreadicus",])
summary(m.Toreadicus.sprint$gam)

m.Toreadicus.sprint.brms <- brm(Veloc~t2(temp) + (1|SGT),
                                 data=dados.completos[dados.completos$especie=="Tropidurus_oreadicus",])
summary(m.Toreadicus.sprint.brms)
plot(m.Toreadicus.sprint.brms)
plot(conditional_effects(m.Toreadicus.sprint.brms), points = T, 
     xlab = "Temperature (°C)", ylab="Predicted Speed (m/s)")

p2 <- conditional_effects(m.Toreadicus.sprint.brms)
#windows(10,10)
quartz(8,8)
ggplot(data = dados.completos[dados.completos$especie=="Tropidurus_oreadicus",],
       aes(x = temp, y = Veloc)) +
  geom_point(alpha = 0.5) +
  geom_line(data = p2$temp, aes(x = temp, y = estimate__), colour = "blue", linewidth = 1.2) +
  geom_ribbon(data = p2$temp, aes(ymin=lower__, ymax=upper__), linetype=3, alpha=.2) +
  labs(y='Predicted Speed (m/s)',x='Temperature (°C)') +
  theme(plot.title = element_text(hjust = 0.5))

#windows(10,10)
quartz(8,8)
plot(m.Toreadicus.sprint$gam,
     main=expression(italic("Tropidurus oreadicus")),
     ylab="Locomotor performance",xlab="Temperature (°C)")


preddata_tpc <- data.frame(temp=seq(10,50,0.1))

## Faz a predicao 
pred_tpc_Toreadicus <- predict(m.Toreadicus.sprint$gam,
                                preddata_tpc,
                                se.fit=T)

## Anexa as predicoes e erros na tabela de dados
preddata_tpc <- rbind(preddata_tpc)
preddata_tpc$predicted <- c(pred_tpc_Toreadicus$fit)
preddata_tpc$se <- c(pred_tpc_Toreadicus$se.fit)

quartz(8,8)
ggplot(preddata_tpc, aes(x=temp, y=predicted)) +  geom_line() + 
  geom_ribbon(aes(ymin=predicted - se, ymax=predicted+se), linetype=3, alpha=.2)+
  geom_point(data = dados.completos[dados.completos$especie=="Tropidurus_oreadicus",],
             aes(x = temp, y = Veloc), alpha = 0.5)+
  labs(y='Predicted Speed (m/s)',x='Temperature (°C)') +
  lims(x=c(10,50),y=c(0,5)) + theme(plot.title = element_text(hjust = 0.5))

#Temperatura ótima
preddata_tpc[preddata_tpc$predicted==max(pred_tpc_Toreadicus$fit),]
p2$temp[p2$temp$estimate__==max(p2$temp$estimate__),]
p2$temp[p2$temp$upper__==max(p2$temp$upper__),]
p2$temp[p2$temp$lower__==max(p2$temp$lower__),]


#30.4 for T. oreadicus

#Prediz com os dados microclimaticos
microclim.SGT <- readRDS("microclimate_EESGT.rds")

env.var.hour <- microclim.SGT %>% 
  group_by(plot, arm.int, campanha, hour, day, month, year) %>% 
  summarise(tmed = mean(temp, na.rm = T),
            rhmed = mean(rh, na.rm = T))

env.var.hour$Ajalapensis_perf <- predict.gam(m.Ajalapensis.sprint$gam,
                                             newdata = data.frame(temp = env.var.hour$tmed))

env.var.hour$Toreadicus_perf <- predict.gam(m.Toreadicus.sprint$gam,
                                             newdata = data.frame(temp = env.var.hour$tmed))

saveRDS(m.Ajalapensis.sprint, "tpc_Ajalapensis.rds")
saveRDS(m.Toreadicus.sprint, "tpc_Toreadicus.rds")

#Temperatura preferencial#
##########################
tpref_SGT <- read.table("Dados_Tpref_SGT.txt", h=T)

tpref_Ajalapensis <- dplyr::filter(tpref_SGT, sp == "A_jalapensis")
table(tpref_Ajalapensis$SGT,tpref_Ajalapensis$sp)
nrow(table(tpref_Ajalapensis$SGT,tpref_Ajalapensis$sp))

tpref_Toreadicus <- dplyr::filter(tpref_SGT, sp == "T_oreadicus")
table(tpref_Toreadicus$SGT,tpref_Toreadicus$sp)
nrow(table(tpref_Toreadicus$SGT,tpref_Toreadicus$sp))

hvtFUN <- function(temp.envr,temp.lab, quantiles, radiation) {
  vtmin <- quantile(temp.lab, quantiles[1], na.rm = TRUE)
  vtmax <- quantile(temp.lab, quantiles[2], na.rm = TRUE)
  hv <- ifelse(temp.envr > vtmin & temp.envr < vtmax, 1, 0)
  hv[radiation == 0] <- 0
  hv
}

(tpref90_Ajalapensis <- quantile(tpref_Ajalapensis$temp, c(0.05, 0.95), na.rm=T))
(tpref50_Ajalapensis <- quantile(tpref_Ajalapensis$temp, c(0.25, 0.75), na.rm=T))

(tpref90_Toreadicus <- quantile(tpref_Toreadicus$temp, c(0.05, 0.95), na.rm=T))
(tpref50_Toreadicus <- quantile(tpref_Toreadicus$temp, c(0.25, 0.75), na.rm=T))

saveRDS(tpref90_Ajalapensis,
        "tpref90_Ajalapensis.rds")

saveRDS(tpref50_Ajalapensis,
        "tpref50_Ajalapensis.rds")

saveRDS(tpref90_Toreadicus,
        "tpref90_Toreadicus.rds")

saveRDS(tpref50_Toreadicus,
        "tpref50_Toreadicus.rds")

env.var.hour$Ajalapensis_ha90 <- hvtFUN(env.var.hour$tmed,
                           tpref_Ajalapensis$temp,
                           c(0.05,0.95),
                           rep(1, nrow(env.var.hour)))

env.var.hour$Ajalapensis_ha50 <- hvtFUN(env.var.hour$tmed,
                           tpref_Ajalapensis$temp,
                           c(0.25,0.75),
                           rep(1, nrow(env.var.hour)))

env.var.hour$Toreadicus_ha90 <- hvtFUN(env.var.hour$tmed,
                                        tpref_Toreadicus$temp,
                                        c(0.05,0.95),
                                        rep(1, nrow(env.var.hour)))

env.var.hour$Toreadicus_ha50 <- hvtFUN(env.var.hour$tmed,
                                        tpref_Toreadicus$temp,
                                        c(0.25,0.75),
                                        rep(1, nrow(env.var.hour)))

ecophys.day <-
  env.var.hour %>%
  dplyr::group_by(plot, arm.int, campanha, day, month, year) %>% 
  dplyr::summarise(Ajalapensis_perf = mean(Ajalapensis_perf),
                   Ajalapensis_ha50 = sum(Ajalapensis_ha50),
                   Ajalapensis_ha90 = sum(Ajalapensis_ha90),
                   Toreadicus_perf = mean(Toreadicus_perf),
                   Toreadicus_ha50 = sum(Toreadicus_ha50),
                   Toreadicus_ha90 = sum(Toreadicus_ha90)
  )

summary(ecophys.day)
pts.arms <- read.table("Pontos_Armadilhas.txt", h = T)
pts.arms.SGT <- pts.arms[pts.arms$local=="EESGT",]

daylength <- geosphere::daylength(lat = mean(pts.arms.SGT$lat),
                                  doy = yday(paste(ecophys.day$year,
                                                   sprintf("%02d", as.numeric(ecophys.day$month)),
                                                   sprintf("%02d", as.numeric(ecophys.day$day)), 
                                                   sep = "-")))

ecophys.day$Ajalapensis_ha50 <- ifelse(ecophys.day$Ajalapensis_ha50 > daylength,
                                       daylength,
                                       ecophys.day$Ajalapensis_ha50)

ecophys.day$Ajalapensis_ha90 <- ifelse(ecophys.day$Ajalapensis_ha90 > daylength,
                                       daylength,
                                       ecophys.day$Ajalapensis_ha90)

ecophys.day$Toreadicus_ha50 <- ifelse(ecophys.day$Toreadicus_ha50 > daylength,
                                       daylength,
                                       ecophys.day$Toreadicus_ha50)

ecophys.day$Toreadicus_ha90 <- ifelse(ecophys.day$Toreadicus_ha90 > daylength,
                                       daylength,
                                       ecophys.day$Toreadicus_ha90)

ecophys.month <-
  ecophys.day %>%
  dplyr::group_by(plot, arm.int, campanha) %>% 
  dplyr::summarise(Ajalapensis_perf = mean(Ajalapensis_perf),
                   Ajalapensis_ha50 = mean(Ajalapensis_ha50),
                   Ajalapensis_ha90 = mean(Ajalapensis_ha90),
                   Toreadicus_perf = mean(Toreadicus_perf),
                   Toreadicus_ha50 = mean(Toreadicus_ha50),
                   Toreadicus_ha90 = mean(Toreadicus_ha90)
  )

ecophys.month$plot <- factor(ecophys.month$plot, levels = c("C2",
                                                  "C1",
                                                  "QP",
                                                  "QT"))

env.vars <- left_join(env.vars, ecophys.month, by = c("plot", "arm.int", "campanha"))
summary(env.vars)

#Open cluster to speed up computation
cl <- parallel::makeCluster(ncol(env.vars)-1,
                            setup_timeout = 0.5)
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()

env.vars$plot.int <- as.integer(env.vars$plot)

#Imputation with missforest
env.vars.imputed <- missForest::missForest(as.matrix(env.vars[,-1]),
                                           verbose = T,
                                           variablewise = T,
                                           maxiter = 100,
                                           ntree = 1000,
                                           parallelize = "variables")
env.vars.imputed$OOBerror

names(env.vars[,-1])[env.vars.imputed$OOBerror>1]#Varables with error higher than 1
names(env.vars[,-1])[env.vars.imputed$OOBerror>5]#Varables with error higher than 5

env.vars.imputed.df <- as.data.frame(env.vars.imputed$ximp)
env.vars.imputed.df <- cbind(env.vars[,1], env.vars.imputed.df)

names(env.vars.imputed.df) <- c("plot", names(env.vars.imputed.df)[-1])
summary(env.vars.imputed.df)
saveRDS(env.vars.imputed.df, "env_vars_imputed_df.rds")

#Which variables better explain the differences among plots?
set.seed(123)
Borutaplot <- Boruta(plot ~., data = env.vars.imputed.df[,-c(2,3,15:21)], doTrace = 2,
                     maxRuns =500)

Borutaplot
saveRDS(Borutaplot, "Borutaplot.rds")
Borutaplot <- readRDS("Borutaplot.rds")
plot(Borutaplot, las = 2)

windows(height = 8, width = 12)
plot(Borutaplot, las = 2, bty = "n", xlab="", xaxt = "n")

lz <- lapply(1:ncol(Borutaplot$ImpHistory), function(i) 
  Borutaplot$ImpHistory[is.finite(Borutaplot$ImpHistory[,i]),i])
names(lz) <- colnames(Borutaplot$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(side = 1, las = 2, labels = names(Labels),
     at = 1:ncol(Borutaplot$ImpHistory), cex.axis = 0.7)

plotImpHistory(Borutaplot)
attStats(Borutaplot)
TentativeRoughFix(Borutaplot)

boxplot(env.vars.imputed.df$pground)

ggpairs.env <- ggpairs(env.vars.imputed.df[env.vars.imputed.df$pground>70,-c(2,3,21)])
quartz(height = 8, width = 12)
ggpairs.env

ggplot(env.vars.imputed.df, 
       aes(x = plot, y = VARI.all)) + 
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
  coord_cartesian( clip = "off")+
  labs(x="", y="Vegetation index")

ggplot(env.vars.imputed.df, 
       aes(x = plot, y = tree.density)) + 
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
  coord_cartesian( clip = "off")+
  labs(x="", y="Tree density")

#Organize fire info
#Time since last fire (TSLF)
fire.arms.df <- read.csv("fire_arms_df.csv")
fire.arms.df <- filter(fire.arms.df,ID>172)
fire.arms.df <- fire.arms.df[,-1]


fire.arms.df$trap <- as.factor(fire.arms.df$regime)
env.vars.imputed.df$trap <- paste(env.vars.imputed.df$plot, 
                                     env.vars.imputed.df$arm.int, sep = "_")
env.vars.imputed.df$trap <- factor(env.vars.imputed.df$trap,
                            levels = levels(fire.arms.df$trap))
levels(fire.arms.df$trap) == levels(env.vars.imputed.df$trap)
env.vars.imputed.df <- env.vars.imputed.df[order(env.vars.imputed.df$campanha),]
env.vars.imputed.df$year <- rep(c(2021,2021,2022,2022),each=48)
env.vars.imputed.df$month <- rep(c(3,7,2,6),each=48)

env.data <- left_join(env.vars.imputed.df, fire.arms.df, by=c("plot", 
                                                           "trap", 
                                                           "year", 
                                                           "month"))
head(env.data)
summary(env.data)
tapply(env.data$TSLF,INDEX = list(env.data$trap,env.data$camp),FUN = min)

pal.plot <- turbo(4)

#MODIS made some errors about fires during our sampling. We will correct with our field information
env.data$TSLF[env.data$plot=="C1" & env.data$campanha==1] <- min(env.data$TSLF[env.data$plot=="C1" & env.data$campanha==1])
env.data$TSLF[env.data$plot=="C1" & env.data$campanha==2] <- min(env.data$TSLF[env.data$plot=="C1" & env.data$campanha==2])
env.data$TSLF[env.data$plot=="C1" & env.data$campanha==3] <- min(env.data$TSLF[env.data$plot=="C1" & env.data$campanha==3])
env.data$TSLF[env.data$plot=="C1" & env.data$campanha==4] <- min(env.data$TSLF[env.data$plot=="C1" & env.data$campanha==4])

#Last burn in C2 was in June 2020
env.data$TSLF[env.data$plot=="C2" & env.data$campanha==1] <- 8
env.data$TSLF[env.data$plot=="C2" & env.data$campanha==2] <- 8+3
env.data$TSLF[env.data$plot=="C2" & env.data$campanha==3] <- 8+12
env.data$TSLF[env.data$plot=="C2" & env.data$campanha==4] <- 8+15

#Last burn in QP was in May 2020
env.data$TSLF[env.data$plot=="QP" & env.data$campanha==1] <- 9
env.data$TSLF[env.data$plot=="QP" & env.data$campanha==2] <- 9+3
env.data$TSLF[env.data$plot=="QP" & env.data$campanha==3] <- 9+12
env.data$TSLF[env.data$plot=="QP" & env.data$campanha==4] <- 1

#Last burn in QP was in May 2020
env.data$TSLF[env.data$plot=="QT" & env.data$campanha==1] <- min(env.data$TSLF[env.data$plot=="QT" & env.data$campanha==1])
env.data$TSLF[env.data$plot=="QT" & env.data$campanha==2] <- min(env.data$TSLF[env.data$plot=="QT" & env.data$campanha==1])+3
env.data$TSLF[env.data$plot=="QT" & env.data$campanha==3] <- min(env.data$TSLF[env.data$plot=="QT" & env.data$campanha==1])+12
env.data$TSLF[env.data$plot=="QT" & env.data$campanha==4] <- 0.5

#Fire severity index
#Fires in mid dry season has weight = 2
#Fires in late dry season has weight = 3
#Other fires = 1
fire.regimes <- read.csv("fire_regimes_arms_df.csv")
seq.fire <- 173:220
(severity.C1 <- rep(mean(fire.regimes$severity[seq.fire[1:12]]),48))
(severity.C2 <- rep(mean(fire.regimes$severity[seq.fire[13:24]]),48))
(severity.QP <- rep(mean(fire.regimes$severity[seq.fire[25:36]]),48))
(severity.QT <- rep(mean(fire.regimes$severity[seq.fire[37:48]]),48))


env.data <- env.data[order(env.vars.imputed.df$plot),]
env.data$plot

env.data$severity <- c(severity.C1,severity.C2,severity.QP,severity.QT)

#Some plots to see relationship of TSLF with environmental data
plot(t.med~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(t.max~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(t.max.abs~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(rh.max.abs~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(rh.min.abs~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(VARI.all~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(tree.density~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])
plot(zentropy~TSLF,data=env.data,col=pal.plot[as.numeric(as.factor(env.data$plot))])


pca.env.all<-prcomp(env.data[,c("TSLF","t.med","t.max","t.max.abs",
                            "rh.min.abs", "rh.max.abs",
                            "VARI.all","zentropy","zkurt","pzabovezmean", "pground","tree.density")],
                scale=T)
summary(pca.env.all)

pca.env.all
plot(pca.env.all)

#windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c("C2",
                                                  "C1",
                                                  "QP",
                                                  "QT"))
env.data$campanha <- as.factor(env.data$campanha)
windows(8,8)
plot.pca.env.all <- autoplot(pca.env.all, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3)
plot.pca.env.all

#Just with less correlated variables
pca.env<-prcomp(env.data[,c("TSLF","t.med","t.max","t.max.abs",
                                "rh.min.abs", "rh.max.abs",
                                "VARI.all", "zentropy", "tree.density")],
                    scale=T)
summary(pca.env)

pca.env
plot(pca.env)

#windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c("C2",
                                                  "C1",
                                                  "QP",
                                                  "QT"))
env.data$campanha <- as.factor(env.data$campanha)
windows(8,8)
plot.pca.env <- autoplot(pca.env, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3)
plot.pca.env

windows(height = 12, width = 8)
gridExtra::grid.arrange(plot.pca.env.all, plot.pca.env)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)
autoplot(pca.env, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3)

windows(height = 8, width = 10)
autoplot(pca.env, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3,x=3,y=4)

#Including ecophysiological variables
pca.env.ecophys <- prcomp(env.data[,c("TSLF","t.med","t.max","t.max.abs",
                            "rh.min.abs", "rh.max.abs",
                            "VARI.all", "zentropy", "tree.density",
                            "Ajalapensis_perf", "Ajalapensis_ha90")],
                scale=T)
summary(pca.env.ecophys)

pca.env
plot(pca.env.ecophys)

#windows(30,30)
env.data$plot <- factor(env.data$plot, levels = c("C2",
                                                  "C1",
                                                  "QP",
                                                  "QT"))
env.data$campanha <- as.factor(env.data$campanha)
windows(8,8)
plot.pca.env.ecophys <- autoplot(pca.env.ecophys, data = env.data, colour = 'plot',shape='campanha',frame=T,
                         loadings=T,loadings.label = TRUE,loadings.label.size = 3)
plot.pca.env.ecophys

windows(height = 12, width = 8)
quartz(height = 12, width = 8)

gridExtra::grid.arrange(plot.pca.env.ecophys, plot.pca.env)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)
autoplot(pca.env.ecophys, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3)

windows(height = 8, width = 10)
quartz(height = 8, width = 10)

autoplot(pca.env.ecophys, data = env.data, colour = 'plot',shape='campanha',frame=T,
         loadings=T,loadings.label = TRUE,loadings.label.size = 3,x=3,y=4)

env.data <- env.data[,c("plot","trap","campanha","severity","TSLF","t.med","t.max","t.max.abs",
                        "rh.min.abs", "rh.max.abs","VARI.all", "zentropy", "tree.density",
                        "Ajalapensis_perf", "Ajalapensis_ha90",
                        "Toreadicus_perf", "Toreadicus_ha90")]

quartz(height = 8, width = 8)
ggplot(env.data, 
       aes(x = plot, y = Ajalapensis_perf)) + 
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
  coord_cartesian( clip = "off")+
  labs(x="", y="Locomotor performance (m/s)")

quartz(height = 8, width = 8)
ggplot(env.data, 
       aes(x = plot, y = Ajalapensis_ha90)) + 
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
  coord_cartesian( clip = "off")+
  labs(x="", y="Hours of activity")

quartz(height = 8, width = 12)
ggpairs.env <- ggpairs(env.data[,-c(2:4)])
ggpairs.env

saveRDS(env.data, "env_data_SGT.rds")
