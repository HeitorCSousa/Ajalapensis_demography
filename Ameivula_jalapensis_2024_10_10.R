
# Set working directory and load packages ---------------------------------

rm(list = ls())

#Heitor
# getwd()
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

install.packages("INLA",repos=c(getOption("repos"),
                                INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)


##############################################
## Importa dados de capturas, recaps e CRC
##############################################

dados <- read.table("Ameivula_jalapensis_04.txt", h=T)
dados$Parcela <- factor(dados$Parcela, levels = c("C2", "C1", "QP", "QT"))
attach(dados)
head(dados)
tail(dados)
str(dados)
# detach(dados)

## Calcula frequencias de capturas e recaps, faz graficos e testes
# *******************************************************************

# Total de capturas e recaps
length(CRC)
table (Sexo)

summary(as.factor(Parcela)) # two observations without "parcela"]

#Capturas e recaps totais por Parcela e especie
(recap.especies<-table(Parcela,Especie))

# Capturas e recaps totais por Parcela
(cap<-table(Parcela))
(cap<-table(Parcela[Recaptura=="N"]))
(recap<-table(Parcela))

# Total de capturas (sem recaps)
sum(cap)
sum(recap)


# Media de recaps por individuo
sum(table(Parcela[Recaptura!="N"]))/length(CRC)
max(table(Numerodetombo[Recaptura=="S"]))
table((table(dados$Numerodetombo[dados$Recaptura=="S"])))


table(Parcela,Recaptura)
a<-table(Parcela,Recaptura)
colnames(a)
colnames(a)<-c("Capturas","Recapturas")
row.names(a) <- c("A2","A1","A3","A4")#renomeia parcelas
a <- (a[order(row.names(a)),])#Reordena as linhas pelas parcelas

a
windows(10,10)
barplot(t(a),beside = TRUE, legend.text=T,args.legend = list(x = "topleft"),
        main=expression(italic("Ameivula jalapensis")),
        ylim=c(0,120),ylab=expression(bold("Frequency")))

#Parcela e campanha
table(Campanha,Parcela)
b<-table(Campanha, Parcela)
colnames(b)<- c("A2","A1","A3","A4")#renomeia parcelas
b <- (b[,order(colnames(b))])#Reordena as linhas pelas parcelas
b

windows(10,10)
barplot(t(b),beside = TRUE, legend.text=T,args.legend = list(x = "topleft"),
        main=expression(italic("Ameivula jalapensis")),
        ylim=c(0,50),ylab=expression(bold("Abundância")))

table(Sexo,Parcela)
sexo<-table(Sexo,Parcela)
colnames(sexo)

colnames(sexo)<-c("A2","A1","A3","A4")
sexo <- (sexo[,order(colnames(sexo))])#Reordena as linhas pelas parcelas
sexo <- sexo[c(2,1,3),]

sexo

windows(20,10)
barplot(sexo,legend.text=T,args.legend = list(x = "topleft"),
        ylim=c(0,60),ylab="Frequ?ncia",cex.axis=1.0,
        cex.lab=1.0,cex.names=1.0,cex.sub=1,beside=TRUE)

(chi.sex.parcela <- chisq.test(sexo[-1,]))
chi.sex.parcela$residuals#A1 tem maiores residuos
(chi.sex.parcela <- chisq.test(sexo[-1,-1]))
prop.table(sexo[-1,],margin = 2)

(sexo.camp.parcela <- table(Sexo, Parcela, Campanha))
(chi.sex.camp1.parcela <- chisq.test(sexo.camp.parcela[-2,,1]))
(chi.sex.camp2.parcela <- chisq.test(sexo.camp.parcela[-2,,2]))
(chi.sex.camp3.parcela <- chisq.test(sexo.camp.parcela[-2,,3]))
(chi.sex.camp4.parcela <- chisq.test(sexo.camp.parcela[-2,,4]))

prop.table(sexo.camp.parcela[-2,,],margin = c(2,3))

chi.sex.camp3.parcela$residuals

# Diferencas no total de capturas e recaps entre tratamentos

recap
(chi.recap<-chisq.test(recap))
chi.recap$residuals # QT tem maior residuo
(chi.recap<-chisq.test(recap[-4]))
chi.recap$residual#QP tem maior residuo
(chi.recap<-chisq.test(recap[-c(3,4)]))#QP e QT possuem mais caps e recaps do que C1 e C2
(chi.recap<-chisq.test(recap[c(3,4)]))#Nao ha diferenca nas caps e recaps entre QP e QT


#N de individuos = caps
cap
(chi.cap<-chisq.test(cap))
chi.cap$residuals # QT tem maior residuo
(chi.cap<-chisq.test(cap[-4]))
chi.cap$residuals # QP tem maior residuo
(chi.cap<-chisq.test(cap[-c(3,4)]))
(chi.cap<-chisq.test(cap[c(3,4)]))



# ## Graficos da variacao temporal do CRC em todas as campanhas
# # *********************************************************

head(dados)
summary(CRC)
sd(CRC, na.rm = T)

windows(10,10)
plot(jitter(Campanha), CRC, ylab="Comprimento rostro-cloacal (mm)", 
     axes=F, cex=1.5, bty="n", ylim=c(30, 70), xlab="Campanha",pch=16,
     col=rgb(0,0,0,.4))
axis(1, 1:4, labels=c("Chuva","Seca","Chuva", "Seca"))
axis(2)
abline(h=40, lty=3)



windows(20,10)
boxplot(CRC~Campanha, axes=F, ylab="Comprimento rostro-cloacal (mm)",
        ylim=c(30, 70),range=0)
axis(1, 1:4, labels=c("Chuva","Seca","Chuva","Seca"))
axis(2)
abline(h=40, lty=3)

windows(10,10)
ggplot(dados, 
       aes(x = as.factor(Campanha), y = CRC)) + 
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
  geom_hline(yintercept = 40,linetype="dashed")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Snout-vent length (mm)")


##Graficos de razao sexual
head(dados)

table(Sexo,Campanha)
sexo<-table(Sexo,Campanha)
colnames(sexo)

colnames(sexo)<-c("Chuva","Seca","Chuva","Seca")
(sexo <- sexo[c(2,1,3),])

windows(10,10)
barplot(sexo,legend.text=T,args.legend = list(x = "topright"), beside=TRUE,
        ylim=c(0,90),ylab="Frequ?ncia",cex.axis=1.0,
        cex.lab=1.0,cex.names=1.0,cex.sub=1)

table(Sexo)

(chi.sex<-chisq.test(table(Sexo)[-2]))
(chi.sex.camp<-chisq.test(sexo[-1,]))
chi.sex.camp$residuals#Ha diferenca na razao sexual entre campanhas
(chi.sex.camp<-chisq.test(sexo[-1,-3]))
prop.table(sexo[-1,],margin = 2)

#Na terceira campanha ha razao sexual diferente do esperado ao acaso

#CRC entre sexos
head(dados)
windows(10,10)
boxplot(CRC~Sexo,data=dados,
        main=expression(italic("Ameivula jalapensis")),
        ylab="Comprimento rostro-cloacal (mm)")

#Testa diferença no CRC entre sexos
shapiro.test(CRC[Sexo!="I"])
bartlett.test(CRC~Sexo,data=dados[Sexo!="I",])

wilcox.test(CRC~Sexo,data=dados[Sexo!="I",])#Nao Significativo

#CRC entre parcelas
windows(10,10)
ggplot(dados, 
       aes(x = as.factor(Parcela), y = CRC)) + 
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
  geom_hline(yintercept = 40,linetype="dashed")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Comprimento rostro-cloacal (mm)")


#Testa diferença no CRC entre parcelas

bartlett.test(CRC~Parcela,data=dados)

r.crc<-rank(CRC)
anova.crc<-aov(r.crc~Parcela,data=dados)
summary(anova.crc)#Significativo
TukeyHSD(anova.crc)

r.crc.sexo<-rank(CRC[Sexo!="I"])
anova.crc.sexo<-aov(r.crc.sexo~Sexo*Parcela,data=dados[Sexo!="I",])
summary(anova.crc.sexo)

windows(10,10)
boxplot(CRC~Sexo*Parcela,data=dados[Sexo!="I",],
        main=expression(italic("Ameivula jalapensis")),
        ylab="Comprimento rostro-cloacal (mm)")

detach(dados)


##########################################################################
## BaSTA
##########################################################################

## Install packages
## ****************


## Clean data
## **********
head(dados)
# detach(dados)

demografia <- dados[dados$Morto!="S", ] #remove dead animals
table(demografia$Campanha)
table(demografia$Numerodetombo)

completos <- complete.cases(demografia[, c("Numerodetombo","Campanha")]) #remove NAs
demografia.planilha <- droplevels(demografia[completos, ])
head(demografia.planilha)
tail(demografia.planilha)
table(demografia.planilha$Recaptura,demografia.planilha$Especie)
str(demografia.planilha)
table(demografia.planilha$Campanha)
table(demografia.planilha$Numerodetombo)

demo.Ameivula<-demografia.planilha[demografia.planilha$Especie=="Ameivula_jalapensis",]

## Prepare input file and run monthly data ("campanha")
## ****************************************************

# Create capture histories

Y.Ameivula <- CensusToCaptHist(ID=demo.Ameivula$Numerodetombo, d=demo.Ameivula$Campanha, timeInt="M")
str(Y.Ameivula)
head(Y.Ameivula)
Y.Ameivula
length(Y.Ameivula$ID)
ch.Ameivula<-data.frame(ch=rep(NA,length(Y.Ameivula$ID)))

for(i in 1:length(Y.Ameivula$ID)){
  ch.Ameivula[i,]<-c(paste0(Y.Ameivula[i,2:ncol(Y.Ameivula)],sep="",collapse = ""))
}
ch.Ameivula

########
#OPENCR#
########

pts.arms<-read.table("Pontos_Armadilhas.txt",h=T)
pts.arms.SGT<-pts.arms[pts.arms$local=="EESGT",]

#(arms.pts.SGT <- subset(arms.pts,local=="SGT"))
coordinates(pts.arms.SGT) <- c("long","lat")
proj4string(pts.arms.SGT) <- CRS("+proj=longlat +datum=WGS84")
pts.arms.SGT.utm<-spTransform(pts.arms.SGT, CRS = CRS("+init=epsg:32722"))

pts.arms.SGT<-pts.arms[pts.arms$local=="EESGT",]
pts.arms.SGT$X <- coordinates(pts.arms.SGT.utm)[,1]
pts.arms.SGT$Y <- coordinates(pts.arms.SGT.utm)[,2]

names(pts.arms.SGT) <- c("Parcela","Gride","Lat","Long","Local",
                         "X","Y")

Ajalapensis.planilha.XY<-left_join(dados,
                               pts.arms.SGT,
                               by = c("Parcela","Gride"))
head(Ajalapensis.planilha.XY)


capts.1<-data.frame(Session = Ajalapensis.planilha.XY$Campanha[Ajalapensis.planilha.XY$Campanha==1],
                  ID = Ajalapensis.planilha.XY$Numerodetombo[Ajalapensis.planilha.XY$Campanha==1],
                  occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==1])),
                  X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Campanha==1],
                  Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Campanha==1],
                  data = Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==1],
                  Parcela = Ajalapensis.planilha.XY$Parcela[Ajalapensis.planilha.XY$Campanha==1])

capts.1$occasion[4:nrow(capts.1)]<-capts.1$occasion[4:nrow(capts.1)]+1
capts.1$occasion[19:nrow(capts.1)]<-capts.1$occasion[19:nrow(capts.1)]+3
capts.1$occasion <- capts.1$occasion+2

capts.2<-data.frame(Session = Ajalapensis.planilha.XY$Campanha[Ajalapensis.planilha.XY$Campanha==2],
                    ID = Ajalapensis.planilha.XY$Numerodetombo[Ajalapensis.planilha.XY$Campanha==2],
                    occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==2])),
                    X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Campanha==2],
                    Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Campanha==2],
                    data = Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==2],
                    Parcela = Ajalapensis.planilha.XY$Parcela[Ajalapensis.planilha.XY$Campanha==2])

capts.2$occasion[32:nrow(capts.2)]<-capts.2$occasion[32:nrow(capts.2)]+1


capts.3<-data.frame(Session = Ajalapensis.planilha.XY$Campanha[Ajalapensis.planilha.XY$Campanha==3],
                    ID = Ajalapensis.planilha.XY$Numerodetombo[Ajalapensis.planilha.XY$Campanha==3],
                    occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==3])),
                    X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Campanha==3],
                    Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Campanha==3],
                    data = Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==3],
                    Parcela = Ajalapensis.planilha.XY$Parcela[Ajalapensis.planilha.XY$Campanha==3])

capts.3$occasion[65:nrow(capts.3)]<-capts.3$occasion[65:nrow(capts.3)]+1


capts.4<-data.frame(Session = Ajalapensis.planilha.XY$Campanha[Ajalapensis.planilha.XY$Campanha==4],
                    ID = Ajalapensis.planilha.XY$Numerodetombo[Ajalapensis.planilha.XY$Campanha==4],
                    occasion = as.integer(as.factor(Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==4])),
                    X = Ajalapensis.planilha.XY$X[Ajalapensis.planilha.XY$Campanha==4],
                    Y = Ajalapensis.planilha.XY$Y[Ajalapensis.planilha.XY$Campanha==4],
                    data = Ajalapensis.planilha.XY$Data[Ajalapensis.planilha.XY$Campanha==4],
                    Parcela = Ajalapensis.planilha.XY$Parcela[Ajalapensis.planilha.XY$Campanha==4])
capts <- rbind(capts.1[,-6],
               capts.2[,-6],
               capts.3[,-6],
               capts.4[,-6])
# crc = oreadicus.planilha.XY$crc,
# massa = oreadicus.planilha.XY$massa,
# sexo = as.factor(oreadicus.planilha.XY$sexo))

completos <- complete.cases(capts[, c("ID", "X","Y")]) #remove NAs
capts <- droplevels(capts[completos, ])

traps.xy<-data.frame(trapID = interaction(pts.arms.SGT$Parcela,pts.arms.SGT$Gride),
                     x = pts.arms.SGT$X,
                     y = pts.arms.SGT$Y)
traps.xy <- read.traps(data = traps.xy[,-1], detector = "multi")

# env.vars.month$t <- rep(1:46,25)

# library(reshape)
# tmed <- cast(env.vars.month, trap ~ t, value="tmed"); tmed <- as.matrix(tmed[,-1])
# tmin2m<- cast(env.vars.month, trap ~ t, value="tmin2m"); tmin2m <- as.matrix(tmin2m[,-1])
# tmedskin<- cast(env.vars.month, trap ~ t, value="tmedskin"); tmedskin <- as.matrix(tmedskin[,-1])
# rhmax<- cast(env.vars.month, trap ~ t, value="rhmax"); rhmax <- as.matrix(rhmax[,-1])
# solmax<- cast(env.vars.month, trap ~ t, value="solmax"); solmax <- as.matrix(solmax[,-1])
# precip<- cast(env.vars.month, trap ~ t, value="precip"); precip<- as.matrix(precip[,-1])
# Toreadicus_perf<- cast(env.vars.month, trap ~ t, value="Toreadicus_perf"); Toreadicus_perf <- as.matrix(Toreadicus_perf[,-1])
# Toreadicus_ha90<- cast(env.vars.month, trap ~ t, value="Toreadicus_ha90"); Toreadicus_ha90 <- as.matrix(Toreadicus_ha90[,-1])
# 
# covariates(traps.xy) <- data.frame(cbind(tmed, tmin2m, tmedskin, rhmax, solmax, precip, Toreadicus_perf, Toreadicus_ha90))
# 
# timevaryingcov(traps.xy) <- list(tmed = 1:46,
#                                  tmin2m = (47:92),
#                                  tmedskin = (93:138),
#                                  rhmax = (139:184),
#                                  solmax = (185:230),
#                                  precip = (231:276),
#                                  Toreadicus_perf = (277:322),
#                                  Toreadicus_ha90 = (323:368))
# timevaryingcov(traps.xy)

Ajalapensis.CHxy<-make.capthist(capts,traps.xy,fmt = "XY",
                               noccasions = 15,
                               covnames = c("Parcela"))

intervals(Ajalapensis.CHxy) <- c(112,245,84)/30 

write.capthist(Ajalapensis.CHxy)

# time.env.cov <- data.frame(tmedJS = colMeans(tmed),
#                            tmin2mJS = colMeans(tmin2m),
#                            tmedskinJS = colMeans(tmedskin),
#                            rhmaxJS = colMeans(rhmax),
#                            solmaxJS = colMeans(solmax),
#                            precipJS = colMeans(precip),
#                            Toreadicus_perfJS = colMeans(Toreadicus_perf),
#                            Toreadicus_ha90JS = colMeans(Toreadicus_ha90))
# 

mesh <- make.mask(traps.xy, spacing = 30, type = "trapbuffer", buffer = 150)

summary(mesh)
summary(Ajalapensis.CHxy)
m.array(Ajalapensis.CHxy)
JS.counts(Ajalapensis.CHxy)
JS.direct(Ajalapensis.CHxy)

summary(traps(Ajalapensis.CHxy))

quartz(h= 12,w = 12)
par(mfrow=c(2,2))
plot(mesh)
plot(Ajalapensis.CHxy$`1`, tracks=T,add=T)

plot(mesh)
plot(Ajalapensis.CHxy$`2`, tracks=T,add=T)

plot(mesh)
plot(Ajalapensis.CHxy$`3`, tracks=T,add=T)

plot(mesh)
plot(Ajalapensis.CHxy$`4`, tracks=T,add=T)

#phi = taxa de sobrevivência
#p = taxa de recaptura/captura
#f = taxa de recrutamento
#CJS non-spatial

#Nulo
fit0<-openCR.fit(Ajalapensis.CHxy,type = "CJS")
summary(fit0)

#Variando pelo tempo
fit1<-openCR.fit(Ajalapensis.CHxy,type = "CJS", model = list(phi ~ -1 + t, 
                                                             p ~ -1 + t))
summary(fit1)

#Variando por parcela
fit2<-openCR.fit(Ajalapensis.CHxy,type = "CJS", model = list(phi ~ -1 + Parcela, 
                                                             p ~ -1 + Parcela))
summary(fit2)
predict(fit2,all=T)

#Variando por pacela e pelo tempo
fit3<-openCR.fit(Ajalapensis.CHxy,type = "CJS", model = list(phi ~ -1 + Parcela + t, 
                                                             p ~ -1 + Parcela + t))
summary(fit3)
predict(fit3,all=T)
#Modelo nao convergiu

AIC(fit0,fit1,fit2,fit3)

#Variando por pacela e pelo tempo
fit3<-openCR.fit(Ajalapensis.CHxy,type = "CJS", 
                 model = list(phi~ -1 + Parcela * t, 
                              p~ -1 + Parcela * t),
                 method ="CG")
summary(fit3)
predict(fit3,all=T)
#Modelo nao convergiu

AIC(fit0,fit1,fit2,fit3)

(parcela.dat <- expand.grid(Parcela = factor(c('C1','C2','QP','QT'), levels = c('C1','C2','QP','QT')), 
                            t = factor(1:4)))

quartz(8,8)
plot(fit3, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o')
plot(fit3, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o')
plot(fit3, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o')
plot(fit3, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o')

#Pradel non-spatial
pradel.fit0<-openCR.fit(Ajalapensis.CHxy,type = "JSSAfCL")
summary(pradel.fit0)

pradel.fit1<-openCR.fit(Ajalapensis.CHxy,type = "JSSAfCL", model = list(phi~t, p~t, f~t),
                        method = "Nelder-Mead", details = list(control = list(maxit = 5000)),
                        start=pradel.fit0)
summary(pradel.fit1)
p.fit1<-predict(pradel.fit1,all=T)

quartz(8,8)
ggplot(p.fit1$phi, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Survival") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit1$p, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Capture") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit1$f, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Recruitment") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))


pradel.fit2<-openCR.fit(Ajalapensis.CHxy,type = "JSSAfCL", model = list(phi~Parcela, p~Parcela, f~Parcela))
summary(pradel.fit2)
p.fit2 <- predict(pradel.fit2,all=T)
p.fit2$phi$Parcela <- factor(p.fit2$phi$Parcela, 
                           levels = c("C2", "C1", "QP", "QT"))

p.fit2$p$Parcela <- factor(p.fit2$p$Parcela, 
                           levels = c("C2", "C1", "QP", "QT"))

p.fit2$f$Parcela <- factor(p.fit2$f$Parcela, 
                           levels = c("C2", "C1", "QP", "QT"))

#--- Steps 1, 2, 3: Extract, Subset, and Rename (No Changes Here) ---
all_coefs_df <- coef(pradel.fit2)
all_vcov <- vcov(pradel.fit2)

phi_names <- grep("^phi", rownames(all_coefs_df), value = TRUE)
phi_coefs <- all_coefs_df[phi_names, "beta"]
names(phi_coefs) <- phi_names
phi_vcov <- all_vcov[phi_names, phi_names]

names(phi_coefs)[names(phi_coefs) == "phi"] <- "(Intercept)"
names(phi_coefs) <- gsub("phi.Parcela", "Parcela", names(phi_coefs))

#--- Step 4: Create the Data for the Grid and Run qdrg() ---

# Define the levels of your factor
factor_levels <- c("C2", "C1", "QP", "QT")

# THE FIX: Create a data frame containing your factor
grid_data <- data.frame(Parcela = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_phi <- qdrg(
  formula = ~ Parcela,
  data = grid_data,          # <-- Add the data argument here
  coef = phi_coefs,
  vcov = phi_vcov
)

#--- Step 5: Calculate and View the Pairwise Contrasts ---
contrasts_phi <- pairs(ref_grid_phi)

print(contrasts_phi)

p_names <- grep("^p", rownames(all_coefs_df), value = TRUE)[1:4]
p_coefs <- all_coefs_df[p_names, "beta"]
names(p_coefs) <- p_names
p_vcov <- all_vcov[p_names, p_names]

names(p_coefs)[names(p_coefs) == "p"] <- "(Intercept)"
names(p_coefs) <- gsub("p.Parcela", "Parcela", names(p_coefs))

#--- Step 4: Create the Data for the Grid and Run qdrg() ---

# Define the levels of your factor
factor_levels <- c("C2", "C1", "QP", "QT")

# THE FIX: Create a data frame containing your factor
grid_data <- data.frame(Parcela = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_p <- qdrg(
  formula = ~ Parcela,
  data = grid_data,          # <-- Add the data argument here
  coef = p_coefs,
  vcov = p_vcov
)

#--- Step 5: Calculate and View the Pairwise Contrasts ---
contrasts_p <- pairs(ref_grid_p)

print(contrasts_p)

f_names <- grep("^f", rownames(all_coefs_df), value = TRUE)[1:4]
f_coefs <- all_coefs_df[f_names, "beta"]
names(f_coefs) <- f_names
f_vcov <- all_vcov[f_names, f_names]

names(f_coefs)[names(f_coefs) == "f"] <- "(Intercept)"
names(f_coefs) <- gsub("f.Parcela", "Parcela", names(f_coefs))

#--- Step 4: Create the Data for the Grid and Run qdrg() ---

# Define the levels of your factor
factor_levels <- c("C1", "C2", "QP", "QT")

# THE FIX: Create a data frame containing your factor
grid_data <- data.frame(Parcela = factor_levels)

# Create the reference grid using the new 'grid_data' object
ref_grid_f <- qdrg(
  formula = ~ Parcela,
  data = grid_data,          # <-- Add the data argument here
  coef = f_coefs,
  vcov = f_vcov
)

#--- Step 5: Calculate and View the Pairwise Contrasts ---
contrasts_f <- pairs(ref_grid_f)

print(contrasts_f)

quartz(8,8)
ggplot(p.fit2$phi[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl, color = Parcela), size=1, show.legend = F) +
  labs(x="Fire severity", y="Survival") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5)) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_phi)

quartz(8,8)
ggplot(p.fit2$p[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl, color = Parcela), size=1, show.legend = F) +
  labs(x="Fire severity", y="Capture probability") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5)) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_p)

quartz(8,8)
ggplot(p.fit2$f[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl, color = Parcela), size=1, show.legend = F) +
  labs(x="Fire severity", y="Recruitment") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5)) +
  theme_minimal() +
  scale_color_manual(values = turbo(4))

print(contrasts_f)

pradel.fit3<-openCR.fit(Ajalapensis.CHxy,type = "JSSAfCL", 
                        model = list(phi~ -1 + Parcela*t, 
                                     p~ -1 + Parcela*t, 
                                     f~ -1 + Parcela*t),
                        start=pradel.fit0)
summary(pradel.fit3)
predict(pradel.fit3,all=T)

AIC(pradel.fit0,pradel.fit1,pradel.fit2,pradel.fit3)


pradel.fit3<-openCR.fit(Ajalapensis.CHxy,type = "JSSAfCL", 
                        model = list(phi~ -1 + Parcela*t, 
                                     p~ -1 + Parcela*t, 
                                     f~ -1 + Parcela*t),
                        method = "Nelder-Mead",
                        details = list(control = list(maxit=20000,
                                                      reltol = 1e-6)),
                        start = pradel.fit0)
summary(pradel.fit3)
predict(pradel.fit3,all=T)

AIC(pradel.fit0,pradel.fit1,pradel.fit2,pradel.fit3)

quartz(8,8)
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o')
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o')
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o')
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o')

quartz(8,8)
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o',par="f",ylim=c(0,2))
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o',par="f")
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o',par="f")
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o',par="f")

quartz(8,8)
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o',par="p",ylim=c(0,.05))
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o',par="p")
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o',par="p")
plot(pradel.fit3, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o',par="p")

d.pradel3 <- derived(pradel.fit3,all=T)
print(d.pradel3,legend=T)

quartz(8,8)
ggplot(d.pradel3[[1]]$estimates, aes(x=session,y=N,group=1)) +
  geom_point(aes(x=session,y=N), size=3,color="forestgreen") +
  geom_line(aes(x=session,y=N), size=1,color="forestgreen") +
  geom_point(data=d.pradel3[[2]]$estimates,aes(x=session,y=N), size=3,color="navyblue") +
  geom_line(data=d.pradel3[[2]]$estimates,aes(x=session,y=N), size=1,color="navyblue") +
  geom_point(data=d.pradel3[[3]]$estimates,aes(x=session,y=N), size=3,color="orange") +
  geom_line(data=d.pradel3[[3]]$estimates,aes(x=session,y=N), size=1,color="orange") +
  geom_point(data=d.pradel3[[4]]$estimates,aes(x=session,y=N), size=3,color="red") +
  geom_line(data=d.pradel3[[4]]$estimates,aes(x=session,y=N), size=1,color="red") +
  labs(x="Campanha", y="Tamanho populacional") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

ggplot(d.pradel3[[1]]$estimates, aes(x=session,y=phi,group=1)) +
  geom_point(aes(x=session,y=phi), size=3,color="forestgreen") +
  geom_line(aes(x=session,y=phi), size=1,color="forestgreen") +
  geom_point(data=d.pradel3[[2]]$estimates,aes(x=session,y=phi), size=3,color="navyblue") +
  geom_line(data=d.pradel3[[2]]$estimates,aes(x=session,y=phi), size=1,color="navyblue") +
  geom_point(data=d.pradel3[[3]]$estimates,aes(x=session,y=phi), size=3,color="orange") +
  geom_line(data=d.pradel3[[3]]$estimates,aes(x=session,y=phi), size=1,color="orange") +
  geom_point(data=d.pradel3[[4]]$estimates,aes(x=session,y=phi), size=3,color="red") +
  geom_line(data=d.pradel3[[4]]$estimates,aes(x=session,y=phi), size=1,color="red") +
  labs(x="Campaign", y="Survival") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

ggplot(d.pradel3[[1]]$estimates, aes(x=session,y=p,group=1)) +
  geom_point(aes(x=session,y=p), size=3,color="forestgreen") +
  geom_line(aes(x=session,y=p), size=1,color="forestgreen") +
  geom_point(data=d.pradel3[[2]]$estimates,aes(x=session,y=p), size=3,color="navyblue") +
  geom_line(data=d.pradel3[[2]]$estimates,aes(x=session,y=p), size=1,color="navyblue") +
  geom_point(data=d.pradel3[[3]]$estimates,aes(x=session,y=p), size=3,color="orange") +
  geom_line(data=d.pradel3[[3]]$estimates,aes(x=session,y=p), size=1,color="orange") +
  geom_point(data=d.pradel3[[4]]$estimates,aes(x=session,y=p), size=3,color="red") +
  geom_line(data=d.pradel3[[4]]$estimates,aes(x=session,y=p), size=1,color="red") +
  labs(x="Campaign", y="Capture probability") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

ggplot(d.pradel3[[1]]$estimates, aes(x=session,y=f,group=1)) +
  geom_point(aes(x=session,y=f), size=3,color="forestgreen") +
  geom_line(aes(x=session,y=f), size=1,color="forestgreen") +
  geom_point(data=d.pradel3[[2]]$estimates,aes(x=session,y=f), size=3,color="navyblue") +
  geom_line(data=d.pradel3[[2]]$estimates,aes(x=session,y=f), size=1,color="navyblue") +
  geom_point(data=d.pradel3[[3]]$estimates,aes(x=session,y=f), size=3,color="orange") +
  geom_line(data=d.pradel3[[3]]$estimates,aes(x=session,y=f), size=1,color="orange") +
  geom_point(data=d.pradel3[[4]]$estimates,aes(x=session,y=f), size=3,color="red") +
  geom_line(data=d.pradel3[[4]]$estimates,aes(x=session,y=f), size=1,color="red") +
  labs(x="Campaign", y="Recruitment") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

#Pradel SCR
pradel.scr.fit.0 <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                  model = list(phi ~ 1,
                               sigma ~ 1,
                               f ~ 1),
                  kernelradius = 5,
                  method = "Nelder-Mead",
                  link = list(phi = "sin"),
                  details = list(control = list(maxit=5000)),
                  start = pradel.fit0,
                  trace = 2,
                  ncores = 6)

summary(pradel.scr.fit.0)
predict(pradel.scr.fit.0,all=T)

pradel.scr.fit.BVN0 <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                               movementmodel = "BVN",
                               model = list(phi ~ 1,
                                            sigma ~ 1,
                                            f ~ 1),
                               kernelradius = 5,
                               method = "Nelder-Mead",
                               link = list(phi = "sin"),
                               details = list(control = list(maxit=5000)),
                               start = pradel.fit0,
                               trace = 2,
                               ncores = 6)

summary(pradel.scr.fit.BVN0)
predict(pradel.scr.fit.BVN0,all=T)


pradel.scr.fit.BVNzi0 <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                                    movementmodel = "BVNzi",
                                    model = list(phi ~ 1,
                                                 sigma ~ 1,
                                                 f ~ 1),
                                    kernelradius = 5,
                                    method = "Nelder-Mead",
                                    link = list(phi = "sin"),
                                    details = list(control = list(maxit=5000)),
                                    start = pradel.fit0,
                                    trace = 2,
                                    ncores = 6)

summary(pradel.scr.fit.BVNzi0)
predict(pradel.scr.fit.BVNzi0,all=T)

pradel.scr.fit.UNIzi0 <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                                    movementmodel = "UNIzi",
                                    model = list(phi ~ 1,
                                                 sigma ~ 1,
                                                 f ~ 1),
                                    kernelradius = 10,
                                    method = "Nelder-Mead",
                                    link = list(phi = "sin", move.a = "sin"),
                                    details = list(control = list(maxit=5000)),
                                    start = pradel.fit0,
                                    trace = 2,
                                    ncores = 6)

summary(pradel.scr.fit.UNIzi0)
predict(pradel.scr.fit.UNIzi0,all=T)

AIC(pradel.scr.fit.0, pradel.scr.fit.BVN0, pradel.scr.fit.BVNzi0, pradel.scr.fit.UNIzi0)

pradel.scr.fit.t <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", 
                             mask = mesh,
                             movementmodel = "UNIzi",
                             kernelradius = 10,
                  model = list(phi ~ -1 + t,
                               lambda0 ~ -1 + t,
                               sigma ~ 1,
                               f ~ -1 + t,
                               move.a ~ -1 + t),
                  link = list(phi = "sin", move.a = "sin"),
                  method = "Nelder-Mead",
                  details = list(control = list(maxit=10000,
                                                reltol = 1e-6)),
                  start = list(pradel.scr.fit.UNIzi0, pradel.fit1),
                  trace = 2,
                  ncores = 6)

summary(pradel.scr.fit.t)

pradel.scr.fit.parcela <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                     model = list(phi ~ -1 + Parcela,
                                  lambda0 ~ -1 + Parcela,
                                  sigma ~ 1,
                                  f ~ -1 + Parcela,
                                  move.a ~ -1 + Parcela),
                     # link = list(phi = "sin", move.a = "sin"),
                     movementmodel = "UNIzi",
                     kernelradius = 10,
                     method = "BFGS",
                     details = list(control = list(maxit=10000,
                                                   reltol = 1e-6)),
                     start = list(pradel.scr.fit.UNIzi0, pradel.fit2),
                     trace = 2,
                     ncores = 6)

summary(pradel.scr.fit.parcela)
predict(pradel.scr.fit.parcela,all=T)

AIC(pradel.scr.fit.parcela,
    pradel.scr.fit.t, pradel.scr.fit.UNIzi0)

pradel.scr.fit.full <- openCR.fit(Ajalapensis.CHxy,type = "JSSAsecrfCL", mask = mesh,
                                     model = list(phi ~ -1 + Parcela * t,
                                                  lambda0 ~ -1 + Parcela * t,
                                                  #sigma ~ -1 +Parcela + t,
                                                  f ~ -1 + Parcela * t),
                                  movementmodel = "UNIzi",
                                  kernelradius = 5,
                                  method = "Nelder-Mead",
                                  details = list(control = list(maxit=10000)),
                                  link = list(phi = "sin", move.b = "sin"),
                                  start = list(pradel.scr.fit.UNIzi0, pradel.fit3),
                                  ncores = 6)

summary(pradel.scr.fit.full)
predict(pradel.scr.fit.full,all=T)

AIC(pradel.scr.fit.full,pradel.scr.fit.parcela,
    pradel.scr.fit.t, pradel.scr.fit.0)

k0 <- make.kernel(pradel.scr.fit.BVNzi0, truncate = 150)
summary(k0)
plot(k0)
plot(k0, type = "Fr")

kC1 <- make.kernel("BVNzi", kernelradius = 5, spacing = 30, move.a = 150, move.b = 0.53, sparsekernel = T, truncate = 150)
summary(kC1)
plot(kC1)
plot(kC1, type = "Fr")

kBP <- make.kernel("BVNzi", kernelradius = 5, spacing = 30, move.a = 150, move.b = 0.7, sparsekernel = T)
summary(kBP)
plot(kBP, contour = T)
plot(kBP, type = "Fr")

kBT <- make.kernel("BVNzi", kernelradius = 5, spacing = 30, move.a = 150, move.b = 0.87, sparsekernel = T)
summary(kBT)
plot(kBT)
plot(kBT, type = "Fr")

ka <- make.kernel("BVN", kernelradius = 5, spacing = 30, move.a = 150, move.b = 0.53)
summary(ka)
plot(ka)
plot(ka, type = "gr")
plot(ka, type = "fr")
plot(ka, type = "Fr")

plot(c(0:150), pkernel(c(0:150), "BVNzi", move.a = 150, move.b = 0.53), 
     type = "l", bty = "n", ylim = c(0.5, 1), ylab = "Cumulative probability", xlab = "Distance (m)")
lines(c(0:150), pkernel(c(0:150), "BVNzi", move.a = 150, move.b = 0.7), type = "l", bty = "n", col = "orange")
lines(c(0:150), pkernel(c(0:150), "BVNzi", move.a = 150, move.b = 0.87), type = "l", bty = "n", col = "red")

plot(c(0:150), qkernel(c(0:150), "BVN", move.a = 150, move.b = 0.53), type = "l", bty = "n")


p.fit.scr.t<-predict(pradel.scr.fit.t,all=T)

quartz(8,8)
ggplot(p.fit.scr.t$phi, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Survival") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit.scr$lambda0, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Capture") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit.scr$f, aes()) +
  geom_pointrange(aes(x=session,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Occasion", y="Recruitment") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

p.fit.scr.parcela<-predict(pradel.scr.fit.parcela,all=T)

quartz(8,8)
ggplot(p.fit.scr.parcela$phi[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Plot", y="Survival") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit.scr.parcela$lambda0[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Plot", y="Capture") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit.scr.parcela$f[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=estimate,ymin=lcl,ymax=ucl), size=1) +
  labs(x="Plot", y="Recruitment") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))

quartz(8,8)
ggplot(p.fit.scr.parcela$move.a[c(1,5,9,13),], aes()) +
  geom_pointrange(aes(x=Parcela,y=1- estimate,ymin=1 - lcl,ymax=1- ucl), size=1) +
  labs(x="Plot", y="Dispersal probability") +
  theme(plot.margin=unit(c(0, 0, 2, 0), "lines"), text=element_text(size=10), legend.position=c(0.1, 0.5))



plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o',par='phi')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o',par='phi')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o',par='phi')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o',par='phi')


plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o',par='f',ylim=c(0,2))
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o',par='f')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o',par='f')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o',par='f')

plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C1',], #add = TRUE, 
     xoffset = 0.1, col = 'forestgreen', pch = 16,type='o',par='lambda0',ylim=c(0,.01))
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='C2',], add = TRUE, 
     xoffset = 0.1, col = 'navyblue', pch = 16,type='o',par='lambda0')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QP',], add = TRUE, 
     xoffset = 0.1, col = 'orange', pch = 16,type='o',par='lambda0')
plot(pradel.scr.fit.full, newdata = parcela.dat[parcela.dat$Parcela=='QT',], add = TRUE, 
     xoffset = 0.1, col = 'red', pch = 16,type='o',par='lambda0')

d.pradel.scr.full <- derived(pradel.scr.fit.full,all=T)
print(d.pradel.scr.full)

saveRDS(fit.t, "best.scr.rds")


ggplot(p.predict.maximiliani$p,aes(time,estimate,ymin=lcl,ymax=ucl))+
  geom_errorbar(width=0.2)+geom_point()+geom_line()+
  xlab("\nMeses")+ylab("Recaptura\n")


# Bayesian CJS/Pradel/Growth Model ---------------------------------------------------

library(runjags)
library(parallel)
library(rjags)
library(reshape)
library(lubridate)

demo.Ameivula$Data <- as.Date(demo.Ameivula$Data, format = "%d/%m/%Y")

demo.Ameivula$TimeWeek <- 1 + lubridate::interval("2021-02-22", 
                                                      demo.Ameivula$Data) %/% 
  weeks(1)

demo.Ameivula$TimeMonth <- 1 + lubridate::interval("2021-02-22", 
                                                  demo.Ameivula$Data) %/% 
  months(1)

demo.Ameivula$TimeMonth[demo.Ameivula$TimeMonth == 14] <- 13
demo.Ameivula$TimeMonth[demo.Ameivula$TimeMonth == 17] <- 16
demo.Ameivula$TimeMonth

## Prepare input file and run monthly data ("TimeMonth")
###################################################

# Create ID
IDENT <- paste0(demo.Ameivula$Parcela, demo.Ameivula$Numerodetombo)
head(IDENT)
demo.Ameivula <- data.frame(IDENT,demo.Ameivula)
rm(IDENT)
str(demo.Ameivula)

# Filter the variables of interest
table1 <- demo.Ameivula[, c("IDENT", "TimeWeek", "Sexo", "CRC", "Parcela")]
str(table1)


# Identifies captures without ID and SVL
table2 <- complete.cases(table1[, c("IDENT", "CRC")])

# Removes captures without ID and SVL
table3 <- table1[table2, ]

# Order the data by ID and TimeMonth
table4 <- table3[order(table3$IDENT, table3$TimeWeek), ]
str(table4)
summary(table4)

# Drop unused levels
table5 <- droplevels(table4)
str(table5)


# Calculates recapture frequencies
recap.table <- data.frame(table(table5$IDENT))
names(recap.table) <- c("ID", "captures")
recap.table
table(recap.table$captures)
table(recap.table$captures) * as.numeric(names(table(recap.table$captures)))
sum(table(recap.table$captures) * as.numeric(names(table(recap.table$captures))))

#Sanity check
sum(table(recap.table$captures) * as.numeric(names(table(recap.table$captures))))==nrow(table5)

## Filters dataset to use in the analyses##
###########################################
head(table5)

Age<-c(rep(NA,nrow(table5)))
Age

datA<-data.frame(table5$TimeWeek,table5$Sexo,table5$IDENT,table5$CRC,Age,table5$Parcela)
names(datA)<-c("TimeWeek","Sex","TrueID","SVL","Age","Site")
datA$Age[datA$SVL<=35]<-0
# datA<-datA[datA$TrueID!="Q3021",]
# datA<-datA[datA$TrueID!="BP1191",]
# datA<-datA[datA$TrueID!="BP2101",]
# datA<-datA[datA$TrueID!="BP3511",]
head(datA)
tail(datA)
str(datA)


##Creates variables for growth model
####################################
###  del is the time period since first individual's capture (0 in first capture)

del<-c()   ### time since first capture

for(i in 1:nrow(datA)){
  del[i] <- datA$TimeWeek[i] - min(datA$TimeWeek[datA$TrueID==datA$TrueID[i]])
}

plot<-cast(datA, TrueID~., value="Site", fun.aggregate=function(x) tail(x,1))  ###determine the site from each individual - filters the last value
plot<-as.character(plot[,2])
plot
#Ordering by fire severity
plot[plot=="C2"] <- 1
plot[plot=="C1"]<-2
plot[plot=="QP"]<-3
plot[plot=="QT"]<-4

plot<-as.numeric(plot)
plot
ind = as.numeric(factor(datA$TrueID))
y = datA$SVL
(n = max(ind) ) ### number of individuals
(m = nrow(datA))### number of observations (captures)

age<- c()  ## idade na primeira captura
for (a in 1:n){ age[a] <- datA$Age[ind==a][1]}

time <- c()
for (a in 1:n){ time[a] <- datA$TimeWeek[ind==a][1]}

head(datA)
tail(datA)
str(datA)

##### #Cria dados de recupera??o de marca para o modelo de sobreviv?ncia## ##################################################################
known.states.cjs<-function(ch){
  state<-ch
  for (i in 1:dim(ch)[1]){
    n1<-min(which(ch[i,]==1))
    n2<-max(which(ch[i,]==1))
    state[i,n1:n2]<-1
    state[i,n1]<-NA
  }
  state[state==0]<-NA
  return(state)
}

cjs.init.z<-function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2<-max(which(ch[i,]==1))
    ch[i,f[i]:n2]<-NA
  }
  for (i in 1:dim(ch)[1])
  { ch[i,1:f[i]]<-NA
  }
  return(ch)
}

eh <- cast(datA,TrueID ~ TimeWeek, fun.aggregate = function(x) as.numeric(length(x) >0),value="SVL");eh <- eh[,2:ncol(eh)]
eh.all <- seq(min(datA$TimeWeek), max(datA$TimeWeek)) #fill ignore time periods
missing <- eh.all[!(eh.all %in% names(eh))]
col=matrix(0,nrow=nrow(eh),ncol=length(missing))
colnames(col) <- missing
eh <- cbind(eh, col)
eh <- eh[,order(as.integer(colnames(eh)))]

#Sanity check
which(colSums(eh)==0) == missing

head(eh)
(f <- apply(eh,1,function(x) which(x==1)[1]))#Time of first capture (individual)
(nind <- nrow(eh))
(n.occasions <- ncol(eh))
m
n

mplot <- data.frame(C = as.numeric(plot==1),
                    Q = as.numeric(plot==2),
                    BP = as.numeric(plot==3),
                    BM = as.numeric(plot==4),
                    BT = as.numeric(plot==5))

# Create matrix X indicating crc
x <- reshape2::acast(datA, TrueID ~ TimeWeek, fun.aggregate = function(x) mean(x, na.rm = T), 
                     value.var = "SVL")

x[is.nan(x)==T] <- NA

x.all <- seq(min(datA$TimeWeek), max(datA$TimeWeek)) #fill ignore time periods
missing.x <- x.all[!(x.all %in% colnames(x))]
missing.x == missing#Sanity check
col=matrix(NA,nrow=nrow(x),ncol=length(missing))
colnames(col) <- missing
x <- cbind(x, col)
x <- x[,order(as.integer(colnames(x)))]

head(x)

ncol(x) == ncol(eh)

missing.cjs<-function(ch, plot){
  state<-ch
  state[state==0]<-NA
  for(s in 1:4){
    state[plot == s, which(colSums(eh[plot == s,]) < 1)] <- 0
  }
  state[is.na(state)] <- 1
  return(state)
}

missing.samp <- missing.cjs(eh, plot)#1 for months samples, 0 for months not sampled
View(missing.samp)

hist(datA$SVL[datA$SVL<=35])
mean(datA$SVL[datA$SVL<=35],na.rm=T)
var(datA$SVL[datA$SVL<=35],na.rm=T)

hist(datA$SVL)
hist(datA$SVL[datA$SVL>55])

mean(datA$SVL[datA$SVL>55],na.rm=T)
var(datA$SVL[datA$SVL>55],na.rm=T)

hist(dados$CRC[dados$Recaptura=="N"])
hist(dados$CRC[dados$Recaptura=="S"])

bugs.data <- list(first = f, nind = dim(eh)[1], n.occasions = dim(eh)[2],
                  y = eh, x = as.matrix(x), z = known.states.cjs(eh), 
                  missing = missing.samp,
                  mu.L0 = mean(datA$SVL[datA$SVL<=35],na.rm=T),
                  tau.L0 = var(datA$SVL[datA$SVL<=35],na.rm=T),
                  AFC = as.numeric(age),
                  #mplot = mplot,
                  plot = plot)
# Valores iniciais
inits <- function (){}

# Defina os parametros a serem monitorados
parameters <- c("alpha.phi","beta.phi", "beta2.phi",
                "alpha.p","beta.p", "beta2.p",
                "p.AFC","r.AFC","var.AFC","mn.AFC",
                "mu.K","mu.LI"
)

# Specify model in BUGS language
sink("cjs-Ajalapensis-crc.jags")
cat("

model {

  #########################
  ## SURVIVAL/GROWTH MODEL#
  #########################

  
   for(i in 1:nind){
    AFC[i] ~ dnegbin(p.AFC , r.AFC)T(0,50)  ###Change trucationfor different species. AFC is age of first capture - known in many cases (for newborns) and estimated when not known
    L0[i] ~ dnorm(mu.L0, tau.L0)  ### draw values for intial size
    LI[i] ~ dnorm(mu.LI,tau.LI)T(0,) ### asymptotic size - taubeta allows for individual variation, while mean size is plot dependent
    newLI[i] ~ dnorm(mu.LI,tau.LI)T(0,)
    LLoldLI[i] <-logdensity.norm(LI[i],mu.LI,tau.LI)
    LLnewLI[i] <-logdensity.norm(newLI[i],mu.LI,tau.LI)
    logit(K[i]) <- K.L[i]  ## mean growth rate is plot dependent with variation defined by xi*theta
    K.L[i] ~ dnorm(mu.K[plot[i]],tau.K)
  }
  
  ### Priors for the growth model
  for(i in 1:4){
    mu.K1[i] ~ dunif(0.5,1)
    mu.K[i] <- log(mu.K1[i]) - log(1-mu.K1[i])
  }
  r.AFC ~ dgamma(0.01,0.01)
  p.AFC <- r.AFC/(r.AFC+mn.AFC)
  mn.AFC ~ dgamma(0.01,0.01)
  var.AFC <- r.AFC*(1-p.AFC)/(p.AFC*p.AFC)
  sd.sample ~ dt(0,0.0004,3)T(0,) ### t priors as in Schofield et al. 2013
  mu.LI ~ dnorm(57,1) #Change for different species
  sd.LI ~ dt(0,0.0004,3)T(0,)
  sd.K ~ dt(0,0.0004,3)T(0,)
  tau.sample <- 1/(sd.sample^2)
  tau.LI <- 1/(sd.LI^2)
  tau.K <- 1/(sd.K^2)    

  
# Priors and constraints
for (i in 1:nind){
  for (t in first[i]:(n.occasions-1)){
  logit(phi[i,t]) <- alpha.phi + beta.phi*x[i,t] + beta2.phi*pow(x[i,t],2)
  logit(p[i,t]) <-  (alpha.p + beta.p*x[i,t] + beta2.p*pow(x[i,t],2)) * missing[i,t]
  } #t
} #i

#### PRIORS
alpha.phi ~ dnorm(0, 0.25) # Prior for survival intercept
beta.phi ~ dnorm(0, 0.25) # Prior for slope parameter
beta2.phi ~ dnorm(0, 0.25) #Prior for quadratic slope parameter
alpha.p ~ dnorm(0, 0.25) #Prior for capture intercept
beta.p ~ dnorm(0, 0.25) # Prior for slope parameter
beta2.p ~ dnorm(0, 0.25) # Prior for quadratic slope parameter



# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,first[i]] <- 1
   for (t in (first[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      x[i,t-1] ~ dnorm(L[i,t-1],tau.sample)T(0,)
      L[i,t-1] <- (L0[i] + (LI[i]-L0[i])*(1-K[i]^(AFC[i]+(t-1)-first[i]))) 
      } #t
   } #i
   


}
",fill = TRUE)
sink()

# MCMC settings
ni <- 50000
nt <- 1
nb <- 90000
nc <- 4
na <- 10000

#Call JAGS from R (BRT 3 min)
bugs.data$y <- as.matrix(bugs.data$y)
bugs.data$z <- as.matrix(bugs.data$z)
bugs.data$missing <- as.matrix(bugs.data$missing)

runjags.options(jagspath = "/usr/local/bin/jags")
# setwd("C:/Users/UFT/Desktop/Ameivula_jalapensis_2022_2023")

cjs.Ajalapensis <- run.jags(data=bugs.data, inits=inits, monitor=parameters, 
                            model="cjs-Ajalapensis-crc.jags",
                            n.chains = nc, adapt = na,thin = nt, sample = ni,
                            burnin = nb,
                        method = "bgparallel", jags.refresh = 10,
                        keep.jags.files = TRUE,
                        summarise = TRUE,
                        modules = c("glm"))

# chain1 <- read.table("C:/Users/UFT/Desktop/Ameivula_jalapensis_2022_2023/runjagsfiles_1/sim.1/CODAchain1.txt")
# chain2 <- read.table("C:/Users/UFT/Desktop/Ameivula_jalapensis_2022_2023/runjagsfiles_1/sim.1/CODAchain2.txt")
# CODAindex <- read.table("C:/Users/UFT/Desktop/Ameivula_jalapensis_2022_2023/runjagsfiles_1/sim.1/CODAindex.txt")
# 
# chain1$par <- NA
# chain2$par <- NA
# for(i in 1:nrow(CODAindex)){
#   chain1$par[CODAindex$V2[i]:CODAindex$V3[i]] <- CODAindex$V1[i]
# }
# 
# for(i in 1:nrow(CODAindex)){
#   chain2$par[CODAindex$V2[i]:CODAindex$V3[i]] <- CODAindex$V1[i]
# }
# 
# names(chain1) <- c("niter", "value", "par")
# names(chain2) <- c("niter", "value", "par")
# 
# m.chain1 <- pivot_wider(data = chain1,
#             names_from = "par",
#             values_from = "value")
# m.chain2 <- pivot_wider(data = chain2,
#                         names_from = "par",
#                         values_from = "value")
# 
# results.cjs.Ajalapensis <- as.mcmc.list(list(as.mcmc(m.chain1[,-1]), as.mcmc(m.chain2[,-1])))

library(MCMCvis)

results.cjs.Ajalapensis <- results.jags(cjs.Ajalapensis)#runjagsfiles

results.cjs.Ajalapensis.extend1 <- extend.jags(results.cjs.Ajalapensis, 
                                               adapt = 10000, burnin = 40000, sample = 50000, 
                                               combine = F, keep.jags.files = TRUE, 
                                               method = "bgparallel")

results.cjs.Ajalapensis <- results.jags(results.cjs.Ajalapensis.extend1)
 #results.cjs.Ajalapensis<- add.summary(results.cjs.Ajalapensis)

results.cjs.Ajalapensis.extend2 <- extend.jags(results.cjs.Ajalapensis, 
                                               adapt = 10000, sample = 100000, 
                                               combine = F, keep.jags.files = TRUE, 
                                               method = "bgparallel")

results.cjs.Ajalapensis <- results.jags(results.cjs.Ajalapensis.extend2)

results.cjs.Ajalapensis.extend3 <- extend.jags(results.cjs.Ajalapensis, 
                                               adapt = 10000, sample = 100000, 
                                               combine = T, keep.jags.files = TRUE, 
                                               method = "bgparallel")

results.cjs.Ajalapensis <- results.jags(results.cjs.Ajalapensis.extend3)

results.cjs.Ajalapensis.extend4 <- extend.jags(results.cjs.Ajalapensis, 
                                               adapt = 10000, burnin = 140000, sample = 300000, 
                                               combine = F, keep.jags.files = TRUE, 
                                               method = "bgparallel")

results.cjs.Ajalapensis <- results.jags(results.cjs.Ajalapensis.extend4)

results.cjs.Ajalapensis.extend5 <- extend.jags(results.cjs.Ajalapensis, 
                                               adapt = 10000,  sample = 200000, 
                                               combine = T, keep.jags.files = TRUE, 
                                               method = "bgparallel")

results.cjs.Ajalapensis <- results.jags(results.cjs.Ajalapensis.extend5)



results.cjs.Ajalapensis.df <- summary(results.cjs.Ajalapensis)
View(results.cjs.Ajalapensis.df)

write.csv(results.cjs.Ajalapensis.df, "results_cjs_Ajalapensis_df.csv")

saveRDS(results.cjs.Ajalapensis, "results_cjs_Ajalapensis.rds")

library(ggmcmc)

S <- ggs(results.cjs.Ajalapensis$mcmc)

ggs_density(S, family = "alpha.p")
ggs_density(S, family = "beta.p")
ggs_density(S, family = "beta2.p")

ggs_traceplot(S, family = "alpha.p")
ggs_traceplot(S, family = "beta.p")
ggs_traceplot(S, family = "beta2.p")

#Plots

#Function to estimate size from age
age_to_size <- function(x,mu.L0,mu.LI,K) mu.L0 + (mu.LI-mu.L0)*(1-plogis(K)^x)

#Function to estimate age from size
size_to_age <- function(x,mu.L0,mu.LI,K) log(1-((x - mu.L0)/(mu.LI - mu.L0)))/log(inv_logit(K))

#Function to estimate size in t1 from size in t0
sizet0_t1 <- function(x,mu.L0,mu.LI,K) age_to_size(size_to_age(x,mu.L0,mu.LI,K)+1,mu.L0,mu.LI,K)

results.cjs.Ajalapensis.df <- read.csv("results_cjs_Ajalapensis_df.csv")

(xx <- 0:12)

(mu.L0 <- min(dados$CRC, na.rm = T))
(mu.LI <- results.cjs.Ajalapensis.df$Mean[grep(pattern = "mu.LI", 
                                                  x = results.cjs.Ajalapensis.df$X)])

(K <- results.cjs.Ajalapensis.df$Mean[grep(pattern = "mu.K", 
                                             x = results.cjs.Ajalapensis.df$X)])

(K.lw <- results.cjs.Ajalapensis.df$Lower95[grep(pattern = "mu.K", 
                                                x = results.cjs.Ajalapensis.df$X)])

(K.up <- results.cjs.Ajalapensis.df$Upper95[grep(pattern = "mu.K", 
                                           x = results.cjs.Ajalapensis.df$X)])


    
curve(mu.L0 + (mu.LI-mu.L0)*(1-plogis(K[1])^x), xlim = c(0, 12), 
      ylab = "Snout-vent length (mm)", xlab = "Time (weeks)", bty = "n", col = turbo(4)[1], lwd = 1.2)
curve(mu.L0 + (mu.LI-mu.L0)*(1-plogis(K[2])^x), xlim = c(0, 48), col = turbo(4)[2], lwd = 1.2, add = T)
curve(mu.L0 + (mu.LI-mu.L0)*(1-plogis(K[3])^x), xlim = c(0, 48), col = turbo(4)[3], lwd = 1.2, add = T)
curve(mu.L0 + (mu.LI-mu.L0)*(1-plogis(K[4])^x), xlim = c(0, 48), col = turbo(4)[4], lwd = 1.2, add = T)

growth.df <- data.frame(plot = as.factor(rep(1:4, each = length(xx))),
                        mean = c(age_to_size(xx, mu.L0, mu.LI, K[1]),
                                 age_to_size(xx, mu.L0, mu.LI, K[2]),
                                 age_to_size(xx, mu.L0, mu.LI, K[3]),
                                 age_to_size(xx, mu.L0, mu.LI, K[4])),
                        lower = c(age_to_size(xx, mu.L0, mu.LI, K.lw[1]),
                                  age_to_size(xx, mu.L0, mu.LI, K.lw[2]),
                                  age_to_size(xx, mu.L0, mu.LI, K.lw[3]),
                                  age_to_size(xx, mu.L0, mu.LI, K.lw[4])),
                        upper = c(age_to_size(xx, mu.L0, mu.LI, K.up[1]),
                                  age_to_size(xx, mu.L0, mu.LI, K.up[2]),
                                  age_to_size(xx, mu.L0, mu.LI, K.up[3]),
                                  age_to_size(xx, mu.L0, mu.LI, K.up[4])),
                        time = rep(xx, 4)
                        )

quartz(width = 10, height = 8)
ggplot(data = growth.df, aes(x = time, y = mean, colour = plot))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = plot), alpha = 0.2, colour = NA) +
  labs(x = "Time (weeks)", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  facet_wrap("plot")

ggplot(data = growth.df, aes(x = time, y = mean, colour = plot))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = plot), alpha = 0.2, colour = NA) +
  labs(x = "Time (weeks)", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(4), name = "Fire severity") +
  scale_fill_manual(values = turbo(4), name = "Fire severity")


# Pradel model ------------------------------------------------------------
# Filter the variables of interest
table1 <- demo.Ameivula[, c("IDENT", "TimeMonth", "Sexo", "Parcela")]
str(table1)
summary(table1)

# Identifies captures without ID and SVL
table2 <- complete.cases(table1[, c("IDENT")])

# Removes captures without ID and SVL
table3 <- table1[table2, ]

# Order the data by ID and TimeMonth
table4 <- table3[order(table3$IDENT, table3$TimeMonth), ]
str(table4)
summary(table4)

# Drop unused levels
table5 <- droplevels(table4)
str(table5)


# Calculates recapture frequencies
recap.table <- data.frame(table(table5$IDENT))
names(recap.table) <- c("ID", "captures")
recap.table
table(recap.table$captures)
table(recap.table$captures) * as.numeric(names(table(recap.table$captures)))
sum(table(recap.table$captures) * as.numeric(names(table(recap.table$captures))))

#Sanity check
sum(table(recap.table$captures) * as.numeric(names(table(recap.table$captures))))==nrow(table5)

## Filters dataset to use in the analyses##
###########################################
head(table5)

datA<-data.frame(table5$TimeMonth, table5$Sexo, table5$IDENT, table5$Parcela)
names(datA)<-c("TimeMonth","Sex","TrueID","Site")

head(datA)
tail(datA)
str(datA)

eh <- cast(datA,TrueID ~ TimeMonth, fun.aggregate = function(x) as.numeric(length(x) >0),value="Site");eh <- eh[,2:ncol(eh)]
eh.all <- seq(min(datA$TimeMonth), max(datA$TimeMonth)) #fill ignore time periods
missing <- eh.all[!(eh.all %in% names(eh))]
col=matrix(0,nrow=nrow(eh),ncol=length(missing))
colnames(col) <- missing
eh <- cbind(eh, col)
eh <- eh[,order(as.integer(colnames(eh)))]

#Sanity check
which(colSums(eh)==0) == missing

head(eh)

plot<-cast(datA, TrueID~., value="Site", fun.aggregate=function(x) tail(x,1))  ###determine the site from each individual - filters the last value
plot<-as.character(plot[,2])
plot
#Ordering by fire severity
plot[plot=="C2"] <- 1
plot[plot=="C1"]<-2
plot[plot=="QP"]<-3
plot[plot=="QT"]<-4

plot<-as.numeric(plot)
plot
(f <- apply(eh,1,function(x) which(x==1)[1]))#Time of first capture (individual)
(nind <- nrow(eh))
(n.occasions <- ncol(eh))
ind = as.numeric(factor(datA$TrueID))

(n = max(ind) ) ### number of individuals
(m = nrow(datA))### number of observations (captures)

#Sanity checks

missing.cjs<-function(ch, plot){
  state<-ch
  state[state==0]<-NA
  for(s in 1:4){
    state[plot == s, which(colSums(eh[plot == s,]) < 1)] <- 0
  }
  state[is.na(state)] <- 1
  return(state)
}

missing.samp <- missing.cjs(eh, plot)#1 for months samples, 0 for months not sampled
View(missing.samp)

# Derivar dados para o modelo:
# e = indice da observacao mais antiga
get.first <- function (x) min(which (x!=0))
e <- apply (eh, 1, get.first)
e <- data.frame(plot=plot,e=e)


# l = indice da ultima observacao
get.last <- function(x) max(which(x!=0))
l <- apply (eh,1, get.last )
l <- data.frame(plot=plot,l=l)


# u = numero de animais observados pela primeira vez em i
u1 <- data.frame(table(e$e[e$plot==1]))
u2 <- data.frame(table(e$e[e$plot==2]))
u3 <- data.frame(table(e$e[e$plot==3]))
u4 <- data.frame(table(e$e[e$plot==4]))

u.all <- seq(1, n.occasions) # Preencher todos os anos ignorados
missing <- u.all[!(u.all %in% u1$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u1) <- names(df)
u1 <- rbind(u1, df)
u1$e<-as.numeric(levels(u1$e))[u1$e]
u1<-u1[order(u1$e), ]
u1<-u1$Freq
u1

missing <- u.all[!(u.all %in% u2$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u2) <- names(df)
u2 <- rbind(u2, df)
u2$e<-as.numeric(levels(u2$e))[u2$e]
u2<-u2[order(u2$e), ]
u2<-u2$Freq
u2

missing <- u.all[!(u.all %in% u3$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u3) <- names(df)
u3 <- rbind(u3, df)
u3$e<-as.numeric(levels(u3$e))[u3$e]
u3<-u3[order(u3$e), ]
u3<-u3$Freq
u3

missing <- u.all[!(u.all %in% u4$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u4) <- names(df)
u4 <- rbind(u4, df)
u4$e<-as.numeric(levels(u4$e))[u4$e]
u4<-u4[order(u4$e), ]
u4<-u4$Freq
u4

(u <- rbind(u1,u2,u3,u4))

# n = numero de animais observados em i
n1 <- colSums(eh[plot==1,])
n2 <- colSums(eh[plot==2,])
n3 <- colSums(eh[plot==3,])
n4 <- colSums(eh[plot==4,])

n<- rbind(n1,n2,n3,n4)

colSums(n) == colSums(eh)

# v = numero de animais observados pela ultima vez em i
v1 <- data.frame(table(l$l[l$plot==1]))
v2 <- data.frame(table(l$l[l$plot==2]))
v3 <- data.frame(table(l$l[l$plot==3]))
v4 <- data.frame(table(l$l[l$plot==4]))

v.all <- seq(1, n.occasions) # Preencher todos os anos ignorados
missing <- v.all[!(v.all %in% v1$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v1) <- names(df)
v1 <- rbind(v1, df)
v1$l<-as.numeric(levels(v1$l))[v1$l]
v1<-v1[order(v1$l), ]
v1<-v1$Freq
v1

missing <- v.all[!(v.all %in% v2$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v2) <- names(df)
v2 <- rbind(v2, df)
v2$l<-as.numeric(levels(v2$l))[v2$l]
v2<-v2[order(v2$l), ]
v2<-v2$Freq
v2

missing <- v.all[!(v.all %in% v3$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v3) <- names(df)
v3 <- rbind(v3, df)
v3$l<-as.numeric(levels(v3$l))[v3$l]
v3<-v3[order(v3$l), ]
v3<-v3$Freq
v3

missing <- v.all[!(v.all %in% v4$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v4) <- names(df)
v4 <- rbind(v4, df)
v4$l<-as.numeric(levels(v4$l))[v4$l]
v4<-v4[order(v4$l), ]
v4<-v4$Freq
v4

v <- rbind(v1,v2,v3,v4)
v

# d = numero de animais removidos da populacao no momento i
d <- rep(0,dim(eh)[2])
d <- rbind(d,d,d,d,d)
d

#Environmental variables
env.data <- readRDS("env_data_SGT.rds")
str(env.data)
unique(datA$TimeMonth)

env.data$TimeMonth[env.data$campanha==1] <- 1
env.data$TimeMonth[env.data$campanha==2] <- 5
env.data$TimeMonth[env.data$campanha==3] <- 13
env.data$TimeMonth[env.data$campanha==4] <- 16
env.data$TimeMonth

env.data <- env.data %>%   
  group_by(plot, TimeMonth) %>% 
  summarise(t.med = mean(t.med),
            t.max = mean(t.max),
            t.max.abs = mean(t.max.abs),
            rh.min.abs = mean(rh.min.abs),
            rh.max.abs = mean(rh.max.abs),
            VARI.all = mean(VARI.all),
            zentropy = mean(zentropy),
            tree.density = mean(tree.density),
            Ajalapensis_perf = mean(Ajalapensis_perf),
            Ajalapensis_ha90 = mean(Ajalapensis_ha90),
            TSLF = mean(TSLF))

TimeMonth.plot <- expand.grid(TimeMonth = 1:16,
                              plot = factor(c("C2", "C1", "QP", "QT"), 
                                            levels = c("C2", "C1", "QP", "QT")))

env.data <- full_join(env.data, TimeMonth.plot, by = c("plot", "TimeMonth"))

env.data <- env.data[order(env.data$TimeMonth),]
View(env.data)

amb <- array(c(rbind((env.data$t.med[env.data$plot=="C2"]-mean(env.data$t.med, na.rm = T))/sd(env.data$t.med, na.rm = T),
                     (env.data$t.med[env.data$plot=="C1"]-mean(env.data$t.med, na.rm = T))/sd(env.data$t.med, na.rm = T),
                     (env.data$t.med[env.data$plot=="QP"]-mean(env.data$t.med, na.rm = T))/sd(env.data$t.med, na.rm = T),
                     (env.data$t.med[env.data$plot=="QT"]-mean(env.data$t.med, na.rm = T))/sd(env.data$t.med, na.rm = T)),
               rbind((env.data$t.max[env.data$plot=="C2"]-mean(env.data$t.max, na.rm = T))/sd(env.data$t.max, na.rm = T),
                     (env.data$t.max[env.data$plot=="C1"]-mean(env.data$t.max, na.rm = T))/sd(env.data$t.max, na.rm = T),
                     (env.data$t.max[env.data$plot=="QP"]-mean(env.data$t.max, na.rm = T))/sd(env.data$t.max, na.rm = T),
                     (env.data$t.max[env.data$plot=="QT"]-mean(env.data$t.max, na.rm = T))/sd(env.data$t.max, na.rm = T)),
               rbind((env.data$t.max.abs[env.data$plot=="C2"]-mean(env.data$t.max.abs, na.rm = T))/sd(env.data$t.max.abs, na.rm = T),
                     (env.data$t.max.abs[env.data$plot=="C1"]-mean(env.data$t.max.abs, na.rm = T))/sd(env.data$t.max.abs, na.rm = T),
                     (env.data$t.max.abs[env.data$plot=="QP"]-mean(env.data$t.max.abs, na.rm = T))/sd(env.data$t.max.abs, na.rm = T),
                     (env.data$t.max.abs[env.data$plot=="QT"]-mean(env.data$t.max.abs, na.rm = T))/sd(env.data$t.max.abs, na.rm = T)),
               rbind((env.data$rh.min.abs[env.data$plot=="C2"]-mean(env.data$rh.min.abs, na.rm = T))/sd(env.data$rh.min.abs, na.rm = T),
                     (env.data$rh.min.abs[env.data$plot=="C1"]-mean(env.data$rh.min.abs, na.rm = T))/sd(env.data$rh.min.abs, na.rm = T),
                     (env.data$rh.min.abs[env.data$plot=="QP"]-mean(env.data$rh.min.abs, na.rm = T))/sd(env.data$rh.min.abs, na.rm = T),
                     (env.data$rh.min.abs[env.data$plot=="QT"]-mean(env.data$rh.min.abs, na.rm = T))/sd(env.data$rh.min.abs, na.rm = T)),
               rbind((env.data$rh.max.abs[env.data$plot=="C2"]-mean(env.data$rh.max.abs, na.rm = T))/sd(env.data$rh.max.abs, na.rm = T),
                     (env.data$rh.max.abs[env.data$plot=="C1"]-mean(env.data$rh.max.abs, na.rm = T))/sd(env.data$rh.max.abs, na.rm = T),
                     (env.data$rh.max.abs[env.data$plot=="QP"]-mean(env.data$rh.max.abs, na.rm = T))/sd(env.data$rh.max.abs, na.rm = T),
                     (env.data$rh.max.abs[env.data$plot=="QT"]-mean(env.data$rh.max.abs, na.rm = T))/sd(env.data$rh.max.abs, na.rm = T)),
               rbind((env.data$VARI.all[env.data$plot=="C2"]-mean(env.data$VARI.all, na.rm = T))/sd(env.data$VARI.all, na.rm = T),
                     (env.data$VARI.all[env.data$plot=="C1"]-mean(env.data$VARI.all, na.rm = T))/sd(env.data$VARI.all, na.rm = T),
                     (env.data$VARI.all[env.data$plot=="QP"]-mean(env.data$VARI.all, na.rm = T))/sd(env.data$VARI.all, na.rm = T),
                     (env.data$VARI.all[env.data$plot=="QT"]-mean(env.data$VARI.all, na.rm = T))/sd(env.data$VARI.all, na.rm = T)),
               rbind((env.data$zentropy[env.data$plot=="C2"]-mean(env.data$zentropy, na.rm = T))/sd(env.data$zentropy, na.rm = T),
                     (env.data$zentropy[env.data$plot=="C1"]-mean(env.data$zentropy, na.rm = T))/sd(env.data$zentropy, na.rm = T),
                     (env.data$zentropy[env.data$plot=="QP"]-mean(env.data$zentropy, na.rm = T))/sd(env.data$zentropy, na.rm = T),
                     (env.data$zentropy[env.data$plot=="QT"]-mean(env.data$zentropy, na.rm = T))/sd(env.data$zentropy, na.rm = T)),
               rbind((env.data$tree.density[env.data$plot=="C2"]-mean(env.data$tree.density, na.rm = T))/sd(env.data$tree.density, na.rm = T),
                     (env.data$tree.density[env.data$plot=="C1"]-mean(env.data$tree.density, na.rm = T))/sd(env.data$tree.density, na.rm = T),
                     (env.data$tree.density[env.data$plot=="QP"]-mean(env.data$tree.density, na.rm = T))/sd(env.data$tree.density, na.rm = T),
                     (env.data$tree.density[env.data$plot=="QT"]-mean(env.data$tree.density, na.rm = T))/sd(env.data$tree.density, na.rm = T)),
               rbind((env.data$Ajalapensis_perf[env.data$plot=="C2"]-mean(env.data$Ajalapensis_perf, na.rm = T))/sd(env.data$Ajalapensis_perf, na.rm = T),
                     (env.data$Ajalapensis_perf[env.data$plot=="C1"]-mean(env.data$Ajalapensis_perf, na.rm = T))/sd(env.data$Ajalapensis_perf, na.rm = T),
                     (env.data$Ajalapensis_perf[env.data$plot=="QP"]-mean(env.data$Ajalapensis_perf, na.rm = T))/sd(env.data$Ajalapensis_perf, na.rm = T),
                     (env.data$Ajalapensis_perf[env.data$plot=="QT"]-mean(env.data$Ajalapensis_perf, na.rm = T))/sd(env.data$Ajalapensis_perf, na.rm = T)),
               rbind((env.data$Ajalapensis_ha90[env.data$plot=="C2"]-mean(env.data$Ajalapensis_ha90, na.rm = T))/sd(env.data$Ajalapensis_ha90, na.rm = T),
                     (env.data$Ajalapensis_ha90[env.data$plot=="C1"]-mean(env.data$Ajalapensis_ha90, na.rm = T))/sd(env.data$Ajalapensis_ha90, na.rm = T),
                     (env.data$Ajalapensis_ha90[env.data$plot=="QP"]-mean(env.data$Ajalapensis_ha90, na.rm = T))/sd(env.data$Ajalapensis_ha90, na.rm = T),
                     (env.data$Ajalapensis_ha90[env.data$plot=="QT"]-mean(env.data$Ajalapensis_ha90, na.rm = T))/sd(env.data$Ajalapensis_ha90, na.rm = T)),
               rbind((env.data$TSLF[env.data$plot=="C2"] - mean(env.data$TSLF, na.rm = T))/sd(env.data$TSLF, na.rm = T),
                     (env.data$TSLF[env.data$plot=="C1"] - mean(env.data$TSLF, na.rm = T))/sd(env.data$TSLF, na.rm = T),
                     (env.data$TSLF[env.data$plot=="QP"] - mean(env.data$TSLF, na.rm = T))/sd(env.data$TSLF, na.rm = T),
                     (env.data$TSLF[env.data$plot=="QT"] - mean(env.data$TSLF, na.rm = T))/sd(env.data$TSLF, na.rm = T))),
             dim = c(4,16,11))
dim(amb)
str(amb)

# Valores iniciais
inits <- function (){ list(tauphib = 1, betaTphi = rep(0,11), varphi = rep(0,11),
                           sigma.phiJS = runif(1, 1, 2),
                           sigma.f = runif(1, 0.5, 1),
                           taufb = 1, betaTf = rep(0,11), varf = rep(0,11),
                           #alpha.pJS = runif(1, -0.5, 0.5),
                           taupb = 1,  betaTp = rep(0,11), varp = rep(0,11),
                           sigma.pJS = runif(1, 0.1, 0.5)
)}

# Defina os parametros a serem monitorados
parameters <- c("phiJS", "alpha.phiJS", "sigma.phiJS",
                "betaphiJS","varphi",
                "f", "alpha.f", "sigma.f",
                "betaf","varf",
                "pJS", "alpha.pJS", "sigma.pJS",
                "betap","varp",
                "rho")




sink("pradel-Ajalapensis-noenv.jags")
cat("

data{
for(j in 1:4){
C[j]<-10000
zeros[j]<-0
}}

model {

 
   #################
   #Pradel JS model#
   #################
   
###########PRIORS#######################
for(j in 1:4){
gamma[j, 1]<-0
phiJS[j, n.occasions]<-0
}

for(t in 1:n.occasions){
for(j in 1:4){
muJS[j,t]~dunif(0,1)
}}

#logit constraint for survival probability(phiJS)
for(j in 1:4){
alpha.phiJS[j] ~ dnorm(0.5,0.01)
mean.phiJS[j] <- 1/(1+exp(-alpha.phiJS[j]))#alpha.phiJS on prob scale
for(t in 1:(n.occasions-1)){
phiJS[j, t] <- 1/(1+exp(-logit.phiJS[j, t]))
logit.phiJS[j, t] <- alpha.phiJS[j] + eps.phiJS[j,t]
eps.phiJS[j,t] ~ dnorm(0,tau.phiJS)
}
}

#log constraint for recruitment rate(f)

for(j in 1:4){
alpha.f[j] ~ dnorm(-0.5,0.01)
mean.f[j] <- exp(alpha.f[j])#alpha.f on prob scale
for(t in 1:(n.occasions-1)){
f[j,t] <- exp(log.f[j,t])
log.f[j,t]<- alpha.f[j] + eps.f[j,t]
eps.f[j,t] ~ dnorm(0,tau.f)
}
}

for(j in 1:4){
alpha.pJS[j] ~ dnorm(-0.5,0.01)
mean.pJS[j] <- 1/(1+exp(-alpha.pJS[j])) #alpha.pJS on prob scale
for(t in 1:n.occasions){
#logit constraint for detectability (p)
pJS[j,t] <- (1/(1+exp(-logit.pJS[j,t]))) * missing[j,t]
logit.pJS[j,t] <- (alpha.pJS[j] + eps.pJS[j,t]) 
eps.pJS[j,t] ~ dnorm(0,tau.pJS)
}}


#temporal random variation
tau.f<-1/(sigma.f*sigma.f)
sigma.f~dunif(0,2)

tau.phiJS<-1/(sigma.phiJS*sigma.phiJS)
sigma.phiJS~dunif(0,2)

tau.pJS<-1/(sigma.pJS*sigma.pJS)
sigma.pJS~dunif(0,2)

###########LIKELIHOOD(ZERO-TRICK)######
for(j in 1:4){
zeros[j]~dpois(zero.mean[j])
zero.mean[j]<--LJS[j]+C[j]
LJS[j]<-sum(l.num[j, 1:n.occasions])-l.denom[j]

#####log-likelihood for the first occasion
l.num[j,1]<-(u[j,1]*log(xi[j,1]))+(n[j,1]*log(pJS[j,1]))+(secondexpo[j,1]*log(1-pJS[j,1]))+
(thirdexpo[j,1]*log(phiJS[j,1]))+(fourthexpo[j,1]*log(muJS[j,1]))+
(d[j,1]*log(1-muJS[j,1]))+(fifthexpo[j,1]*log(1-(pJS[j,1]*(1-muJS[j,1]))))+
(sixthexpo[j,1]*log(chi[j,1]))
xi[j,1]<-1
secondexpo_a[j,1]<-sum(u[j, 1:1])
secondexpo_b[j,1]<-0
secondexpo[j,1]<-secondexpo_a[j,1]-secondexpo_b[j,1]-n[j,1]
thirdexpo[j,1]<-sum(v[j,2:n.occasions])
fourthexpo[j,1]<-n[j,1]-d[j,1]
fifthexpo[j,1]<-sum(u[j,2:n.occasions])
sixthexpo[j,1]<-v[j,1]-d[j,1]

#####log-likelihood for the last occasion
l.num[j,n.occasions]<-(u[j,n.occasions]*log(xi[j,n.occasions]))+(firstexpo[j,n.occasions]*(log(phiJS[j,n.occasions-1])-log(phiJS[j,n.occasions-1]+f[j,n.occasions-1])))+
(n[j,n.occasions]*log(pJS[j,n.occasions]))+(secondexpo[j,n.occasions]*log(1-pJS[j,n.occasions]))+
(fourthexpo[j,n.occasions]*log(muJS[j,n.occasions]))+(d[j,n.occasions]*log(1-muJS[j,n.occasions]))+
(fifthexpo[j,n.occasions]*log(1-(pJS[j,n.occasions]*(1-muJS[j,n.occasions]))))+
(sixthexpo[j,n.occasions]*log(chi[j,n.occasions]))
chi[j,n.occasions]<-1

firstexpo[j,n.occasions]<-sum(u[j,1:(n.occasions-1)])
secondexpo_a[j,n.occasions]<-sum(u[j,1:n.occasions])
secondexpo_b[j,n.occasions]<-sum(v[j,1:(n.occasions-1)])
secondexpo[j,n.occasions]<-secondexpo_a[j,n.occasions]-secondexpo_b[j,n.occasions]-n[j,n.occasions]
fourthexpo[j,n.occasions]<-n[j,n.occasions]-d[j,n.occasions]
fifthexpo[j,n.occasions]<-0
sixthexpo[j,n.occasions]<-v[j,n.occasions]-d[j,n.occasions]
}

#####likelihood from occasion 2 to n.occasions-1
for(j in 1:4){
for(i in 2:(n.occasions-1)){
l.num[j,i]<-(u[j,i]*log(xi[j,i]))+(firstexpo[j,i]*(log(phiJS[j,i-1])-log(phiJS[j,i-1]+f[j,i-1])))+
(n[j,i]*log(pJS[j,i]))+(secondexpo[j,i]*log(1-pJS[j,i]))+
(thirdexpo[j,i]*log(phiJS[j,i]))+(fourthexpo[j,i]*log(muJS[j,i]))+
(d[j,i]*log(1-muJS[j,i]))+(fifthexpo[j,i]*log(1-(pJS[j,i]*(1-muJS[j,i]))))+
(sixthexpo[j,i]*log(chi[j,i]))

#first exponent
firstexpo[j,i]<-sum(u[j,1:(i-1)])

#second exponent
secondexpo_a[j,i]<-sum(u[j,1:i])
secondexpo_b[j,i]<-sum(v[j,1:(i-1)])
secondexpo[j,i]<-secondexpo_a[j,i]-secondexpo_b[j,i]-n[j,i]

#third exponent
thirdexpo[j,i]<-sum(v[j,(i+1):n.occasions])

#fourth exponent
fourthexpo[j,i]<-n[j,i]-d[j,i]

#fifth exponent
fifthexpo[j,i]<-sum(u[j,(i+1):n.occasions])

#sixth exponent
sixthexpo[j,i]<-v[j,i]-d[j,i]
}
}

#####likelihood denominator
#1st product
PROD1.1[1]<-1
PROD1.2[1]<-1
PROD1.3[1]<-1
PROD1.4[1]<-1

for(j in 1:(n.occasions-1)){
PROD1_tmp1[1,j]<-0
PROD1_tmp2[1,j]<-0
PROD1_tmp3[1,j]<-0
PROD1_tmp4[1,j]<-0
}

#fill part of PROD1_tmp
for(i in 2:(n.occasions-1)){
for(j in i:(n.occasions-1)){
PROD1_tmp1[i,j]<-0
PROD1_tmp2[i,j]<-0
PROD1_tmp3[i,j]<-0
PROD1_tmp4[i,j]<-0
}
}

for(i in 2:n.occasions){
for(j in 1:(i-1)){
PROD1_tmp1[i,j]<-phiJS[1,j]*(1-(pJS[1,j]*(1-muJS[1,j])))
PROD1_tmp2[i,j]<-phiJS[2,j]*(1-(pJS[2,j]*(1-muJS[2,j])))
PROD1_tmp3[i,j]<-phiJS[3,j]*(1-(pJS[3,j]*(1-muJS[3,j])))
PROD1_tmp4[i,j]<-phiJS[4,j]*(1-(pJS[4,j]*(1-muJS[4,j])))
}
}


PROD1.1[2]<-PROD1_tmp1[2,1]
PROD1.2[2]<-PROD1_tmp2[2,1]
PROD1.3[2]<-PROD1_tmp3[2,1]
PROD1.4[2]<-PROD1_tmp4[2,1]

for(i in 3:n.occasions){
PROD1.1[i]<-prod(PROD1_tmp1[i,1:(i-1)])
PROD1.2[i]<-prod(PROD1_tmp2[i,1:(i-1)])
PROD1.3[i]<-prod(PROD1_tmp3[i,1:(i-1)])
PROD1.4[i]<-prod(PROD1_tmp4[i,1:(i-1)])
}

#2nd product
PROD2.1[n.occasions]<-1
PROD2.2[n.occasions]<-1
PROD2.3[n.occasions]<-1
PROD2.4[n.occasions]<-1

for(i in 1:(n.occasions-1)){
for(j in (i+1):n.occasions){
PROD2_tmp1[i,j]<-gamma[1,j]
PROD2_tmp2[i,j]<-gamma[2,j]
PROD2_tmp3[i,j]<-gamma[3,j]
PROD2_tmp4[i,j]<-gamma[4,j]
}
}

#fill part of PROD2_tmp
for(i in 1:(n.occasions-1)){
for(j in 1:i){
PROD2_tmp1[i,j]<-0
PROD2_tmp2[i,j]<-0
PROD2_tmp3[i,j]<-0
PROD2_tmp4[i,j]<-0
}
}

PROD2.1[n.occasions-1]<-PROD2_tmp1[(n.occasions-1),n.occasions]
PROD2.2[n.occasions-1]<-PROD2_tmp2[(n.occasions-1),n.occasions]
PROD2.3[n.occasions-1]<-PROD2_tmp3[(n.occasions-1),n.occasions]
PROD2.4[n.occasions-1]<-PROD2_tmp4[(n.occasions-1),n.occasions]

for(i in 1:(n.occasions-2)){
PROD2.1[i]<-prod(PROD2_tmp1[i,(i+1):n.occasions])
PROD2.2[i]<-prod(PROD2_tmp2[i,(i+1):n.occasions])
PROD2.3[i]<-prod(PROD2_tmp3[i,(i+1):n.occasions])
PROD2.4[i]<-prod(PROD2_tmp4[i,(i+1):n.occasions])

}
for(i in 1:n.occasions){
denom_base_tmp1[i]<-xi[1,i]*PROD1.1[i]*PROD2.1[i]*pJS[1,i]
denom_base_tmp2[i]<-xi[2,i]*PROD1.2[i]*PROD2.2[i]*pJS[2,i]
denom_base_tmp3[i]<-xi[3,i]*PROD1.3[i]*PROD2.3[i]*pJS[3,i]
denom_base_tmp4[i]<-xi[4,i]*PROD1.4[i]*PROD2.4[i]*pJS[4,i]

}

denom_base1 <- sum(denom_base_tmp1[])
denom_base2 <- sum(denom_base_tmp2[])
denom_base3 <- sum(denom_base_tmp3[])
denom_base4 <- sum(denom_base_tmp4[])

denom_expo1 <- sum(u[1,1:n.occasions])
denom_expo2 <- sum(u[2,1:n.occasions])
denom_expo3 <- sum(u[3,1:n.occasions])
denom_expo4 <- sum(u[4,1:n.occasions])

l.denom[1] <- denom_expo1 * log(denom_base1)
l.denom[2] <- denom_expo2 * log(denom_base2)
l.denom[3] <- denom_expo3 * log(denom_base3)
l.denom[4] <- denom_expo4 * log(denom_base4)


#################Define xi and chi
for(i in 2:n.occasions){
for(j in 1:4){
xi.tmp[j,i]<-(1-gamma[j,i])+
(gamma[j,i]*((1-pJS[j,i-1])/(1-(pJS[j,i-1]*(1-muJS[j,i-1]))))*xi[j,i-1])
xi[j,i]<-max(xi.tmp[j,i],0.00001)
}
}

for(i in 1:(n.occasions-1)){
for(j in 1:4){
chi[j,i]<-(1-phiJS[j,i])+(phiJS[j,i]*(1-pJS[j,i+1])*chi[j,i+1])
}
}

#################Gamma and rho as derived parameter
for(i in 2:n.occasions){
for(j in 1:4){
rho[j,i]<-phiJS[j,i-1]+f[j,i-1]
}
}

for(i in 2:n.occasions){
for(j in 1:4){
gamma[j,i]<-phiJS[j,i-1]/(phiJS[j,i-1]+f[j,i-1])
}
}

}
",fill = TRUE)
sink()

# MCMC settings
ni <- 100000
nt <- 1
nb <- 200000
nc <- 4
na <- 50000

#Call JAGS from R (BRT 3 min)
bugs.data <- list(u = u, n = n, v = v, d = d, n.occasions = dim (eh)[2],
                  amb = amb, missing = missing.samp[1:4,])

bugs.data$missing <- as.matrix(bugs.data$missing)
runjags.options(jagspath = "/usr/local/bin/jags")

pradel.Ajalapensis <- run.jags(data=bugs.data, inits=inits, monitor=parameters, model="pradel-Ajalapensis-noenv.jags",
                         n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                         method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                         summarise = FALSE,
                         modules = c("glm"))

results.pradel.Ajalapensis <- results.jags(pradel.Ajalapensis)

results.pradel.Ajalapensis.df <- summary(results.pradel.Ajalapensis)
View(results.pradel.Ajalapensis.df)

write.csv(results.pradel.Ajalapensis.df, "results_pradel_Ajalapensis_df.csv")
saveRDS(results.pradel.Ajalapensis, "results_pradel_Ajalapensis.rds")

# Specify model in BUGS language
sink("pradel-Ajalapensis-env.jags")
cat("

data{
for(j in 1:4){
C[j]<-10000
zeros[j]<-0
}}

model {

 
   #################
   #Pradel JS model#
   #################
   
###########PRIORS#######################
for(j in 1:4){
gamma[j, 1]<-0
phiJS[j, n.occasions]<-0
}

for(t in 1:n.occasions){
for(j in 1:4){
muJS[j,t]~dunif(0,1)
}}

#logit constraint for survival probability(phiJS)
for(j in 1:4){
alpha.phiJS[j] ~ dnorm(0.5,0.01)
mean.phiJS[j] <- 1/(1+exp(-alpha.phiJS[j]))#alpha.phiJS on prob scale
for(t in 1:(n.occasions-1)){
phiJS[j, t] <- 1/(1+exp(-logit.phiJS[j, t]))
logit.phiJS[j, t] <- alpha.phiJS[j] + eps.phiJS[j,t] + inprod(amb[j,t,],betaphiJS)
eps.phiJS[j,t] ~ dnorm(0,tau.phiJS)
}
}

for(j in 1:11){
varphi[j]~dbern(0.5)
betaTphi[j]~dnorm(0,tauphib)
betaphiJS[j]<-varphi[j]*betaTphi[j]
}

#environmental parameters
tauphib~dgamma(1,0.001)

#log constraint for recruitment rate(f)

for(j in 1:4){
alpha.f[j] ~ dnorm(-0.5,0.01)
mean.f[j] <- exp(alpha.f[j])#alpha.f on prob scale
for(t in 1:(n.occasions-1)){
f[j,t] <- exp(log.f[j,t])
log.f[j,t]<- alpha.f[j] + eps.f[j,t]+ inprod(amb[j,t,],betaf)
eps.f[j,t] ~ dnorm(0,tau.f)
}
}

for(j in 1:11){
varf[j]~dbern(0.5)
betaTf[j]~dnorm(0,taufb)
betaf[j]<-varf[j]*betaTf[j]
}

#environmental parameters
taufb~dgamma(1,0.001)

for(j in 1:4){
alpha.pJS[j] ~ dnorm(-0.5,0.01)
mean.pJS[j] <- 1/(1+exp(-alpha.pJS[j])) #alpha.pJS on prob scale
for(t in 1:n.occasions){
#logit constraint for detectability (p)
pJS[j,t] <- (1/(1+exp(-logit.pJS[j,t])))
logit.pJS[j,t] <- alpha.pJS[j] + eps.pJS[j,t] + inprod(amb[j,t,],betap) 
eps.pJS[j,t] ~ dnorm(0,tau.pJS)
}}

for(j in 1:11){
varp[j]~dbern(0.5)
betaTp[j]~dnorm(0,taupb)
betap[j]<-varp[j]*betaTp[j]
}

#environmental parameters
taupb~dgamma(1,0.001)

#temporal random variation
tau.f<-1/(sigma.f*sigma.f)
sigma.f~dunif(0,2)

tau.phiJS<-1/(sigma.phiJS*sigma.phiJS)
sigma.phiJS~dunif(0,2)

tau.pJS<-1/(sigma.pJS*sigma.pJS)
sigma.pJS~dunif(0,2)

###########LIKELIHOOD(ZERO-TRICK)######
for(j in 1:4){
zeros[j]~dpois(zero.mean[j])
zero.mean[j]<--LJS[j]+C[j]
LJS[j]<-sum(l.num[j, 1:n.occasions])-l.denom[j]

#####log-likelihood for the first occasion
l.num[j,1]<-(u[j,1]*log(xi[j,1]))+(n[j,1]*log(pJS[j,1]))+(secondexpo[j,1]*log(1-pJS[j,1]))+
(thirdexpo[j,1]*log(phiJS[j,1]))+(fourthexpo[j,1]*log(muJS[j,1]))+
(d[j,1]*log(1-muJS[j,1]))+(fifthexpo[j,1]*log(1-(pJS[j,1]*(1-muJS[j,1]))))+
(sixthexpo[j,1]*log(chi[j,1]))
xi[j,1]<-1
secondexpo_a[j,1]<-sum(u[j, 1:1])
secondexpo_b[j,1]<-0
secondexpo[j,1]<-secondexpo_a[j,1]-secondexpo_b[j,1]-n[j,1]
thirdexpo[j,1]<-sum(v[j,2:n.occasions])
fourthexpo[j,1]<-n[j,1]-d[j,1]
fifthexpo[j,1]<-sum(u[j,2:n.occasions])
sixthexpo[j,1]<-v[j,1]-d[j,1]

#####log-likelihood for the last occasion
l.num[j,n.occasions]<-(u[j,n.occasions]*log(xi[j,n.occasions]))+(firstexpo[j,n.occasions]*(log(phiJS[j,n.occasions-1])-log(phiJS[j,n.occasions-1]+f[j,n.occasions-1])))+
(n[j,n.occasions]*log(pJS[j,n.occasions]))+(secondexpo[j,n.occasions]*log(1-pJS[j,n.occasions]))+
(fourthexpo[j,n.occasions]*log(muJS[j,n.occasions]))+(d[j,n.occasions]*log(1-muJS[j,n.occasions]))+
(fifthexpo[j,n.occasions]*log(1-(pJS[j,n.occasions]*(1-muJS[j,n.occasions]))))+
(sixthexpo[j,n.occasions]*log(chi[j,n.occasions]))
chi[j,n.occasions]<-1

firstexpo[j,n.occasions]<-sum(u[j,1:(n.occasions-1)])
secondexpo_a[j,n.occasions]<-sum(u[j,1:n.occasions])
secondexpo_b[j,n.occasions]<-sum(v[j,1:(n.occasions-1)])
secondexpo[j,n.occasions]<-secondexpo_a[j,n.occasions]-secondexpo_b[j,n.occasions]-n[j,n.occasions]
fourthexpo[j,n.occasions]<-n[j,n.occasions]-d[j,n.occasions]
fifthexpo[j,n.occasions]<-0
sixthexpo[j,n.occasions]<-v[j,n.occasions]-d[j,n.occasions]
}

#####likelihood from occasion 2 to n.occasions-1
for(j in 1:4){
for(i in 2:(n.occasions-1)){
l.num[j,i]<-(u[j,i]*log(xi[j,i]))+(firstexpo[j,i]*(log(phiJS[j,i-1])-log(phiJS[j,i-1]+f[j,i-1])))+
(n[j,i]*log(pJS[j,i]))+(secondexpo[j,i]*log(1-pJS[j,i]))+
(thirdexpo[j,i]*log(phiJS[j,i]))+(fourthexpo[j,i]*log(muJS[j,i]))+
(d[j,i]*log(1-muJS[j,i]))+(fifthexpo[j,i]*log(1-(pJS[j,i]*(1-muJS[j,i]))))+
(sixthexpo[j,i]*log(chi[j,i]))

#first exponent
firstexpo[j,i]<-sum(u[j,1:(i-1)])

#second exponent
secondexpo_a[j,i]<-sum(u[j,1:i])
secondexpo_b[j,i]<-sum(v[j,1:(i-1)])
secondexpo[j,i]<-secondexpo_a[j,i]-secondexpo_b[j,i]-n[j,i]

#third exponent
thirdexpo[j,i]<-sum(v[j,(i+1):n.occasions])

#fourth exponent
fourthexpo[j,i]<-n[j,i]-d[j,i]

#fifth exponent
fifthexpo[j,i]<-sum(u[j,(i+1):n.occasions])

#sixth exponent
sixthexpo[j,i]<-v[j,i]-d[j,i]
}
}

#####likelihood denominator
#1st product
PROD1.1[1]<-1
PROD1.2[1]<-1
PROD1.3[1]<-1
PROD1.4[1]<-1

for(j in 1:(n.occasions-1)){
PROD1_tmp1[1,j]<-0
PROD1_tmp2[1,j]<-0
PROD1_tmp3[1,j]<-0
PROD1_tmp4[1,j]<-0
}

#fill part of PROD1_tmp
for(i in 2:(n.occasions-1)){
for(j in i:(n.occasions-1)){
PROD1_tmp1[i,j]<-0
PROD1_tmp2[i,j]<-0
PROD1_tmp3[i,j]<-0
PROD1_tmp4[i,j]<-0
}
}

for(i in 2:n.occasions){
for(j in 1:(i-1)){
PROD1_tmp1[i,j]<-phiJS[1,j]*(1-(pJS[1,j]*(1-muJS[1,j])))
PROD1_tmp2[i,j]<-phiJS[2,j]*(1-(pJS[2,j]*(1-muJS[2,j])))
PROD1_tmp3[i,j]<-phiJS[3,j]*(1-(pJS[3,j]*(1-muJS[3,j])))
PROD1_tmp4[i,j]<-phiJS[4,j]*(1-(pJS[4,j]*(1-muJS[4,j])))
}
}


PROD1.1[2]<-PROD1_tmp1[2,1]
PROD1.2[2]<-PROD1_tmp2[2,1]
PROD1.3[2]<-PROD1_tmp3[2,1]
PROD1.4[2]<-PROD1_tmp4[2,1]

for(i in 3:n.occasions){
PROD1.1[i]<-prod(PROD1_tmp1[i,1:(i-1)])
PROD1.2[i]<-prod(PROD1_tmp2[i,1:(i-1)])
PROD1.3[i]<-prod(PROD1_tmp3[i,1:(i-1)])
PROD1.4[i]<-prod(PROD1_tmp4[i,1:(i-1)])
}

#2nd product
PROD2.1[n.occasions]<-1
PROD2.2[n.occasions]<-1
PROD2.3[n.occasions]<-1
PROD2.4[n.occasions]<-1

for(i in 1:(n.occasions-1)){
for(j in (i+1):n.occasions){
PROD2_tmp1[i,j]<-gamma[1,j]
PROD2_tmp2[i,j]<-gamma[2,j]
PROD2_tmp3[i,j]<-gamma[3,j]
PROD2_tmp4[i,j]<-gamma[4,j]
}
}

#fill part of PROD2_tmp
for(i in 1:(n.occasions-1)){
for(j in 1:i){
PROD2_tmp1[i,j]<-0
PROD2_tmp2[i,j]<-0
PROD2_tmp3[i,j]<-0
PROD2_tmp4[i,j]<-0
}
}

PROD2.1[n.occasions-1]<-PROD2_tmp1[(n.occasions-1),n.occasions]
PROD2.2[n.occasions-1]<-PROD2_tmp2[(n.occasions-1),n.occasions]
PROD2.3[n.occasions-1]<-PROD2_tmp3[(n.occasions-1),n.occasions]
PROD2.4[n.occasions-1]<-PROD2_tmp4[(n.occasions-1),n.occasions]

for(i in 1:(n.occasions-2)){
PROD2.1[i]<-prod(PROD2_tmp1[i,(i+1):n.occasions])
PROD2.2[i]<-prod(PROD2_tmp2[i,(i+1):n.occasions])
PROD2.3[i]<-prod(PROD2_tmp3[i,(i+1):n.occasions])
PROD2.4[i]<-prod(PROD2_tmp4[i,(i+1):n.occasions])

}
for(i in 1:n.occasions){
denom_base_tmp1[i]<-xi[1,i]*PROD1.1[i]*PROD2.1[i]*pJS[1,i]
denom_base_tmp2[i]<-xi[2,i]*PROD1.2[i]*PROD2.2[i]*pJS[2,i]
denom_base_tmp3[i]<-xi[3,i]*PROD1.3[i]*PROD2.3[i]*pJS[3,i]
denom_base_tmp4[i]<-xi[4,i]*PROD1.4[i]*PROD2.4[i]*pJS[4,i]

}

denom_base1 <- sum(denom_base_tmp1[])
denom_base2 <- sum(denom_base_tmp2[])
denom_base3 <- sum(denom_base_tmp3[])
denom_base4 <- sum(denom_base_tmp4[])

denom_expo1 <- sum(u[1,1:n.occasions])
denom_expo2 <- sum(u[2,1:n.occasions])
denom_expo3 <- sum(u[3,1:n.occasions])
denom_expo4 <- sum(u[4,1:n.occasions])

l.denom[1] <- denom_expo1 * log(denom_base1)
l.denom[2] <- denom_expo2 * log(denom_base2)
l.denom[3] <- denom_expo3 * log(denom_base3)
l.denom[4] <- denom_expo4 * log(denom_base4)


#################Define xi and chi
for(i in 2:n.occasions){
for(j in 1:4){
xi.tmp[j,i]<-(1-gamma[j,i])+
(gamma[j,i]*((1-pJS[j,i-1])/(1-(pJS[j,i-1]*(1-muJS[j,i-1]))))*xi[j,i-1])
xi[j,i]<-max(xi.tmp[j,i],0.00001)
}
}

for(i in 1:(n.occasions-1)){
for(j in 1:4){
chi[j,i]<-(1-phiJS[j,i])+(phiJS[j,i]*(1-pJS[j,i+1])*chi[j,i+1])
}
}

#################Gamma and rho as derived parameter
for(i in 2:n.occasions){
for(j in 1:4){
rho[j,i]<-phiJS[j,i-1]+f[j,i-1]
}
}

for(i in 2:n.occasions){
for(j in 1:4){
gamma[j,i]<-phiJS[j,i-1]/(phiJS[j,i-1]+f[j,i-1])
}
}

}
",fill = TRUE)
sink()

bugs.data <- list(u = u[,c(1,5,13,16)], n = n[,c(1,5,13,16)], v = v[,c(1,5,13,16)], d = d[,c(1,5,13,16)], n.occasions = dim (eh[,c(1,5,13,16)])[2],
                  amb = amb[,c(1,5,13,16),])

pradel.Ajalapensis <- run.jags(data=bugs.data, inits=inits, monitor=parameters, model="pradel-Ajalapensis-env.jags",
                               n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                               method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                               summarise = FALSE,
                               modules = c("glm"))

results.pradel.Ajalapensis <- results.jags(pradel.Ajalapensis)

results.pradel.Ajalapensis.df <- summary(results.pradel.Ajalapensis)
View(results.pradel.Ajalapensis.df)

write.csv(results.pradel.Ajalapensis.df, "results_pradel_env_Ajalapensis_df.csv")
saveRDS(results.pradel.Ajalapensis, "results_pradel_env_Ajalapensis.rds")

#Plots
pradel.Ajalapensis.df <- read.csv("results_pradel_Ajalapensis_df.csv")

f.pradel <- pradel.Ajalapensis.df[grep(pattern = "f", x = pradel.Ajalapensis.df$X)[1:60],]

phi.pradel <- pradel.Ajalapensis.df[grep(pattern = "phi", x = pradel.Ajalapensis.df$X)[1:60],]

p.pradel <- pradel.Ajalapensis.df[grep(pattern = "pJS", x = pradel.Ajalapensis.df$X)[1:64],]


f.pradel$plot <- rep(1:4,15)
f.pradel$time <- rep(1:15, each = 4)

phi.pradel$plot <- rep(1:4,15)
phi.pradel$time <- rep(1:15, each = 4)

p.pradel$plot <- rep(1:4,16)
p.pradel$time <- rep(1:16, each = 4)


phi.pradel$plot <- as.factor(phi.pradel$plot)
f.pradel$plot <- as.factor(f.pradel$plot)
p.pradel$plot <- as.factor(p.pradel$plot)


ggplot(phi.pradel, aes(x = time, y = Mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = Mean,colour=plot), alpha=0.5, linewidth = 2) +
  geom_ribbon(aes(x = time, ymin = Lower95, ymax = Upper95, fill=plot), colour = NA, alpha = 0.2)+
  #ylim(c(0.75,1))+
  scale_color_manual(values=turbo(4))+
    scale_fill_manual(values=turbo(4))+
  labs(x = "Months", y = "Survival")

ggplot(p.pradel, aes(x = time, y = Mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = Mean,colour=plot), alpha=0.5, linewidth = 2) +
  geom_ribbon(aes(x = time, ymin = Lower95, ymax = Upper95, fill=plot), colour = NA, alpha = 0.2)+
  #ylim(c(0.75,1))+
  scale_color_manual(values=turbo(4))+
  scale_fill_manual(values=turbo(4))+
  labs(x = "Months", y = "Capture probability")

ggplot(f.pradel, aes(x = time, y = Mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = Mean,colour=plot), alpha=0.5, linewidth = 2) +
  geom_ribbon(aes(x = time, ymin = Lower95, ymax = Upper95, fill=plot), colour = NA, alpha = 0.2)+
  #ylim(c(0.75,1))+
  scale_color_manual(values=turbo(4))+
  scale_fill_manual(values=turbo(4))+
  labs(x = "Months", y = "Recruitment")

f.pradel.mean <- pradel.Ajalapensis.df[grep(pattern = "f", x = pradel.Ajalapensis.df$X)[61:64],]

phi.pradel.mean <- pradel.Ajalapensis.df[grep(pattern = "phi", x = pradel.Ajalapensis.df$X)[65:68],]

p.pradel.mean <- pradel.Ajalapensis.df[grep(pattern = "pJS", x = pradel.Ajalapensis.df$X)[65:68],]

f.pradel.mean$plot <- as.factor(1:4)
phi.pradel.mean$plot <- as.factor(1:4)
p.pradel.mean$plot <- as.factor(1:4)

#Rescale
f.pradel.mean$Mean <- exp(f.pradel.mean$Mean)
f.pradel.mean$Lower95 <- exp(f.pradel.mean$Lower95)
f.pradel.mean$Upper95 <- exp(f.pradel.mean$Upper95)

phi.pradel.mean$Mean <- plogis(phi.pradel.mean$Mean)
phi.pradel.mean$Lower95 <- plogis(phi.pradel.mean$Lower95)
phi.pradel.mean$Upper95 <- plogis(phi.pradel.mean$Upper95)

p.pradel.mean$Mean <- plogis(p.pradel.mean$Mean)
p.pradel.mean$Lower95 <- plogis(p.pradel.mean$Lower95)
p.pradel.mean$Upper95 <- plogis(p.pradel.mean$Upper95)

quartz(height = 8, width = 8)
ggplot(phi.pradel.mean, aes(x = plot, colour = plot))+
  geom_pointrange(aes(x = plot, y = Mean, ymin = Lower95, ymax = Upper95, colour=plot), size = 1.5, linewidth = 1.2)+
  scale_color_manual(values=turbo(4), name = "Fire severity")+
  labs(x = "", y = "Survival")

quartz(height = 8, width = 8)
ggplot(f.pradel.mean, aes(x = plot, colour = plot))+
  geom_pointrange(aes(x = plot, y = Mean, ymin = Lower95, ymax = Upper95, colour=plot), size = 1.5, linewidth = 1.2)+
scale_color_manual(values=turbo(4), name = "Fire severity")+
  labs(x = "", y = "Recruitment")

quartz(height = 8, width = 8)
ggplot(p.pradel.mean, aes(x = plot, colour = plot))+
  geom_pointrange(aes(x = plot, y = Mean, ymin = Lower95, ymax = Upper95, colour=plot), size = 1.5, linewidth = 1.2)+
  scale_color_manual(values=turbo(4), name = "Fire severity")+
  labs(x = "", y = "Capture probability")


# SVL with environment ----------------------------------------------------

env.data <- readRDS("env_data_SGT.rds")

svl.data <- dados[,c("Campanha", "Data", "Mes", "Ano", "Parcela", "Gride", "Numerodetombo", "Recaptura", 
                     "Massa", "CRC", "CC", "BC", "Quebrada", "Sexo", "Ovos", "Morto")]

names(svl.data) <- c("campanha", "date", "month", "year", "plot", "trap", "ID", "recapture", 
                     "mass", "svl", "tl", "tb", "broken", "sex", "eggs", "dead")

svl.data$trap <- as.factor(paste(svl.data$plot, svl.data$trap, sep = "_"))
svl.data$campanha <- as.factor(svl.data$campanha)

svl.env <- left_join(svl.data, env.data, by = c("campanha", "plot", "trap"))
svl.env[,c(18:28)] <- scale(svl.env[,c(18:28)])

head(svl.env)

hist((svl.env$svl))

mfull.svl <- brm(svl ~ plot + TSLF + t.med + t.max + t.max.abs + rh.min.abs + rh.max.abs + 
                VARI.all + zentropy + tree.density + Ajalapensis_perf + Ajalapensis_ha90,
                data = svl.env, cores = 4)

summary(mfull.svl)
plot(mfull.svl)
bayes_R2(mfull.svl)

# perform variable selection without cross-validation
vs <- cv_varsel(mfull.svl, validate_search = F)
summary(vs)
summary(vs, deltas = T)

suggest_size(vs)
ranking(vs)

quartz(height = 8, width = 8)
plot(vs, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.svl)

# perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits <- run_cvfun(
  refm_obj,
  ### Only for the sake of speed (not recommended in general):
  K = 10
  ###
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

quartz(height = 8, width = 8)
plot(cv_proportions(ranking(cvvs), cumulate = T))

quartz(height = 8, width = 8)
plot(cvvs, deltas = T, stats = "mlpd") + theme_minimal()


msel.svl.re <- brm(svl ~ Ajalapensis_perf + (1|campanha/plot/trap),
                   data = svl.env, cores = 4, iter = 4000,
                   control = list(adapt_delta = 0.99))

summary(msel.svl.re)
bayes_R2(msel.svl.re)

#Diagnostic plots
plot(msel.svl.re)
pp_check(msel.svl.re)

#Conditional effects plot
p1 <- conditional_effects(msel.svl.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$Ajalapensis_perf, aes(x = Ajalapensis_perf, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = svl.env, aes(x = Ajalapensis_perf, y = svl, color = plot), alpha = 0.5, size = 2)+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Locomotor performance", y="Snout-vent length (mm)") +
  theme_minimal()
  
#Testing differences among plots
msel.svl.plot.re <- brm(svl ~ plot + (1|campanha/trap),
                        data = svl.env, cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.99))

summary(msel.svl.plot.re)
bayes_R2(msel.svl.plot.re)

#Diagnostic plots
plot(msel.svl.plot.re)
pp_check(msel.svl.plot.re)

#Conditional effects plot
p1 <- conditional_effects(msel.svl.plot.re)

quartz(height = 8, width = 8)
ggplot(svl.env, 
       aes(x = plot, y = svl, colour = plot)) + 
  ggdist::stat_halfeye(aes(fill = plot),
    adjust = 5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA,
    alpha = 0.5
  ) +
  gghalves::geom_half_point(aes(colour = plot),
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  geom_pointrange(data = p1$plot, aes(x = plot, y = estimate__, 
                                      ymin = lower__, ymax = upper__, 
                                      colour=plot), 
                  size = 1.2, linewidth = 1.2) +

  geom_hline(yintercept = 40,linetype="dashed")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Snout-vent length (mm)")+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  scale_fill_manual(values = turbo(4), name = "Fire severity") + 
  theme_minimal()


emmeans::emmeans(msel.svl.plot.re, pairwise ~ plot)


# Abundance with environment -----------------------------------------------

str(svl.env)


svl.env <- svl.env %>%
  group_by(campanha) %>%
  mutate(sampling_day = dense_rank(date)) %>%  # Rank the dates within each month
  ungroup()

summary(svl.env$sampling_day)

trap_locations <- fire.regimes.arms %>%
  distinct(plot, trap, treatment)

trap_locations$campaign <- substr(trap_locations$trap,1,2)
trap_locations$campaign[trap_locations$plot=="C1"] <- NA
trap_locations$campaign[trap_locations$plot=="C2"] <- NA
trap_locations$campaign[trap_locations$plot=="QP"] <- NA
trap_locations$campaign[trap_locations$plot=="QT"] <- NA

trap_locations_heitor <- trap_locations[is.na(trap_locations$campaign),]
trap_locations_bruna <- trap_locations[!is.na(trap_locations$campaign),]

trap_locations_heitor <- expand.grid.df(trap_locations_heitor[,-4], 
                                        data.frame(campaign = as.character(c(1:4))))

trap_locations <- rbind(trap_locations_bruna, trap_locations_heitor)

# Now, create the full dataset of all possible sampling events
# `crossing` will combine every row from `trap_locations` with every sampling day.
all_sampling_events <- crossing(trap_locations, sampling_day = 1:15)

# Count captures only for existing combinations in your data
capture_counts <- svl.env %>%
  count(campanha, plot, trap, sampling_day, name = "freq")

names(capture_counts) <- c("campaign", "plot", "trap", "sampling_day", "freq")
str(capture_counts)
str(all_sampling_events)

capture_counts$campaign <- as.character(capture_counts$campaign)
capture_counts$plot <- as.character(capture_counts$plot)
capture_counts$trap <- as.character(capture_counts$trap)

# Join the actual captures to the master list of all sampling events
Ajalapensis.captures.day <- all_sampling_events %>%
  left_join(capture_counts, by = c("campaign", "plot", "trap", "sampling_day")) %>%
  # If a trap on a given day had no captures, its freq is NA. Change that to 0.
  mutate(freq = replace_na(freq, 0))

# View the result
print(Ajalapensis.captures.day)

summary(Ajalapensis.captures.day)
sum(Ajalapensis.captures.day$freq)


Ajalapensis.captures.day <- Ajalapensis.captures.day %>%  
  pivot_wider(names_from = sampling_day, 
              values_from = freq, 
              names_prefix = "day_")

summary(Ajalapensis.captures.day)
str(Ajalapensis.captures.day)

env.data$plot <- as.character(env.data$plot)
env.data$trap <- as.character(env.data$trap)
env.data$campaign <- as.character(env.data$campanha)

Ajalapensis.captures.day.env <- left_join(env.data,
                                          Ajalapensis.captures.day, 
                                          by = c("campaign", "plot", "trap"))
str(env.data)
str(Ajalapensis.captures.day.env)
summary(Ajalapensis.captures.day.env)

## brms------------------------------------------------------------------------

capts.env <- data.frame(Ajalapensis.captures.day.env[,-c(20:34)],
                   capts = rowSums(Ajalapensis.captures.day.env[,c(20:34)]))

capts.env[,c(4:15)] <- scale(capts.env[,c(4:15)])


mfull.capts <- brm(capts ~ plot + TSLF + t.med + t.max + t.max.abs + rh.min.abs + rh.max.abs + 
                   VARI.all + zentropy + tree.density + Ajalapensis_perf + Ajalapensis_ha90,
                   data = capts.env, cores = 4, family = "poisson")

summary(mfull.capts)
plot(mfull.capts)
bayes_R2(mfull.capts)

# perform variable selection without cross-validation
vs_capts <- cv_varsel(mfull.capts, validate_search = F)
summary(vs_capts)
summary(vs_capts, deltas = T)

suggest_size(vs_capts)

quartz(height = 8, width = 8)
plot(vs_capts, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.capts)

# perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_capts <- run_cvfun(
  refm_obj,
  K = 10
  ###
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

quartz(height = 8, width = 8)
plot(cv_proportions(ranking(cvvs_capts), cumulate = T))

quartz(height = 8, width = 8)
plot(cvvs_capts, deltas = T, stats = "mlpd") + theme_minimal()


msel.capts.re <- brm(capts ~ Ajalapensis_ha90 + (1|campaign/plot/trap),
                     data = capts.env, cores = 4, family = "poisson",
                     iter = 5000, control = list(adapt_delta = 0.999))

summary(msel.capts.re)
bayes_R2(msel.capts.re)

plot(msel.capts.re)

pp_check(msel.capts.re)

#Conditional effects plot
p1 <- conditional_effects(msel.capts.re)

capts.env$plot <- factor(capts.env$plot,
                         levels = c("C2", "C1", "QP", "QT"))

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$Ajalapensis_ha90, aes(x = Ajalapensis_ha90, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = capts.env, aes(x = Ajalapensis_ha90, y = capts, color = plot), alpha = 0.5, size = 2)+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Hours of activity", y="Captures") + 
  theme_minimal()



msel.capts.plot.re <- brm(capts ~ plot + (1|campaign/trap),
                        data = capts.env, cores = 4, family = "poisson",
                        iter = 4000, control = list(adapt_delta = 0.99))

summary(msel.capts.plot.re)
plot(msel.capts.plot.re)

pp_check(msel.capts.plot.re)
bayes_R2(msel.capts.plot.re)


p1 <- conditional_effects(msel.capts.plot.re)

quartz(height = 8, width = 8)
ggplot(capts.env, 
       aes(x = plot, y = capts, colour = plot)) + 
  ggdist::stat_halfeye(aes(fill = plot),
                       adjust = 5,
                       width = .5,
                       .width = 0,
                       justification = -.2,
                       point_colour = NA,
                       alpha = 0.5
  ) +
  gghalves::geom_half_point(aes(colour = plot),
                            ## draw jitter on the left
                            side = "l", 
                            ## control range of jitter
                            range_scale = .5, 
                            ## add some transparency
                            alpha = .3
  ) +
  geom_pointrange(data = p1$plot, aes(x = plot, y = estimate__, 
                                      ymin = lower__, ymax = upper__, 
                                      colour=plot), 
                  size = 1.2, linewidth = 1.2) +
  coord_cartesian( clip = "off")+
  labs(x="", y="Captures")+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  theme_minimal()

emmeans::emmeans(msel.capts.plot.re, pairwise ~ plot)




## INLA------------------------------------------------------------------------

(y.mat <- as.matrix(Ajalapensis.captures.day.env[,c(20:34)]))
counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    Ajalapensis.captures.day.env$campaign,
                                    paste(Ajalapensis.captures.day.env$campaign, 
                                          Ajalapensis.captures.day.env$plot, sep = "_"),
                                    paste(Ajalapensis.captures.day.env$campaign, 
                                          Ajalapensis.captures.day.env$plot, 
                                          Ajalapensis.captures.day.env$trap, sep = "_"),
                                    scale(Ajalapensis.captures.day.env$Ajalapensis_ha90))

print(counts.and.count.covs)

# Only intercept
out.inla.env0 <- inla(counts.and.count.covs[,c(1:16)] ~ 1,
                   data = list(counts.and.count.covs = counts.and.count.covs[,c(1:16)]),
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                                      theta2 = list(prior = "flat",
                                                                    param = numeric()))),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.env0, digits = 3)


plot(out.inla.env0,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env0,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env0,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = T)


# Only random effects
out.inla.env1 <- inla(counts.and.count.covs[,c(1:16)] ~ 1 + 
                        f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs[,c(1:16)],
                               campaign = counts.and.count.covs$X2,
                               plot_id = counts.and.count.covs$X3),
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.1,
                                        prec.intercept = 0.1),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, .1)),
                                                      theta2 = list(param = c(0, .1)), 
                                                      theta3 = list(prior = "gamma", param = c(2, 0.5)),
                                                      theta4 = list(prior = "pc.prec", param = c(1, 0.1)),
                                                      theta5 = list(param = c(0, .1))
                                                      
                   )
                   ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.env1, digits = 3)


plot(out.inla.env1,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

# Hours of activity effect

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    scale(Ajalapensis.captures.day.env$Ajalapensis_ha90),
                                    paste(Ajalapensis.captures.day.env$campaign, 
                                          Ajalapensis.captures.day.env$plot, sep = "_"))

print(counts.and.count.covs)

out.inla.env2 <- inla(counts.and.count.covs ~ 1 + ha_90 + f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs,
                               plot_id = counts.and.count.covs$X3,
                               ha_90 = counts.and.count.covs$X2),
                   
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                      theta2 = list(param = c(0, 0.1)), 
                                                      theta3 = list(param = c(0, 0.1)),
                                                      theta4 = list(param = c(0, 0.1)),
                                                      theta5 = list(param = c(0, 0.1)),
                                                      theta6 = list(param = c(0, 0.1))
                   )
                   ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.env2, digits = 3)


plot(out.inla.env2,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.env2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE) 

# Plots effects
dummy_vars <- model.matrix(~ factor(Ajalapensis.captures.day.env$plot, 
                                    levels = c("C2", "C1","QP","QT")))

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1, 
                                    dummy_vars[, 2], 
                                    dummy_vars[, 3], 
                                    dummy_vars[, 4]
)

print(counts.and.count.covs)




out.inla.env3 <- inla(counts.and.count.covs ~ 1 + plot + f(campaign, model = "iid"),
                      data = list(
                        counts.and.count.covs = counts.and.count.covs,
                        plot = factor(Ajalapensis.captures.day.env$plot, 
                                      levels = c("C2", "C1","QP","QT")),          # 'plot' factor for the detection model
                        campaign = Ajalapensis.captures.day.env$campaign   # 'campaign' for the detection model
                      ),
                      
                      family = "nmixnb",
                      #Priors for detection parameters
                      control.fixed = list(mean = 0, 
                                           mean.intercept = 0, 
                                           prec = 0.01,
                                           prec.intercept = 0.01),
                      #Priors for abundance and overdispersion parameters
                      control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                         theta2 = list(param = c(0, 0.1)), 
                                                         theta3 = list(param = c(0, 0.1)),
                                                         theta4 = list(param = c(0, 0.1)),
                                                         theta5 = list(param = c(0, 0.1)),
                                                         theta6 = list(param = c(0, 0.1))
                      )
                      ),
                      verbose = TRUE,
                      control.inla = list(int.strategy = "eb"),
                      control.compute=list(config = TRUE, 
                                           waic = TRUE,
                                           residuals = TRUE,
                                           cpo = TRUE))
summary(out.inla.env3, digits = 3)


plot(out.inla.env3,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.env3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE) 

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1)

print(counts.and.count.covs)

out.inla.env4 <- inla(counts.and.count.covs ~ 1 + f(campaign, model = "iid"),
                      data = list(
                        counts.and.count.covs = counts.and.count.covs,
                        campaign = Ajalapensis.captures.day.env$campaign   # 'campaign' for the detection model
                      ),
                      
                      family = "nmixnb",
                      #Priors for detection parameters
                      control.fixed = list(mean = 0, 
                                           mean.intercept = 0, 
                                           prec = 0.01,
                                           prec.intercept = 0.01),
                      #Priors for abundance and overdispersion parameters
                      control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                         theta2 = list(param = c(0, 0.1)), 
                                                         theta3 = list(param = c(0, 0.1)),
                                                         theta4 = list(param = c(0, 0.1)),
                                                         theta5 = list(param = c(0, 0.1)),
                                                         theta6 = list(param = c(0, 0.1))
                      )
                      ),
                      verbose = TRUE,
                      control.inla = list(int.strategy = "eb"),
                      control.compute=list(config = TRUE, 
                                           waic = TRUE,
                                           residuals = TRUE,
                                           cpo = TRUE))
summary(out.inla.env4, digits = 3)


plot(out.inla.env4,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.env4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.env4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE) 


out.inla.env0$waic$waic
out.inla.env1$waic$waic
out.inla.env2$waic$waic
out.inla.env3$waic$waic
out.inla.env4$waic$waic

inla_env_models <- list(out.inla.env0, 
                        out.inla.env1, 
                        out.inla.env2, 
                        out.inla.env3,
                        out.inla.env4)

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
  Weight = round(w_akaike,2)
)

waic_table$formula[1] <- "(lambda ~ 1), (p ~ 1)"
waic_table$formula[2] <- "(lambda ~ 1), (p ~ (1|campaign:plot))"
waic_table$formula[3] <- "(lambda ~ ha90), (p ~ ha90 + (1|campaign:plot))"
waic_table$formula[4] <- "(lambda ~ plot), (p ~ plot + (1|campaign))"
waic_table$formula[5] <- "(lambda ~ 1), (p ~ 1 + (1|campaign))"

# Sort by WAIC (best model first)
waic_table <- waic_table[order(waic_table$WAIC), ]

# Print table
print(waic_table)

summary(out.inla.env2)
summary(out.inla.env3)

summary(plogis(as.matrix(out.inla.env2$summary.linear.predictor[,c("mean", 
                                                                "sd", 
                                                                "0.025quant", 
                                                                "0.5quant",
                                                                "0.975quant",
                                                                "mode")])))
summary(out.inla.env2$summary.fitted.values)


# Hours of activity effects

ha_xx <- seq(min(scale(Ajalapensis.captures.day.env$Ajalapensis_ha90)), 
               max(scale(Ajalapensis.captures.day.env$Ajalapensis_ha90)), 
               by = 0.1)

# Draw 1000 samples from the model's posterior
posterior_samples <- inla.posterior.sample(10000, out.inla.env2)

# The 'contents' table describes the structure of the model's output
contents <- out.inla.env2$misc$configs$contents

# Find the numeric index for the intercept
intercept_idx <- contents$start[contents$tag == "(Intercept)"]

# Find the numeric index for the slope
slope_idx <- contents$start[contents$tag == "ha_90"]

# Step 3: Run the corrected sapply loop
predicted_lines <- sapply(posterior_samples, function(s) {
  # Access the column by its position [ ,1] instead of by name
  intercept <- s$latent[intercept_idx, 1]
  slope <- s$latent[slope_idx, 1]
  
  # The rest of the code is the same
  linear_predictor <- intercept + slope * freq_xx
  return(plogis(linear_predictor))
})

rel_p_ha_90 <- data.frame(ha_90 = ha_xx,
                         median = apply(predicted_lines, 1, quantile, 0.5),
                         lower = apply(predicted_lines, 1, quantile, 0.025),
                         upper = apply(predicted_lines, 1, quantile, 0.975))

pred_df <- data.frame(ha_90 = scale(Ajalapensis.captures.day.env$Ajalapensis_ha90),
                      trap = Ajalapensis.captures.day.env$trap,
                      plot = factor(Ajalapensis.captures.day.env$plot,
                                    levels = c("C2", "C1", "QP", "QT")),
                      p = plogis(out.inla.env2$summary.linear.predictor$`0.5quant`),
                      ncaps = rowSums(y.mat))

quartz(h = 8, w = 8)
ggplot(pred_df, aes(x = ha_90, y = p)) +
  geom_point(data = pred_df, aes(x = ha_90, y = p, colour = plot), size = 3, alpha = 0.3) +
  scale_color_manual(values = turbo(4)) +
  geom_line(data = rel_p_ha_90, aes(x = ha_90, y = median), color = "blue", linewidth = 2) +
  geom_ribbon(data = rel_p_ha_90, aes(ymin = lower, ymax = upper), alpha = 0.3) +
  labs(x = "Hours of activity", y = "Capture") +
  theme_minimal()

ggplot(pred_df, aes(x = ncaps, y = p)) +
  geom_point(size = 3, alpha = 0.3)

ggplot(pred_df, aes(x = ha_90, y = ncaps)) +
  geom_point(size = 3, alpha = 0.3)

# Plots effects

# Capture

#--- Step 1: Get posterior samples for the fixed effects ---
# We'll use the robust index method to avoid errors.
post_samples_fixed <- inla.posterior.sample(10000, out.inla.env3)
contents <- out.inla.env3$misc$configs$contents

#--- Step 2: Get the numeric index for each fixed effect ---
idx_intercept <- contents$start[contents$tag == "(Intercept)"]
idx_C1 <- contents$start[contents$tag == "plotC1"]
idx_QP <- contents$start[contents$tag == "plotQP"]
idx_QT <- contents$start[contents$tag == "plotQT"]

#--- Step 3: Extract the posterior samples for each coefficient ---
# This loop gets the values for each sample and stores them in vectors
b_intercept <- sapply(post_samples_fixed, function(s) s$latent[idx_intercept, 1])
b_C1 <- sapply(post_samples_fixed, function(s) s$latent[idx_C1, 1])
b_QP <- sapply(post_samples_fixed, function(s) s$latent[idx_QP, 1])
b_QT <- sapply(post_samples_fixed, function(s) s$latent[idx_QT, 1])

#--- Step 4: Calculate all 6 pairwise contrasts on the logit scale ---
# The coefficients themselves are already differences relative to the intercept.
# To compare other levels, we subtract their respective coefficients.
contrasts_p <- list(
  `C1 - C2` = b_C1,
  `QP - C2` = b_QP,
  `QT - C2` = b_QT,
  `C1 - QP` = b_C1 - b_QP,
  `C1 - QT` = b_C1 - b_QT,
  `QP - QT` = b_QP - b_QT
)

#--- Step 5: Summarize the results into a data frame ---
contrast_summary_p <- t(sapply(contrasts_p, quantile, probs = c(0.025, 0.5, 0.975)))

contrast_summary_p_df <- as.data.frame(contrast_summary_p)
colnames(contrast_summary_p_df) <- c("Lower.HPD", "Median", "Upper.HPD")

print(contrast_summary_p_df, digits = 3)

#--- Step 2: Calculate detection probability for each plot type ---
# The model is: logit(p) = Intercept + effect_C1 + effect_QP + effect_QT
# We use plogis() to transform back to the probability scale


# Loop through samples to get posteriors for p
p_C2_post <- sapply(post_samples_fixed, function(s) plogis(s$latent[idx_intercept, 1]))
p_C1_post <- sapply(post_samples_fixed, function(s) plogis(s$latent[idx_intercept, 1] + s$latent[idx_C1, 1]))
p_QP_post <- sapply(post_samples_fixed, function(s) plogis(s$latent[idx_intercept, 1] + s$latent[idx_QP, 1]))
p_QT_post <- sapply(post_samples_fixed, function(s) plogis(s$latent[idx_intercept, 1] + s$latent[idx_QT, 1]))

#--- Step 3: Summarize the posteriors for plotting ---

p_post <- data.frame(
  C2 = p_C2_post,
  C1 = p_C1_post,
  QP = p_QP_post,
  QT = p_QT_post
)


p_post_long <- p_post %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Plot", 
    values_to = "p"
  ) %>% 
  mutate(Plot = factor(Plot, levels = c("C2", "C1", "QP", "QT")))

quartz(h = 8, w = 8)
ggplot(p_post_long, aes(x = Plot, y = p, fill = Plot, color = Plot)) +
  
  # FIX 1: Use the correct function name 'stat_halfeye'
  ggdist::stat_halfeye(
    adjust = 1,
    width = .8,
    justification = -.1,
    alpha = 0.5,
    .width = 0.95,  # Sets the credible interval to 95%
    point_interval = "median_hdi" # Display the median and Highest Density Interval
  ) +
  
  # FIX 2: Wrap the labels argument in a function
  scale_y_log10() +
  
  # Add other labels and scales
  labs(
    y = "Capture probability",
    x = "Fire regime"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  guides(fill = "none") # Hide the redundant legend

# Abundance

#--- Step 1: Get samples from the posterior hyperparameters ---
# Increase sample size for more stable estimates
post_samples <- inla.hyperpar.sample(10000, out.inla.env3)

# Extract the posterior samples for each beta coefficient for abundance.
# Use the exact names from your model's summary output.
beta1 <- post_samples[, "beta[1] for NMixNB observations"] # Intercept (C2)
beta2 <- post_samples[, "beta[2] for NMixNB observations"] # C1 vs C2
beta3 <- post_samples[, "beta[3] for NMixNB observations"] # QP vs C2
beta4 <- post_samples[, "beta[4] for NMixNB observations"] # QT vs C2

#--- Step 2: Calculate all 6 pairwise contrasts ---
contrasts <- list(
  `C1 - C2` = beta2,
  `QP - C2` = beta3,
  `QT - C2` = beta4,
  `C1 - QP` = beta2 - beta3,
  `C1 - QT` = beta2 - beta4,
  `QP - QT` = beta3 - beta4
)

#--- Step 3: Summarize the results into a data frame ---
contrast_summary <- t(sapply(contrasts, quantile, probs = c(0.025, 0.5, 0.975)))

# Make it a clean data frame and print it
contrast_summary_df <- as.data.frame(contrast_summary)
colnames(contrast_summary_df) <- c("Lower.HPD", "Median", "Upper.HPD")

print(contrast_summary_df, digits = 3)

#--- Step 1: Calculate lambda for each plot type using the posterior samples ---
# These are on the log scale, so we use exp() to transform them back
lambda_post <- data.frame(
  C2 = exp(beta1),
  C1 = exp(beta1 + beta2),
  QP = exp(beta1 + beta3),
  QT = exp(beta1 + beta4)
)


lambda_post_long <- lambda_post %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Plot", 
    values_to = "lambda"
  ) %>% 
  mutate(Plot = factor(Plot, levels = c("C2", "C1", "QP", "QT")))

y_limits <- quantile(lambda_post_long$lambda, c(0.025, 0.975))
y_min <- y_limits[1]
y_max <- y_limits[2]

quartz(h = 8, w = 8)
ggplot(lambda_post_long, aes(x = Plot, y = lambda, fill = Plot, color = Plot)) +
  
  # FIX 1: Use the correct function name 'stat_halfeye'
  ggdist::stat_halfeye(
    adjust = 1,
    width = .8,
    justification = -.1,
    alpha = 0.5,
    .width = 0.95,  # Sets the credible interval to 95%
    point_interval = "median_hdi" # Display the median and Highest Density Interval
  ) +
  
  # FIX 2: Wrap the labels argument in a function
  scale_y_log10(
    labels = function(x) number(x, accuracy = 1)
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
  guides(fill = "none") # Hide the redundant legend


# Sex with environment -----------------------------------------------
head(svl.env)

svl.env$sex[svl.env$sex=="I"] <- NA

table(svl.env$sex)

mfull.sex <- brm(sex ~ plot + TSLF + t.med + t.max + t.max.abs + rh.min.abs + rh.max.abs + 
                     VARI.all + zentropy + tree.density + Ajalapensis_perf + Ajalapensis_ha90,
                   data = svl.env, cores = 4, family = "bernoulli")

summary(mfull.sex)
plot(mfull.sex)
# plot(conditional_effects(mfull.sex))

# perform variable selection without cross-validation
vs_sex <- cv_varsel(mfull.sex, validate_search = F)
summary(vs_sex)
summary(vs_sex, deltas = T)

suggest_size(vs_sex)

quartz(height = 8, width = 8)
plot(vs_sex, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.sex)

# perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_sex <- run_cvfun(
  refm_obj,
  ### Only for the sake of speed (not recommended in general):
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

quartz(height = 8, width = 8)
plot(cv_proportions(ranking(cvvs_sex), cumulate = T))

quartz(height = 8, width = 8)
plot(cvvs_sex, deltas = T, stats = "mlpd") + theme_minimal()


msel.sex.re <- brm(sex ~ Ajalapensis_perf + (1|campanha/plot/trap),
                   data = svl.env, cores = 4, family = "bernoulli",
                   iter = 4000, control = list(adapt_delta = 0.99))

summary(msel.sex.re)
bayes_R2(msel.sex.re)

msel.null.re <- brm(sex ~ 1 + (1|campanha/plot/trap),
                   data = svl.env, cores = 4, family = "bernoulli",
                   iter = 4000, control = list(adapt_delta = 0.99))

summary(msel.null.re)
bayes_R2(msel.null.re)

msel.sex.plot.re <- brm(sex ~ plot + (1|campanha/trap),
                        data = svl.env, cores = 4, family = "bernoulli",
                        iter = 4000, control = list(adapt_delta = 0.99))

summary(msel.sex.plot.re)
plot(msel.sex.plot.re)

pp_check(msel.sex.plot.re)
bayes_R2(msel.sex.plot.re)

loo(msel.sex.re, msel.null.re, msel.sex.plot.re)

conditional_effects(msel.sex.re)
par(mfrow=c(1,1))
cdplot(as.factor(sex) ~ Ajalapensis_perf, data = svl.env)

svl.env$sex <- as.integer(as.factor(svl.env$sex))-1


# p1 <- conditional_effects(msel.sex.re)
# 
# 
# 
# #Salvar!
# windows(height = 8, width = 8)
# quartz(height = 8, width = 8)
# ggplot(p1$Ajalapensis_perf, aes(x = Ajalapensis_perf, y  = estimate__)) +
#   geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
#   geom_line(color = "blue", linewidth = 1) +
#   geom_point(data = svl.env, aes(x = Ajalapensis_perf, y = sex), alpha = 0.5, size = 2)+
#   # scale_colour_manual(values = turbo(4), name = "Fire severity")+
#   labs(x="Locomotor performance", y="Proportion of males")

p1 <- conditional_effects(msel.sex.plot.re)


quartz(height = 8, width = 8)
ggplot(svl.env, 
       aes(x = plot, y = sex, colour = plot)) + 
  # gghalves::geom_half_point(aes(y = jitter(sex, 0.2),colour = plot),
  #                           ## draw jitter on the left
  #                           side = "l", 
  #                           ## control range of jitter
  #                           range_scale = .5, 
  #                           ## add some transparency
  #                           alpha = .3
  # ) +
  ggdist::geom_dotsinterval(aes(fill = plot))+
# ggdist::stat_histinterval(point_color = NA, interval_color = NA, 
#                           aes(fill = plot), 
#                          justification = -0.05)+
  geom_pointrange(data = p1$plot, aes(x = plot, y = estimate__, 
                                      ymin = lower__, ymax = upper__, 
                                      colour=plot), 
                  size = 1.2, linewidth = 1.2) +
    coord_cartesian(clip = "off")+
  labs(x="", y="Proportion of males")+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  scale_fill_manual(values = turbo(4), name = "Fire severity")

emmeans::emmeans(msel.sex.plot.re, pairwise ~ plot)

table(svl.env$sex, svl.env$plot)


# Body condition with environment -----------------------------------------------

# Remove casos onde massa = NA
completos <- complete.cases(svl.env[,c("mass")])
condicao.dados <- droplevels(svl.env[completos,])

# Visualiza dados (crc, massa)
plot(condicao.dados$svl, condicao.dados$mass, las=1, bty="n")
plot(log(condicao.dados$svl), log(condicao.dados$mass), las=1, bty="n")


# Calcula a condicao corporal como "scaled mass index" (Peig & Green, 2009)
# install.packages("lmodel2", dependencies=T)
library(lmodel2)

sma.jalapensis <- lmodel2(log(mass+1) ~ log(svl), data=condicao.dados, "relative", "relative", 99)
sma.jalapensis
par(mfrow = c(1,2))

plot(sma.jalapensis, "OLS")
plot(sma.jalapensis, "SMA")
plot(sma.jalapensis, "MA")
plot(sma.jalapensis, "RMA")

par(mfrow = c(1,1))

attach(condicao.dados)
smi <- mass*((mean(svl)/svl)^2.027180)  #2.027180 = slope of SMA regression
boxplot(smi)
plot(svl, smi, las=1, bty="n")
boxplot(smi~sex)
detach(condicao.dados)
condicao.dados$smi <- smi
rm(smi)
detach("package:lmodel2", unload = TRUE)

# Testa efeitos das variaveis ambientais sobbre a condicao corporal

mfull.smi <- brm(smi ~ plot + TSLF + t.med + t.max + t.max.abs + rh.min.abs + rh.max.abs + 
                   VARI.all + zentropy + tree.density + Ajalapensis_perf + Ajalapensis_ha90,
                 data = condicao.dados, cores = 4)

summary(mfull.smi)
plot(mfull.smi)
bayes_R2(mfull.smi)

# perform variable selection without cross-validation
vs_smi <- cv_varsel(mfull.smi, validate_search = F)
summary(vs_smi)
summary(vs_smi, deltas = T)

suggest_size(vs_smi)

quartz(height = 8, width = 8)
plot(vs_smi, deltas = T, stats = "mlpd") + theme_minimal()

refm_obj <- get_refmodel(mfull.smi)

# perform variable selection with cross-validation
# Refit the reference model K times:
cv_fits_smi <- run_cvfun(
  refm_obj,
  ### Only for the sake of speed (not recommended in general):
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

quartz(height = 8, width = 8)
plot(cv_proportions(ranking(cvvs_smi), cumulate = T))

quartz(height = 8, width = 8)
plot(cvvs_smi, deltas = T, stats = "mlpd") + theme_minimal()


msel.smi.re <- brm(smi ~ Ajalapensis_perf + (1|campanha/plot/trap),
                   data = condicao.dados, cores = 4, iter = 5000,
                   control = list(adapt_delta = 0.999),
                   silent = 0)

summary(msel.smi.re)
plot(msel.smi.re)

pp_check(msel.smi.re)
bayes_R2(msel.smi.re)

#Conditional effects plot
p1 <- conditional_effects(msel.smi.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$Ajalapensis_perf, aes(x = Ajalapensis_perf, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = condicao.dados, aes(x = Ajalapensis_perf, y = smi, color = plot), alpha = 0.5, size = 2)+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Locomotor performance", y="Scaled mass index") +
  theme_minimal()



msel.smi.plot.re <- brm(smi ~ plot + (1|campanha/trap),
                        data = condicao.dados, cores = 4, iter = 5000,
                        control = list(adapt_delta = 0.999))

summary(msel.smi.plot.re)
plot(msel.smi.plot.re)

pp_check(msel.smi.plot.re)
bayes_R2(msel.smi.plot.re)

p1 <- conditional_effects(msel.smi.plot.re)

quartz(height = 8, width = 8)
ggplot(condicao.dados, 
       aes(x = plot, y = smi, colour = plot)) + 
  ggdist::stat_halfeye(aes(fill = plot),
                       adjust = 5,
                       width = .5,
                       .width = 0,
                       justification = -.2,
                       point_colour = NA,
                       alpha = 0.5
  ) +
  gghalves::geom_half_point(aes(colour = plot),
                            ## draw jitter on the left
                            side = "l", 
                            ## control range of jitter
                            range_scale = .5, 
                            ## add some transparency
                            alpha = .3
  ) +
  geom_pointrange(data = p1$plot, aes(x = plot, y = estimate__, 
                                      ymin = lower__, ymax = upper__, 
                                      colour=plot), 
                  size = 1.2, linewidth = 1.2) +
  
  coord_cartesian( clip = "off")+
  labs(x="", y="Scaled mass index")+
  scale_colour_manual(values = turbo(4), name = "Fire severity")+
  scale_fill_manual(values = turbo(4), name = "Fire severity") +
  theme_minimal()

emmeans::emmeans(msel.smi.plot.re, pairwise ~ plot)


# Fire regime components --------------------------------------------------

brunadata <- readxl::read_excel("Ameivula_jalapensis_EESGT_BrunaGomes.xlsx", na = "NA")
glimpse(brunadata)

brunadata <- dplyr::rename(brunadata, "plot" = "sampling_point")
brunadata$trap <- paste0(brunadata$campaign, brunadata$Y)
glimpse(brunadata)
summary(brunadata)
brunadata$plot <- as.character(brunadata$plot)

fire.regimes.arms <- read.csv("fire_regimes_arms_df.csv")

brunadata <- left_join(brunadata, fire.regimes.arms[,-1], by = c("plot", "trap", "treatment"))
summary(brunadata)

## SVL -------------------------------------------------------------------


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

svlbruna <- brunadata[,c("campaign", "date", "plot", "aiq_code", "TUF", "treatment", "treat_code",
                         "recapture", "svl_mm","weight_animal", "cicatriz_umbilical", "sex","ovada",
                         "trap", "freq", "MeanTSLF", "MFRI","severity", "TSLF")]

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

names(svlbruna) <- c("campaign", "date", "plot", "aiq_code", "TUF", "treatment", "treat_code",
                     "recapture", "svl","mass", "umb_scar", "sex","eggs",
                     "trap", "freq", "MeanTSLF", "MFRI", "severity", "TSLF")

# Number of recaptures
table(svlbruna$recapture)

# Number of captures (individuals)
nrow(svlbruna) - table(svlbruna$recapture)

# Mean recapture rate
table(svlbruna$recapture)/nrow(svlbruna)

# Sex
table(svlbruna$sex)

svl.data <- dados[,c("Campanha", "Data", "Mes", "Ano", "Parcela", "Gride", "Numerodetombo", "Recaptura", 
                     "CRC", "Massa",  "Sexo", "Ovos")]

names(svl.data) <- c("campaign", "date", "month", "year", "plot", "trap",  
                     "ID", "recapture", "svl", "mass", "sex", "eggs")

svl.data$trap <- as.factor(paste(svl.data$plot, svl.data$trap, sep = "_"))
svl.data$campaign <- as.factor(svl.data$campaign)

svl.fire <- left_join(svl.data, fire.regimes.arms[,-c(1:2, 10)], by = c("plot", "trap"))
env.data$campaign <- env.data$campanha
svl.fire <- left_join(svl.fire, env.data[,c(1:3,5,18)], by = c("plot", "trap", "campaign"))

svl.fire$date <- as.POSIXct(svl.fire$date, format = "%d/%m/%Y", tz = "UTC")
summary(svl.fire)
unique(svl.fire$severity)
plot(svl ~ severity, data = svl.fire)

svl.merged <- dplyr::full_join(svl.fire, svlbruna)
svl.merged$month <- month(svl.merged$date)
svl.merged$year <- year(svl.merged$date)

summary(svl.merged)
table(svl.merged$TUF == svl.merged$TSLF)

svl.merged$TSLF[!is.na(svl.merged$TUF)] <- svl.merged$TSLF[!is.na(svl.merged$TUF)]*12
svl.merged$TSLF[svl.merged$treatment=="Q"] <- 0

ggplot(svl.merged, aes(y = svl, x = jitter(severity))) +
  geom_point(alpha = 0.5) +
  geom_smooth()+ 
  geom_quantile(color = "red", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9))

ggplot(svl.merged, aes(y = svl, x = jitter(freq))) +
  geom_point(alpha = 0.5) +
  geom_smooth()+ 
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

quartz(h = 8, w = 12)
ggplot(svl.merged, aes(y = svl, x = as.Date(date))) +
  labs(x = "Date", y = "Snout-vent length (mm)") +
  scale_x_date(date_minor_breaks = "1 month") +
  theme_minimal() +
  geom_hline(aes(yintercept = 40), linetype = 2, alpha = 0.5)+
  geom_point(alpha = 0.5) 

quartz(h = 8, w = 8)
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


quartz(h = 8, w = 8)
ggpairs(svl.merged[,c(13:16,19)])
usdm::vifstep(svl.merged[,c(13:16,19)], keep = "severity", th = 2)
usdm::vif(svl.merged[,c(13:16,19)])

svl.merged[,c(13:16,19)] <- scale(svl.merged[,c(13:16,19)])

# Testing the effects of fire components
msel.svl.MeanTSLF.re <- brm(svl ~ MeanTSLF + (1|campaign/plot/trap),
                            data = svl.merged, 
                            cores = 4)

summary(msel.svl.MeanTSLF.re)
bayes_R2(msel.svl.MeanTSLF.re)

#Diagnostic plots
plot(msel.svl.MeanTSLF.re)
pp_check(msel.svl.MeanTSLF.re)

msel.svl.TSLF.re <- brm(svl ~ TSLF + (1|campaign/plot/trap),
                            data = svl.merged, 
                            cores = 4)

summary(msel.svl.TSLF.re)
bayes_R2(msel.svl.TSLF.re)

#Diagnostic plots
plot(msel.svl.TSLF.re)
pp_check(msel.svl.TSLF.re)

msel.svl.severity.re <- brm(svl ~ severity + (1|campaign/plot/trap),
                        data = svl.merged, 
                        cores = 4)

summary(msel.svl.severity.re)
bayes_R2(msel.svl.severity.re)

#Diagnostic plots
plot(msel.svl.severity.re)
pp_check(msel.svl.severity.re)

msel.svl.firenull.re <- brm(svl ~ 1 + (1|campaign/plot/trap),
                        data = svl.merged, 
                        cores = 4)

summary(msel.svl.firenull.re)
bayes_R2(msel.svl.firenull.re)

#Diagnostic plots
plot(msel.svl.firenull.re)
pp_check(msel.svl.firenull.re)

loo(msel.svl.MeanTSLF.re, msel.svl.TSLF.re, 
    msel.svl.severity.re, msel.svl.firenull.re)


#Conditional effects plot
p1 <- conditional_effects(msel.svl.MeanTSLF.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$MeanTSLF, aes(x = MeanTSLF, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = 40), linetype = 2) +
  geom_point(data = svl.merged, aes(x = MeanTSLF, y = svl), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Mean fire interval", y="Snout-vent length (mm)") +
  theme_minimal()

#Conditional effects plot
p1 <- conditional_effects(msel.svl.severity.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$severity, aes(x = severity, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = 40), linetype = 2) +
  geom_point(data = svl.merged, aes(x = severity, y = svl), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Mean fire severity", y="Snout-vent length (mm)") +
  theme_minimal()


## Body condition -----------------------------------------------

# Remove casos onde massa = NA
completos <- complete.cases(svl.merged[,c("mass")])
condicao.dados <- droplevels(svl.merged[completos,])

# Visualiza dados (crc, massa)
plot(condicao.dados$svl, condicao.dados$mass, las=1, bty="n")
plot(log(condicao.dados$svl), log(condicao.dados$mass), las=1, bty="n")


mod1 <- lm(condicao.dados$mass ~ condicao.dados$svl)
summary(mod1)

mod2 <- lm(sqrt(condicao.dados$mass) ~ condicao.dados$svl)
summary(mod2)

mod3 <- lm((condicao.dados$mass ^ (1/3)) ~ condicao.dados$svl)
summary(mod3)

mod4 <- lm(log(condicao.dados$mass ^ (1/3)) ~ log(condicao.dados$svl))
summary(mod4)

AIC(mod1, mod2, mod3, mod4)

plot(condicao.dados$svl,
     log(condicao.dados$mass^(1/3)),
     pch = 21,
     cex = 1,
     col = rgb(217, 72, 1, 100, maxColorValue = 255),
     bg = rgb(253, 174, 107, 100, maxColorValue = 255),
     las = 1,
     bty = "n",
     ylab = "log(Massa ^ 1/3 (g))",
     xlab = "Comprimento Rostro-Cloacal (mm)") +
abline(a = lm(log(condicao.dados$mass^(1/3)) ~ condicao.dados$svl), col = "darkred")

zcrit <- qnorm(0.9999)
summary(rstandard(mod4))
st_residuals <- condicao.dados %>% 
  filter(abs(rstandard(mod4)) > zcrit)
st_residuals

# Remove outliers
condicao.dados <- condicao.dados %>% 
  filter(abs(rstandard(mod4)) < zcrit)

mod4 <- lm(log(condicao.dados$mass ^ (1/3)) ~ log(condicao.dados$svl))
summary(mod4)

plot(condicao.dados$svl,
     log(condicao.dados$mass^(1/3)),
     pch = 21,
     cex = 1,
     col = rgb(217, 72, 1, 100, maxColorValue = 255),
     bg = rgb(253, 174, 107, 100, maxColorValue = 255),
     las = 1,
     bty = "n",
     ylab = "log(Massa ^ 1/3 (g))",
     xlab = "Comprimento Rostro-Cloacal (mm)") +
  abline(a = lm(log(condicao.dados$mass^(1/3)) ~ condicao.dados$svl), col = "darkred")

st_residuals <- condicao.dados %>% 
  filter(abs(rstandard(mod4)) > zcrit)
st_residuals

p <- length(coefficients(mod4))
n <- length(fitted(mod4))
ratio <- p / n
leverage <- condicao.dados %>% 
  filter(hatvalues(mod4) > (3 * ratio))
leverage

cutoff <- 4 / (dim(condicao.dados)[1])
cooks_dis <- condicao.dados %>% 
  filter(cooks.distance(mod4) > cutoff) 
cooks_dis

bonf <- car::outlierTest(mod4,
                         n.max = 9999)

condicao.dados[names(bonf$rstudent),]

# Calcula a condicao corporal como "scaled mass index" (Peig & Green, 2009)
# install.packages("lmodel2", dependencies=T)
library(lmodel2)

sma.jalapensis <- lmodel2(log(mass+1) ~ log(svl), data=condicao.dados, "relative", "relative", 99)
sma.jalapensis
par(mfrow = c(1,2))

plot(sma.jalapensis, "OLS")
plot(sma.jalapensis, "SMA")
plot(sma.jalapensis, "MA")
plot(sma.jalapensis, "RMA")

par(mfrow = c(1,1))

attach(condicao.dados)
smi <- mass*((mean(svl)/svl)^1.929048)  #1.929048= slope of SMA regression
boxplot(smi)
plot(svl, smi, las=1, bty="n")
boxplot(smi~sex)
detach(condicao.dados)
condicao.dados$smi <- smi
rm(smi)
detach("package:lmodel2", unload = TRUE)

# Testing the effects of fire components
msel.smi.MeanTSLF.re <- brm(smi ~ MeanTSLF + (1|campaign/plot/trap),
                            data = condicao.dados, 
                            cores = 4)

summary(msel.smi.MeanTSLF.re)
bayes_R2(msel.smi.MeanTSLF.re)

#Diagnostic plots
plot(msel.smi.MeanTSLF.re)
pp_check(msel.smi.MeanTSLF.re)

msel.smi.TSLF.re <- brm(smi ~ TSLF + (1|campaign/plot/trap),
                        data = condicao.dados, 
                        cores = 4)

summary(msel.smi.TSLF.re)
bayes_R2(msel.smi.TSLF.re)

#Diagnostic plots
plot(msel.smi.TSLF.re)
pp_check(msel.smi.TSLF.re)

msel.smi.severity.re <- brm(smi ~ severity + (1|campaign/plot/trap),
                            data = condicao.dados, 
                            cores = 4)

summary(msel.smi.severity.re)
bayes_R2(msel.smi.severity.re)

#Diagnostic plots
plot(msel.smi.severity.re)
pp_check(msel.smi.severity.re)

msel.smi.firenull.re <- brm(smi ~ 1 + (1|campaign/plot/trap),
                            data = condicao.dados, 
                            cores = 4)

summary(msel.smi.firenull.re)
bayes_R2(msel.smi.firenull.re)

#Diagnostic plots
plot(msel.smi.firenull.re)
pp_check(msel.smi.firenull.re)

loo(msel.smi.MeanTSLF.re, msel.smi.TSLF.re,
    msel.smi.severity.re, msel.smi.firenull.re)


#Conditional effects plot
p1 <- conditional_effects(msel.smi.TSLF.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$TSLF, aes(x = TSLF, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = condicao.dados, aes(x = TSLF, y = smi), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Time since last fire", y="Scaled mass index") +
  theme_minimal()

#Conditional effects plot
p1 <- conditional_effects(msel.smi.severity.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$severity, aes(x = severity, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = condicao.dados, aes(x = severity, y = smi), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Mean fire severity", y="Scaled mass index") +
  theme_minimal()

## Sex ratio -----------------------------------------------
head(svl.merged)
table(svl.merged$sex)

svl.merged$sex[svl.merged$sex=="I"] <- NA
svl.merged$sex[svl.merged$sex=="j"] <- NA
svl.merged$sex[svl.merged$sex=="J"] <- NA
svl.merged$sex[svl.merged$sex=="f"] <- "F"
svl.merged$sex[svl.merged$sex=="m"] <- "M"

table(svl.merged$sex)

cdplot(as.factor(sex) ~ TSLF, data = svl.merged)
cdplot(as.factor(sex) ~ MeanTSLF, data = svl.merged)
cdplot(as.factor(sex) ~ severity, data = svl.merged)
cdplot(as.factor(sex) ~ freq, data = svl.merged)

# Testing the effects of fire components
msel.sex.MeanTSLF.re <- brm(sex ~ MeanTSLF + (1|campaign/plot/trap),
                            data = svl.merged, family = "bernoulli", 
                            cores = 4)

summary(msel.sex.MeanTSLF.re)
bayes_R2(msel.sex.MeanTSLF.re)

#Diagnostic plots
plot(msel.sex.MeanTSLF.re)
pp_check(msel.sex.MeanTSLF.re)

msel.sex.TSLF.re <- brm(sex ~ TSLF + (1|campaign/plot/trap),
                        data = svl.merged, family = "bernoulli", 
                        cores = 4)

summary(msel.sex.TSLF.re)
bayes_R2(msel.sex.TSLF.re)

#Diagnostic plots
plot(msel.sex.TSLF.re)
pp_check(msel.sex.TSLF.re)

msel.sex.severity.re <- brm(sex ~ severity + (1|campaign/plot/trap),
                            data = svl.merged, family = "bernoulli", 
                            cores = 4)

summary(msel.sex.severity.re)
bayes_R2(msel.sex.severity.re)

#Diagnostic plots
plot(msel.sex.severity.re)
pp_check(msel.sex.severity.re)


msel.sex.firenull.re <- brm(sex ~ 1 + (1|campaign/plot/trap),
                            data = svl.merged, family = "bernoulli", 
                            cores = 4)

summary(msel.sex.firenull.re)
bayes_R2(msel.sex.firenull.re)

#Diagnostic plots
plot(msel.sex.firenull.re)
pp_check(msel.sex.firenull.re)

loo(msel.sex.MeanTSLF.re, msel.sex.TSLF.re, 
    msel.sex.severity.re, msel.sex.firenull.re)


#Conditional effects plot
p1 <- conditional_effects(msel.sex.TSLF.re)

# Not saved!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$TSLF, aes(x = TSLF, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = svl.merged, aes(x = TSLF, y = as.integer(sex)-1), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Time since last fire", y="Proportion of males") +
  theme_minimal()

## Abundance------------------------------------------------------------------------

str(svl.merged)


svl.merged <- svl.merged %>%
  group_by(campaign) %>%
  mutate(sampling_day = dense_rank(date)) %>%  # Rank the dates within each month
  ungroup()

summary(svl.merged$sampling_day)

trap_locations <- fire.regimes.arms %>%
  distinct(plot, trap, treatment)

trap_locations$campaign <- substr(trap_locations$trap,1,2)
trap_locations$campaign[trap_locations$plot=="C1"] <- NA
trap_locations$campaign[trap_locations$plot=="C2"] <- NA
trap_locations$campaign[trap_locations$plot=="QP"] <- NA
trap_locations$campaign[trap_locations$plot=="QT"] <- NA

trap_locations_heitor <- trap_locations[is.na(trap_locations$campaign),]
trap_locations_bruna <- trap_locations[!is.na(trap_locations$campaign),]

trap_locations_heitor <- expand.grid.df(trap_locations_heitor[,-4], 
                                        data.frame(campaign = as.character(c(1:4))))

trap_locations <- rbind(trap_locations_bruna, trap_locations_heitor)

# Now, create the full dataset of all possible sampling events
# `crossing` will combine every row from `trap_locations` with every sampling day.
all_sampling_events <- crossing(trap_locations, sampling_day = 1:15)

# Count captures only for existing combinations in your data
capture_counts <- svl.merged %>%
  count(campaign, plot, trap, sampling_day, name = "freq")

# Join the actual captures to the master list of all sampling events
Ajalapensis.captures.day <- all_sampling_events %>%
  left_join(capture_counts, by = c("campaign", "plot", "trap", "sampling_day")) %>%
  # If a trap on a given day had no captures, its freq is NA. Change that to 0.
  mutate(freq = replace_na(freq, 0))

# View the result
print(Ajalapensis.captures.day)


summary(Ajalapensis.captures.day)
sum(Ajalapensis.captures.day$freq)


Ajalapensis.captures.day <- Ajalapensis.captures.day %>%  
  pivot_wider(names_from = sampling_day, 
              values_from = freq, 
              names_prefix = "day_")

summary(Ajalapensis.captures.day)
str(Ajalapensis.captures.day)

Ajalapensis.captures.day.fire <- left_join(Ajalapensis.captures.day,
                                           fire.regimes.arms[,-c(1,2)],
                                           by = c("plot","trap","treatment"))

Ajalapensis.captures.day.fire <- left_join(Ajalapensis.captures.day.fire,
                                           env.data[, c(1, 2, 5, 18)],
                                           by = c("plot","trap","campaign"))
summary(Ajalapensis.captures.day.fire)

Ajalapensis.captures.day.fire$TSLF.y[is.na(Ajalapensis.captures.day.fire$TSLF.y)] <- Ajalapensis.captures.day.fire$TSLF.x[is.na(Ajalapensis.captures.day.fire$TSLF.y)]*12

Ajalapensis.captures.day.fire$TSLF.y[Ajalapensis.captures.day.fire$treatment=="Q"] <- 0
Ajalapensis.captures.day.fire <- Ajalapensis.captures.day.fire[,-23]
Ajalapensis.captures.day.fire <- dplyr::rename(Ajalapensis.captures.day.fire, "TSLF" = "TSLF.y")

summary(Ajalapensis.captures.day.fire)

# View(Ajalapensis.captures.day.fire)

## brms------------------------------------------------------------------------
Ajalapensis.captures.fire <- Ajalapensis.captures.day.fire[,c(1:4,20:23)]

Ajalapensis.captures.fire$capts <- rowSums(Ajalapensis.captures.day.fire[,c(5:19)])
summary(Ajalapensis.captures.fire)
hist(Ajalapensis.captures.fire$capts)

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
msel.capts.MeanTSLF.re <- brm(capts ~ MeanTSLF + (1|campaign/plot/trap),
                            data = Ajalapensis.captures.fire, family = "poisson", 
                            cores = 4)

summary(msel.capts.MeanTSLF.re)
bayes_R2(msel.capts.MeanTSLF.re)

#Diagnostic plots
plot(msel.capts.MeanTSLF.re)
pp_check(msel.capts.MeanTSLF.re)
plot(conditional_effects(msel.capts.MeanTSLF.re), points = T)

msel.capts.TSLF.re <- brm(capts ~ TSLF + (1|campaign/plot/trap),
                          data = Ajalapensis.captures.fire, family = "poisson", 
                          cores = 4)

summary(msel.capts.TSLF.re)
bayes_R2(msel.capts.TSLF.re)
plot(conditional_effects(msel.capts.TSLF.re), points = T)


#Diagnostic plots
plot(msel.capts.TSLF.re)
pp_check(msel.capts.TSLF.re)

msel.capts.severity.re <- brm(capts ~ severity + (1|campaign/plot/trap),
                              data = Ajalapensis.captures.fire, family = "poisson", 
                              cores = 4, control = list(adapt_delta = 0.9))

summary(msel.capts.severity.re)
bayes_R2(msel.capts.severity.re)

plot(conditional_effects(msel.capts.severity.re), points = T)


#Diagnostic plots
plot(msel.capts.severity.re)
pp_check(msel.capts.severity.re)

msel.capts.firenull.re <- brm(capts ~ 1 + (1|campaign/plot/trap),
                              data = Ajalapensis.captures.fire, family = "poisson", 
                              cores = 4, control = list(adapt_delta = 0.9))

summary(msel.capts.firenull.re)
bayes_R2(msel.capts.firenull.re)

#Diagnostic plots
plot(msel.capts.firenull.re)
pp_check(msel.capts.firenull.re)

loo(msel.capts.MeanTSLF.re, msel.capts.TSLF.re, 
    msel.capts.severity.re, msel.capts.firenull.re)


#Conditional effects plot
p1 <- conditional_effects(msel.capts.severity.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$severity, aes(x = severity, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = Ajalapensis.captures.fire, aes(x = jitter(severity), y = capts), alpha = 0.4, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Mean fire severity", y="Number of captures") +
  theme_minimal()

p1 <- conditional_effects(msel.capts.MeanTSLF.re)

#Salvar!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$MeanTSLF, aes(x = MeanTSLF, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = Ajalapensis.captures.fire, aes(x = MeanTSLF, y = capts), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Mean fire interval", y="Number of captures") +
  theme_minimal()

p1 <- conditional_effects(msel.capts.TSLF.re)

# Not saved!
windows(height = 8, width = 8)
quartz(height = 8, width = 8)
ggplot(p1$TSLF, aes(x = TSLF, y  = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3)+
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = Ajalapensis.captures.fire, aes(x = TSLF, y = capts), alpha = 0.5, size = 2)+
  # scale_colour_manual(values = turbo(4), name = "Fire severity")+
  labs(x="Time since last fire", y="Number of captures")

### INLA------------------------------------------------------------------------

(y.mat <- as.matrix(Ajalapensis.captures.day.fire[,c(5:19)]))
counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    Ajalapensis.captures.day.fire$campaign,
                                    paste(Ajalapensis.captures.day.fire$campaign, 
                                          Ajalapensis.captures.day.fire$plot, sep = "_"),
                                    paste(Ajalapensis.captures.day.fire$campaign, 
                                          Ajalapensis.captures.day.fire$plot, 
                                          Ajalapensis.captures.day.fire$trap, sep = "_"),
                                    scale(Ajalapensis.captures.day.fire$freq), 
                                    scale(Ajalapensis.captures.day.fire$MeanTSLF),
                                    scale(Ajalapensis.captures.day.fire$severity),
                                    scale(Ajalapensis.captures.day.fire$TSLF))

print(counts.and.count.covs)

# Only intercept
out.inla.0 <- inla(counts.and.count.covs[,c(1:16)] ~ 1,
                   data = list(counts.and.count.covs = counts.and.count.covs[,c(1:16)]),
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                                      theta2 = list(prior = "flat",
                                                                    param = numeric()))),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.0, digits = 3)


plot(out.inla.0,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.0,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.0,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = T)


# Only random effects
out.inla.1 <- inla(counts.and.count.covs[,c(1:18)] ~ 1 + f(campaign) + f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs[,c(1:18)],
                               campaign = counts.and.count.covs$X2, 
                               plot_id = counts.and.count.covs$X3),
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                      theta2 = list(param = c(0, 0.1)), 
                                                      theta3 = list(param = c(0, 0.1)),
                                                      theta4 = list(param = c(0, 0.1)),
                                                      theta5 = list(param = c(0, 0.1))
                                                      
                                                      )
                                         ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.1, digits = 3)


plot(out.inla.1,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.1,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 


# Mean fire interval effect

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    scale(Ajalapensis.captures.day.fire$MeanTSLF),
                                    Ajalapensis.captures.day.fire$campaign,
                                    paste(Ajalapensis.captures.day.fire$campaign, 
                                          Ajalapensis.captures.day.fire$plot, sep = "_"))

print(counts.and.count.covs)

out.inla.2 <- inla(counts.and.count.covs ~ 1 + MeanTSLF  + f(campaign) + f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs,
                               campaign = counts.and.count.covs$X3, 
                               plot_id = counts.and.count.covs$X4,
                               MeanTSLF = counts.and.count.covs$X2),
                   
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                      theta2 = list(param = c(0, 0.1)), 
                                                      theta3 = list(param = c(0, 0.1)),
                                                      theta4 = list(param = c(0, 0.1)),
                                                      theta5 = list(param = c(0, 0.1)),
                                                      theta6 = list(param = c(0, 0.1))
                   )
                   ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.2, digits = 3)


plot(out.inla.2,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.2,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE) 

# Mean fire severity effect

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    scale(Ajalapensis.captures.day.fire$severity),
                                    Ajalapensis.captures.day.fire$campaign,
                                    paste(Ajalapensis.captures.day.fire$campaign, 
                                          Ajalapensis.captures.day.fire$plot, sep = "_"))

print(counts.and.count.covs)

out.inla.3 <- inla(counts.and.count.covs ~ 1 + severity  + f(campaign) + f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs,
                               campaign = counts.and.count.covs$X3, 
                               plot_id = counts.and.count.covs$X4,
                               severity = counts.and.count.covs$X2),
                   
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                      theta2 = list(param = c(0, 0.1)), 
                                                      theta3 = list(param = c(0, 0.1)),
                                                      theta4 = list(param = c(0, 0.1)),
                                                      theta5 = list(param = c(0, 0.1)),
                                                      theta6 = list(param = c(0, 0.1))
                   )
                   ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.3, digits = 3)


plot(out.inla.3,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.3,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE) 

# Time since last fire effect

counts.and.count.covs <- inla.mdata(y.mat, 
                                    1,
                                    scale(Ajalapensis.captures.day.fire$TSLF),
                                    Ajalapensis.captures.day.fire$campaign,
                                    paste(Ajalapensis.captures.day.fire$campaign, 
                                          Ajalapensis.captures.day.fire$plot, sep = "_"))

print(counts.and.count.covs)

out.inla.4 <- inla(counts.and.count.covs ~ 1 + TSLF  + f(campaign) + f(plot_id),
                   data = list(counts.and.count.covs = counts.and.count.covs,
                               campaign = counts.and.count.covs$X3, 
                               plot_id = counts.and.count.covs$X4,
                               TSLF = counts.and.count.covs$X2),
                   
                   family = "nmixnb",
                   #Priors for detection parameters
                   control.fixed = list(mean = 0, 
                                        mean.intercept = 0, 
                                        prec = 0.01,
                                        prec.intercept = 0.01),
                   #Priors for abundance and overdispersion parameters
                   control.family = list(hyper = list(theta1 = list(param = c(0, 0.1)),
                                                      theta2 = list(param = c(0, 0.1)), 
                                                      theta3 = list(param = c(0, 0.1)),
                                                      theta4 = list(param = c(0, 0.1)),
                                                      theta5 = list(param = c(0, 0.1)),
                                                      theta6 = list(param = c(0, 0.1))
                   )
                   ),
                   verbose = TRUE,
                   control.inla = list(int.strategy = "eb"),
                   control.compute=list(config = TRUE, 
                                        waic = TRUE,
                                        residuals = TRUE,
                                        cpo = TRUE))
summary(out.inla.4, digits = 3)


plot(out.inla.4,
     plot.fixed.effects = TRUE,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = TRUE,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = TRUE,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = F)

plot(out.inla.4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = TRUE,
     plot.q = F,
     plot.cpo = F) 

plot(out.inla.4,
     plot.fixed.effects = F,
     plot.lincomb = F,
     plot.random.effects = F,
     plot.hyperparameters = F,
     plot.predictor = F,
     plot.q = F,
     plot.cpo = TRUE)

out.inla.0$waic$waic
out.inla.1$waic$waic
out.inla.2$waic$waic
out.inla.3$waic$waic
out.inla.4$waic$waic

inla_models <- list(out.inla.0, 
                    out.inla.1, 
                    out.inla.2, 
                    out.inla.3, 
                    out.inla.4)

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
  Weight = round(w_akaike,2)
)

waic_table$formula[1] <- "(lambda ~ 1), (p ~ 1)"
waic_table$formula[2] <- "(lambda ~ 1), (p ~ (1|campaign/plot))"
waic_table$formula[3] <- "(lambda ~ MeanTSLF), (p ~ MeanTSLF + (1|campaign/plot))"
waic_table$formula[4] <- "(lambda ~ FireSev), (p ~ FireSev + (1|campaign/plot))"
waic_table$formula[5] <- "(lambda ~ TSLF), (p ~ TSLF + (1|campaign/plot))"



# Sort by WAIC (best model first)
waic_table <- waic_table[order(waic_table$WAIC), ]

# Print table
print(waic_table)

summary(out.inla.3)
summary(out.inla.4)

summary(plogis(as.matrix(out.inla.3$summary.linear.predictor[,c("mean", 
                                                                "sd", 
                                                                "0.025quant", 
                                                                "0.5quant",
                                                                "0.975quant",
                                                                "mode")])))
summary(out.inla.3$summary.fitted.values)

inla.nmix.lambda.fitted <- function (result, sample.size = 1000, return.posteriors = FALSE, 
                                     scale = "exp") 
{
  fam <- result$.args$family
  if (length(grep(pattern = "nmix", x = fam)) == 0) {
    stop("This function is only for models with 'nmix' or 'nmixnb' likelihoods")
  }
  if (missing(result)) {
    stop("Please specify a model result")
  }
  s.check <- as.numeric(scale == "exp") + as.numeric(scale == 
                                                       "log")
  if (s.check == 0) {
    stop("Scale must be set to 'exp' or 'log'")
  }
  if (sample.size < 500) {
    warning("Please increase the sample size")
  }
  mdata.obj <- result$.args$data[[1]]
  counts <- as.data.frame(mdata.obj[grep(pattern = "Y", names(mdata.obj))])
  lambda.covs <- as.data.frame(mdata.obj[grep(pattern = "X", 
                                              names(mdata.obj))])
  lambda.covs <- as.matrix(lambda.covs[,-c(3,4)])
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
  fitted.meds <- round(apply(fitted.posteriors, 1, median), 
                       4)
  fitted.means <- round(apply(fitted.posteriors, 1, mean), 
                        4)
  fitted.sds <- round(apply(fitted.posteriors, 1, sd), 4)
  fitted.q025 <- round(apply(fitted.posteriors, 1, quantile, 
                             probs = 0.025), 4)
  fitted.q500 <- round(apply(fitted.posteriors, 1, quantile, 
                             probs = 0.5), 4)
  fitted.q975 <- round(apply(fitted.posteriors, 1, quantile, 
                             probs = 0.975), 4)
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  fitted.modes <- round(apply(fitted.posteriors, 1, Mode), 
                        4)
  fitted.summary <- data.frame(mean.lambda = fitted.means, 
                               sd.lambda = fitted.sds, quant025.lambda = fitted.q025, 
                               median.lambda = fitted.q500, quant975.lambda = fitted.q975, 
                               mode.lambda = fitted.modes)
  fitted.summary <- cbind(index, fitted.summary)
  fitted.posteriors <- cbind(index, fitted.posteriors)
  if (return.posteriors == TRUE) {
    out <- list(fitted.summary = fitted.summary, fitted.posteriors = fitted.posteriors)
  }
  else {
    out <- list(fitted.summary = fitted.summary)
  }
  return(out)
}


# Severity effect
# Frequency effects
out.inla.3.lambda.fits <- inla.nmix.lambda.fitted(result = out.inla.3,
                                                  sample.size = 10000, 
                                                  return.posteriors = FALSE)$fitted.summary

ggplot(stack(out.inla.3.lambda.fits[,c("mean.lambda", "median.lambda", "mode.lambda")]), 
       aes(x = values, fill = ind)) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "Abundance")+
  theme_classic()

ggplot(out.inla.3.lambda.fits, aes(x = median.lambda)) +
  geom_density(fill = "lightblue", alpha = 0.7, color = "black") +
  labs(
    y = "Density", 
    x = "Abundance"
  ) +
  theme_classic()

ggplot(stack(log(out.inla.3.lambda.fits[,c("mean.lambda", "median.lambda", "mode.lambda")])), 
       aes(x = values, fill = ind)) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "log(Abundance)")+
  theme_classic()

severity_xx <- seq(min(scale(Ajalapensis.captures.day.fire$severity)), 
                   max(scale(Ajalapensis.captures.day.fire$severity)), 
                   by = 0.1)

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

# Step 3: Run the corrected sapply loop
predicted_lines <- sapply(posterior_samples, function(s) {
  # Access the column by its position [ ,1] instead of by name
  intercept <- s$latent[intercept_idx, 1]
  slope <- s$latent[slope_idx, 1]
  
  # The rest of the code is the same
  linear_predictor <- intercept + slope * severity_xx
  return(plogis(linear_predictor))
})

rel_p_severity <- data.frame(severity = severity_xx,
                         median = apply(predicted_lines, 1, quantile, 0.5),
                         lower = apply(predicted_lines, 1, quantile, 0.025),
                         upper = apply(predicted_lines, 1, quantile, 0.975))

pred_df <- data.frame(severity = scale(Ajalapensis.captures.day.fire$severity),
                      trap = Ajalapensis.captures.day.fire$trap,
                      N = out.inla.4.lambda.fits$median.lambda,
                      Nlow = out.inla.4.lambda.fits$quant025.lambda,
                      Nhigh = out.inla.4.lambda.fits$quant975.lambda,
                      p = plogis(out.inla.4$summary.linear.predictor$`0.5quant`),
                      ncaps = rowSums(y.mat))

quartz(h = 8, w = 8)
ggplot(pred_df, aes(x = severity, y = p)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_line(data = rel_p_severity, aes(x = severity, y = median), col = "blue", linewidth = 2) +
  geom_ribbon(data = rel_p_severity, aes(ymin = lower, ymax = upper), alpha = 0.3) +
  labs(x = "Mean fire severity", y = "Capture probability") +
  scale_y_log10(labels = number)+
  theme_minimal()

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
monthly_freq_stack <- tapp(fire_binary_stack, index = raster_months, fun = "sum", na.rm = TRUE)

plot(monthly_freq_stack)

#--- Step 2: Apply Seasonal Weights ---
# Define the weights for each of the 12 months
seasonal_weights <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3)

# Multiply the 12-layer frequency stack by the 12 weights
monthly_severity_stack <- monthly_freq_stack * seasonal_weights

plot(monthly_severity_stack)
#--- Step 3: Calculate the Final Mean Severity ---
# As in your script, calculate the mean of the 12 monthly severity scores for each pixel
severity_raster <- app(monthly_severity_stack, "mean", na.rm = TRUE)
names(severity_raster) <- "severity"

# Optional: Save your new raster
writeRaster(severity_raster, "SGT_Severity.tif", overwrite = TRUE)

# Plot to check
plot(severity_raster)

#Read and plot points coordinates
arms.pts <- read.table("Pontos_Arms_Bruna_Heitor.txt",h=T)
summary(arms.pts)
arms.pts$sampdes <- c(rep("Rotating", 172),
                      rep("Fixed", 48))

coordinates(arms.pts) <- c("Longitude","Latitude")
proj4string(arms.pts) <- CRS("+proj=longlat +datum=WGS84")
arms.pts
arms.pts <- st_as_sf(arms.pts)

plot(severity_raster)
plot(arms.pts, add = T)

#Read biomes
biomes <- read_biomes()
biomes

Cerrado <- biomes[biomes$name_biome=="Cerrado",]

#Read Conservation Units
CUs <- read_conservation_units()

SGT <- CUs[grep("SERRA GERAL DO TOCANTINS", CUs$name_conservation_unit),]

#Read states
states <- read_state()

#Read South America
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
  geom_rect(aes(xmin = -47.3, xmax = -45.8, ymin = -11.5, ymax = -10.4), color = "red", fill = NA) +
  theme_test()  +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "lightblue")) +
  coord_sf(xlim = c(-74, -35), ylim = c(-34, 6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

inset.sgt <- ggplot() +
  geom_sf(data = SGT, color = "forestgreen", fill = NA, linewidth = 1.25) +
  # Add your trap locations
  geom_sf(
    data = arms.pts, 
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("Fixed" = 22, "Rotating" = 21), name = "Sampling Design") +
  geom_rect(aes(xmin = -47.3, xmax = -45.8, ymin = -11.5, ymax = -10.4), color = "red", fill = NA) +
  theme_test()  +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "lightgray")) +
  coord_sf(xlim = c(-47.3, -45.8), ylim = c(-11.5, -10.4), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

# dowload tiles and compose raster (SpatRaster)
bbox.arms <- st_bbox(arms.pts)
bbox.arms[1:4] <- c(bbox.arms[1] - 0.05, 
                    bbox.arms[2] - 0.01, 
                    bbox.arms[3] + 0.05, 
                    bbox.arms[4] + 0.01)

imagery <- get_tiles(bbox.arms, crop = TRUE, provider = "Esri.WorldImagery")
plot(imagery)
plot(st_geometry(arms.pts), col = "red", pch = 21, add=T)

map.arms.sgt <- ggplot() +
  geom_spatraster_rgb(data = imagery, interpolate = F) +
  # Add your trap locations
  geom_sf(
    data = arms.pts, 
    size = 3, 
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3
  ) +
  scale_shape_manual(values = c("Fixed" = 22, "Rotating" = 21), name = "Sampling Design") +
  xlim(bbox.arms[1], bbox.arms[3]) +
  ylim(bbox.arms[2], bbox.arms[4]) +
  theme_test() +
  guides(size = "none") +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white"),
    pad_y = unit(0.2, "in")
  ) +
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0.53, "in"), pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20")
  )+
  coord_sf(xlim = c(bbox.arms[1], bbox.arms[3]), ylim = c(bbox.arms[2], bbox.arms[4]), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")


quartz(width = 12, height = 10)
map.arms.sgt

# Combining both maps
print(inset.cerrado, vp = viewport(0.7, 0.8, width = 0.2, height = 0.2))
print(inset.sgt, vp = viewport(0.68, 0.6, width = 0.2, height = 0.2))

# Map with fire severity

#--- Create the Final Map ---
map.arms.severity <- ggplot() +
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
    data = arms.pts, 
    size = 3, 
    aes(shape = sampdes), # Use shape to distinguish sampling design
    fill = "white",
    color = "black",
    alpha = 0.3
  ) +
  scale_shape_manual(values = c("Fixed" = 22, "Rotating" = 21), name = "Sampling Design") +
  
  # Add scale bar and north arrow
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white")
  ) +
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0.53, "in"), pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20")
  )+
  
  # Set coordinates and labels
  coord_sf(xlim = c(bbox.arms[1], bbox.arms[3]), ylim = c(bbox.arms[2], bbox.arms[4]), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()

# Display the map
print(map.arms.severity)

quartz(width = 12, height = 10)
print(map.arms.severity)

# Combining both maps
print(inset.cerrado, vp = viewport(0.7, 0.8, width = 0.2, height = 0.2))
print(inset.sgt, vp = viewport(0.68, 0.6, width = 0.2, height = 0.2))
