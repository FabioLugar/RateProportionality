setwd("~/ProportionalityIII/")
library(MCMCglmm)
library(treeplyr)
library(plyr)
library(dplyr)
library(evolqg)
library(ggplot2)
library(reshape)
library(psych)
set.seed(293923)
{traitnames <- c("m1.bl", "m2.bl", "m3.bl", "m1.md",  "m2.md", "m3.md")
  GcorR <- GcorL <- matrix(NA, nrow=6, ncol=6)
  diag(GcorR) <- diag(GcorL) <- 1
  rownames(GcorR) <- rownames(GcorL) <- colnames(GcorL) <- colnames(GcorR) <- traitnames
  GcorR["m1.md", "m1.bl"] <- GcorR["m1.bl", "m1.md"]  <- 0.403
  GcorR["m2.md", "m2.bl"] <- GcorR["m2.bl", "m2.md"]  <- 0.706
  GcorR["m3.md", "m3.bl"] <- GcorR["m3.bl", "m3.md"]  <- 0.654
  GcorR["m1.md", "m2.md"] <- GcorR["m2.md", "m1.md"]  <- 0.876
  GcorR["m1.md", "m3.md"] <- GcorR["m3.md", "m1.md"]  <- 0.744
  GcorR["m2.md", "m3.md"] <- GcorR["m3.md", "m2.md"]  <- 0.920
  GcorR["m1.bl", "m2.bl"] <- GcorR["m2.bl", "m1.bl"]  <- 0.903
  GcorR["m1.bl", "m3.bl"] <- GcorR["m3.bl", "m1.bl"]  <- 0.583
  GcorR["m2.bl", "m3.bl"] <- GcorR["m3.bl", "m2.bl"]  <- 0.850
  GcorR["m1.md", "m2.bl"] <- GcorR["m2.bl", "m1.md"]  <- 0.525
  GcorR["m1.md", "m3.bl"] <- GcorR["m3.bl", "m1.md"]  <- 0.358
  GcorR["m2.md", "m3.bl"] <- GcorR["m3.bl", "m2.md"]  <- 0.590
  GcorR["m1.bl", "m2.md"] <- GcorR["m2.md", "m1.bl"]  <- 0.478
  GcorR["m1.bl", "m3.md"] <- GcorR["m3.md", "m1.bl"]  <- 0.227
  GcorR["m2.bl", "m3.md"] <- GcorR["m3.md", "m2.bl"]  <- 0.601 #missing from the dataset! This is the left value.
  
  GcorL["m1.md", "m1.bl"] <- GcorL["m1.bl", "m1.md"]  <- 0.484
  GcorL["m2.md", "m2.bl"] <- GcorL["m2.bl", "m2.md"]  <- 0.491
  GcorL["m3.md", "m3.bl"] <- GcorL["m3.bl", "m3.md"]  <- 0.738
  GcorL["m1.md", "m2.md"] <- GcorL["m2.md", "m1.md"]  <- 0.818
  GcorL["m1.md", "m3.md"] <- GcorL["m3.md", "m1.md"]  <- 0.777
  GcorL["m2.md", "m3.md"] <- GcorL["m3.md", "m2.md"]  <- 0.902
  GcorL["m1.bl", "m2.bl"] <- GcorL["m2.bl", "m1.bl"]  <- 0.899
  GcorL["m1.bl", "m3.bl"] <- GcorL["m3.bl", "m1.bl"]  <- 0.790
  GcorL["m2.bl", "m3.bl"] <- GcorL["m3.bl", "m2.bl"]  <- 0.876
  GcorL["m1.md", "m2.bl"] <- GcorL["m2.bl", "m1.md"]  <- 0.319
  GcorL["m1.md", "m3.bl"] <- GcorL["m3.bl", "m1.md"]  <- 0.436
  GcorL["m2.md", "m3.bl"] <- GcorL["m3.bl", "m2.md"]  <- 0.494
  GcorL["m1.bl", "m2.md"] <- GcorL["m2.md", "m1.bl"]  <- 0.199
  GcorL["m1.bl", "m3.md"] <- GcorL["m3.md", "m1.bl"]  <- 0.663
  GcorL["m2.bl", "m3.md"] <- GcorL["m3.md", "m2.bl"]  <- 0.761 #missing from the dataset! This is the left value.
  
  Gcor <- (GcorL+GcorR)/2
  eigen(Gcor)
}

# importing data----
Gdata <- read.csv("./hluskoG.csv")
editnames <- read.csv("./EditedNames_key.csv",as.is=T)
tree<-read.nexus("./tree_2.nex")
dat<-read.table("./Plavcan_Lower Molar_All.txt",header=T,sep="\t",as.is = T)
dat$sex[dat$sex=="?"]<-NA
dat$sex[dat$sex==""]<-NA
for(i in 1:nrow(editnames)){
  dat$Species[dat$Species==editnames[i,1]] <- editnames[i, 2]
}
dat<-dplyr::select(dat, Species:m3.md, m1.bl.tal, m2.bl.tal, m3.bl.tal)
colnames(dat)<-sub(".tal","",colnames(dat))

# logtransforming----
dat[,4:9] <- log(dat[,4:9])

# Pooled within-group P and SDs for scaling----
Plik<-evolqg::CalculateMatrix(lm(cbind(m1.md,m2.md,m3.md,m1.bl,m2.bl,m3.bl)~sex*Species, data = dat))
sds<-sqrt(diag(Plik))

dat_scaled<-cbind(apply(dat[,c(1,3)],2,factor),sweep(dat[,-c(1:3)],2,sds,FUN = "/"))
dat_scaled<-subset(dat_scaled,!is.na(sex))

#running glmm with babbon G as a strong prior----
# prior <- list(R=list(V=nearPD(Gcor)$mat, nu=1000), G=list(G1=list(V=diag(6), nu=0.002)))
prior <- list(R=list(V=diag(6), nu=0.002), G=list(G1=list(V=diag(6), nu=0.002)))

mglm_out <- MCMCglmm(cbind(m1.md,m2.md,m3.md,m1.bl,m2.bl,m3.bl) ~ trait + trait:sex - 1, 
                     random=~us(trait):Species, rcov= ~us(trait):units, 
                     prior=prior, data=dat_scaled, 
                     family=c("gaussian", "gaussian", "gaussian",
                              "gaussian", "gaussian", "gaussian"),
                     nitt=26000, thin=20, burnin=6000)

# generating heritabilities ----
h2s<-adply(1:6,1,function(i){
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  
  betaPars<-estBetaParams(mu = mean(Gdata$resid.h2[i],Gdata$resid.h2[i+6]),
                          var= mean(Gdata$resid.h2SE[i],Gdata$resid.h2SE[i+6])^2)
  rbeta(1000,betaPars$alpha,betaPars$beta)
})[,-1]

Gs <- laply(1:1000,
            function(i) {
              cov2cor(ExtendMatrix(Gcor,ret.dim = 5)$ExtMat) * tcrossprod(sqrt(h2s[,i])*sds)
            }) %>%
  aperm(., c(2,3,1))

Ps <- laply(1:1000,
            function(i) {
              M<-matrix(mglm_out$VCV[i,37:72], ncol=6, nrow=6) 
              M * tcrossprod(sds)
            }) %>%
  aperm(., c(2,3,1))

# write.csv(as.matrix(laply(1:1000, function(i) KrzCor(cov2cor(Ps[,,i]),Gcor)) %>% summary), file="results/Krz_cor.csv")

propCors.plot<-
  laply(1:1000, function(i) c(cov2cor(Ps[,,i]))) %>% t %>% data.frame(.,G=c(Gcor)) %>% melt(., id.vars="G") %>% unique %>% subset(., G<1) %>%
  ggplot(., aes(G,value))+
  geom_abline(aes(slope=1, intercept=0))+
  stat_summary(fun = median, fun.min = min, fun.max = max)+
  geom_smooth(method="lm", color="black", linetype=2)+
  xlab("Genetic correlations")+
  ylab("Phenotypic correlations")
ggsave("results/propCors.pdf",propCors.plot,width = 5,height = 5)

colnames(Ps)<-rownames(Ps)<-colnames(Gs)<-rownames(Gs)<-traitnames

spdata <- dat %>% group_by(Species) %>% 
  dplyr::summarize(SE.m1.md=sd(m1.md,na.rm=TRUE)/sqrt(n()),
                   SE.m2.md=sd(m2.md,na.rm=TRUE)/sqrt(n()),
                   SE.m3.md=sd(m3.md,na.rm=TRUE)/sqrt(n()),
                   SE.m1.bl=sd(m1.bl,na.rm=TRUE)/sqrt(n()),
                   SE.m2.bl=sd(m2.bl,na.rm=TRUE)/sqrt(n()),
                   SE.m3.bl=sd(m3.bl,na.rm=TRUE)/sqrt(n()),
                   m1.md=mean(m1.md,na.rm=TRUE),
                   m2.md=mean(m2.md,na.rm=TRUE),
                   m3.md=mean(m3.md,na.rm=TRUE),
                   m1.bl=mean(m1.bl,na.rm=TRUE),
                   m2.bl=mean(m2.bl,na.rm=TRUE),
                   m3.bl=mean(m3.bl,na.rm=TRUE))

td <- make.treedata(tree, spdata)

save.image("data_imported.Rdata")
