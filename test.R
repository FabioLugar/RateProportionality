setwd("~/ProportionalityII/")
library(treeplyr)
library(plyr)
library(dplyr)
library(evolqg)
library(ggplot2)
library(ggtree)
library(reshape)
library(PCMFit)
library(PCMkappa)
library(PCMBaseCpp)
library(abind)
library(cowplot)
library(doParallel)
library(phytools)
library(MuMIn)
library(psych)
registerDoParallel(cores = 50)
options(PCMBase.Threshold.EV = 1e-8)
source('parametersLimits.R', local=FALSE)
load("data_imported.Rdata")

X<-t(td$dat %>% select(., m1.md:m3.bl))
SEsp<-t(td$dat %>% select(., SE.m1.md:SE.m3.bl)) 
SEsp[is.na(SEsp)]<-0
colnames(X)<- colnames(SEsp)<- td$phy$tip.label

X0<-aaply(X,1,function(x){
  td <- make.treedata(tree, c(na.omit(x)))
  fastAnc(td$phy,td$dat$trait)[1]
})

Dpic<-aaply(X[,!is.na(X[3,])],1,function(x){
  td <- make.treedata(tree, x)
  pic(td$dat$trait,multi2di(td$phy))
}) %>% tcrossprod

kappa=1

BM_name<-paste0("BM",
                "__Global_X0",
                "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x",
                "__Omitted_Sigmae_x")

BMp_name<-paste0("BMkappa",
                 "__Global_X0",
                 "__NonNegative_kappa",
                 "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Fixed_Sigma_x",
                 "__Omitted_Sigmae_x")

modelBM<-PCM(BM_name, k = 6,
             params = list(X0=X0,
                           Sigma_x = abind(UpperTriFactor(Dpic),along=3)))
fitBM<-
  PCMFit(model = modelBM, tree = td$phy, X = X,
         metaI = PCMInfoCpp,
         SE = SEsp, doParallel =T)

X0<-fitBM$modelOptim$X0
kappa<-tr(tcrossprod(fitBM$modelOptim$Sigma_x[,,1]))/tr(Ps[,,1])
modelBMkappa<-PCM(BMp_name, k = 6,
                  params = list(X0=X0,
                                kappa = kappa,
                                Sigma_x = abind(UpperTriFactor(Ps[,,1]),along=3)))

fitBMkappaP<-
  PCMFit(model = modelBMkappa, tree = td$phy, X = X,
         metaI = PCMInfoCpp,
         SE = SEsp, doParallel =T)


PCMLik(model = modelBMkappa, tree = td$phy, X = X,SE = SEsp)


########sim

sim<-PCMSim(tree=td$phy, model=fitBMkappaP$modelOptim, X0=X0, SE=SEsp)

fitBM<-
  PCMFit(model = modelBM, tree = td$phy, X = sim[,td$phy$tip.label],
         metaI = PCMInfoCpp,
         SE = SEsp, doParallel =T)
fitBMkappaP<-
  PCMFit(model = modelBMkappa, tree = td$phy, X = sim[,td$phy$tip.label],
         metaI = PCMInfoCpp,
         SE = SEsp, doParallel =T)

fittedBM<-fitBM$modelOptim
fittedBMkappa<-PCMApplyTransformation(fitBMkappaP$modelOptim)

SigmaFull  <- tcrossprod(fittedBM$Sigma_x[,,1]) 
Sigmakappa <- tcrossprod(fittedBMkappa$Sigma_x[,,1])
P<-Ps[,,1]

tr(SigmaFull)/tr(P)
tr(Sigmakappa)/tr(P)
tr(P)

