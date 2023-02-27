setwd("~/ProportionalityIII/")
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
registerDoParallel(cores = 25)
options(PCMBase.Threshold.EV = 1e-8)
source('parametersLimits.R', local=FALSE)
load("data_imported.Rdata")

#importing trees and data----

treeBlock<-read.nexus("treeBlock.nex")
treeBlock<-drop.tip.multiPhylo(treeBlock,treeBlock[[1]]$tip.label[!treeBlock[[1]]$tip.label %in% td$phy$tip.label])
for(i in 1:length(treeBlock)) {
  treeBlock[[i]]<-di2multi(treeBlock[[i]],tol = 0.05)
  treeBlock[[i]]<-force.ultrametric(treeBlock[[i]])
  treeBlock[[i]]<-PCMTree(treeBlock[[i]])
}

X<-t(td$dat %>% select(., m1.md:m3.bl))
SEsp<-t(td$dat %>% select(., SE.m1.md:SE.m3.bl)) 
SEsp[is.na(SEsp)]<-0
colnames(X)<- colnames(SEsp)<- treeBlock[[1]]$tip.label

# Runing a constrasts and ancs for initial values----
X0<-aaply(X,1,function(x){
  td <- make.treedata(tree, c(na.omit(x)))
  fastAnc(td$phy,td$dat$trait)[1]
})

Dpic<-aaply(X[,!is.na(X[3,])],1,function(x){
  td <- make.treedata(tree, x)
  pic(td$dat$trait,multi2di(td$phy))
}) %>% tcrossprod

#Seting model names----
BM_name<-paste0("BM",
                "__Global_X0",
                "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x",
                "__Omitted_Sigmae_x")

BMp_name<-paste0("BMkappa",
                 "__Global_X0",
                 "__NonNegative_kappa",
                 "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Fixed_Sigma_x",
                 "__Omitted_Sigmae_x")

modelFits<-list(BM=vector("list",length(treeBlock)),
                BMkappaP=vector("list",length(treeBlock)),
                BMkappaG=vector("list",length(treeBlock)))


for(i in 1:length(treeBlock)){
  print(paste0("Iter=",i))
  ##############
  print("Fitting full BM")
  
  model<-PCM(BM_name, k = 6,
             params = list(X0=X0,
                           Sigma_x = abind(UpperTriFactor(Dpic),along=3)))
  modelFits$BM[[i]]<-
    PCMFit(model = model, tree = treeBlock[[i]], X = X,
           metaI = PCMInfoCpp,
           SE = SEsp, doParallel =T)
  print(paste0("loglik=",round(modelFits$BM[[i]]$logLikOptim,3)))
  print("Fitting kappa P")
  X0<-modelFits$BM[[i]]$modelOptim$X0
  kappa<-tr(tcrossprod(modelFits$BM[[i]]$modelOptim$Sigma_x[,,1]))/tr(Ps[,,i])
  model<-PCM(BMp_name, k = 6,
             params = list(X0=X0,
                           kappa = kappa,
                           Sigma_x = abind(UpperTriFactor(Ps[,,i]),along=3)))
  
  modelFits$BMkappaP[[i]]<-
    PCMFit(model = model, tree = treeBlock[[i]], X = X,
           metaI = PCMInfoCpp,
           SE = SEsp, doParallel =T)
  print(paste0("loglik=",round(modelFits$BMkappaP[[i]]$logLikOptim,3)))
  print("Fitting kappa G")
  
  model<-PCM(BMp_name, k = 6,
             params = list(X0=modelFits$BM[[i]]$modelOptim$X0,
                           kappa = kappa,
                           Sigma_x = abind(UpperTriFactor(Gs[,,i]),along=3)))
  
  modelFits$BMkappaG[[i]]<-
    PCMFit(model = model, tree = treeBlock[[i]], X = X,
           metaI = PCMInfoCpp,
           SE = SEsp, doParallel =T)
  print(paste0("loglik=",round(modelFits$BMkappaG[[i]]$logLikOptim,3)))
  if(i %in% seq(50,1000, by=50)) {
    save.image("fitDist.Rdata")
    print("Saved a copy")}
  print(paste0("Iter=",i," DONE"))
}

save.image("fitDist.Rdata")