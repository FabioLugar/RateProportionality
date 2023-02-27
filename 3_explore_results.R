setwd("~/ProportionalityIII/")
load("fitDist.Rdata")
library(plyr)
library(dplyr)
library(ggplot2)
library(doParallel)
library(reshape)
library(ggbeeswarm)
library(evolqg)
library(PCMBase)
library(PCMFit)
library(PCMBaseCpp)
library(PCMkappa)
library(cowplot)
library(emorph2)
library(geomorph)
library(treeplyr)
registerDoParallel(cores = 50)
theme_set(theme_bw())

spmean<-t(X)
spse<-t(SEsp)
for(j in seq_along(spse[1,])){
  spse[spse[,j]==0 & !is.na(spmean[,j]),j]<-max(spse[,j])
}
spse<-spse[!is.na(spmean[,3]),]
spmean<-spmean[!is.na(spmean[,3]),]

physignalK<-ldply(seq_along(treeBlock), function(i) {
  spmean.x<-ldply(seq_along(spmean[,1]), function(j){
    rnorm(6,spmean[j,],spse[j,])
  })
  rownames(spmean.x)<-rownames(spmean)
  drop<-treeBlock[[i]]$tip.label[!treeBlock[[i]]$tip.label %in% rownames(spmean.x)]
  tr<- drop.tip(treeBlock[[i]],drop)
  class(tr)<-class(tr)[-1]
  p<-physignal(as.matrix(spmean.x), tr,iter = 1)
  data.frame(obs=p$phy.signal,rand=p$random.K[2])
},.parallel = T)

kmult.plot<-
  physignalK %>% subset(., obs>0.00001) %>%
  ggplot(., aes(obs))+
  geom_density(fill="gray",adjust = 2)+
  geom_vline(aes(xintercept=quantile(rand, 0.95)), linetype=2)+
  xlab("K-mult")
  
ggsave("results/kmult.pdf",
       kmult.plot,
       width = 5,
       height = 4)


i <- sum(!laply(modelFits$BM, is.null))
likBMmodel <- laply(modelFits$BM[1:i], function(X)
  X$logLikOptim[1])
likBMkappaPmodel <-
  laply(modelFits$BMkappaP[1:i], function(X)
    X$logLikOptim[1])
likBMkappaGmodel <-
  laply(modelFits$BMkappaG[1:i], function(X)
    X$logLikOptim[1])
n <- dim(X)[2]

modelCompare <- rbind(
  data.frame(
    model = "BM",
    logLik = likBMmodel,
    k = PCMParamCount(modelFits$BM[[1]]$modelOptim)
  ),
  data.frame(
    model = "BMkappaP",
    logLik = likBMkappaPmodel,
    k = PCMParamCount(modelFits$BMkappaP[[1]]$modelOptim)
  ),
  data.frame(
    model = "BMkappaG",
    logLik = likBMkappaGmodel,
    k = PCMParamCount(modelFits$BMkappaG[[1]]$modelOptim)
  )
) %>%
  mutate(
    .,
    AICc = (2 * k - 2 * logLik) + (2 * k) * (k + 1) / (n - k - 1),
    BIC = (log(n) * k - 2 * logLik)
  )

# I <- likBMmodel>0 & likBMkappaGmodel>0 & likBMkappaPmodel>0
# I <- likBMmodel>500 & likBMkappaGmodel>500 & likBMkappaPmodel>500
I <-
  likBMmodel > likBMkappaPmodel &
  likBMmodel > likBMkappaGmodel & 
  likBMkappaGmodel > 570

d2p <- modelCompare %>% select(., -k) %>%
  subset(., c(I, I, I)) %>%
  melt
CPCOLS <- c("#87CEEB", "#EEE685", "#436EEE")
p1 <-
  ggplot(d2p, aes(model, value)) +
  facet_wrap(variable~. , scales = "free_y") +
  geom_violin(scale = 3,
              alpha = 1,
              show.legend = FALSE) +
  geom_quasirandom(
    aes(color = model),
    alpha = 1,
    size = 0.1,
    show.legend = F
  ) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = CPCOLS) +
  # coord_flip()+
  scale_x_discrete(labels = c("BM", expression(Sigma %prop% P), expression(Sigma %prop% G)))
p1
ggsave("results/model_compare.pdf",
       p1,
       width = 10,
       height = 4)

{
  modeldifs <-
    rbind(
      data.frame(
        model = "P",
        loglik = modelCompare$logLik[1001:2000] - modelCompare$logLik[1:1000],
        AICc = modelCompare$AICc[1001:2000] - modelCompare$AICc[1:1000],
        BIC = modelCompare$BIC[1001:2000] - modelCompare$BIC[1:1000]
      )[I, ],
      data.frame(
        model = "G",
        loglik = modelCompare$logLik[2001:3000] - modelCompare$logLik[1:1000],
        AICc = modelCompare$AICc[2001:3000] - modelCompare$AICc[1:1000],
        BIC = modelCompare$BIC[2001:3000] - modelCompare$BIC[1:1000]
      )[I, ]
    ) #%>% subset(., abs(loglik) < 5)
  
  dens.temp <- apply(subset(modeldifs, model == "P")[, -1], 2, density)
  dens.temp <- ldply(dens.temp, function(x)
    data.frame(x = x$x, y = x$y))
  names(dens.temp)[1] <- "variable"
  dens.temp$model <- "P"
  dens <- dens.temp
  dens.temp <- apply(subset(modeldifs, model == "G")[, -1], 2, density)
  dens.temp <- ldply(dens.temp, function(x)
    data.frame(x = x$x, y = x$y))
  names(dens.temp)[1] <- "variable"
  dens.temp$model <- "G"
  dens <- rbind(dens, dens.temp)
  
  dens.temp <- apply(subset(modeldifs, model == "P")[, -1], 2, density)
  dens.temp <- ldply(dens.temp, function(x) {
    ql <- quantile(x$x, 0.025)
    qu <- quantile(x$x, 0.975)
    i <- x$x < qu & x$x > ql
    data.frame(x = x$x[i], y = x$y[i])
  })
  names(dens.temp)[1] <- "variable"
  dens.temp$model <- "P"
  densCI <- dens.temp
  dens.temp <- apply(subset(modeldifs, model == "G")[, -1], 2, density)
  dens.temp <- ldply(dens.temp, function(x) {
    ql <- quantile(x$x, 0.025)
    qu <- quantile(x$x, 0.975)
    i <- x$x < qu & x$x > ql
    data.frame(x = x$x[i], y = x$y[i])
  })
  names(dens.temp)[1] <- "variable"
  dens.temp$model <- "G"
  densCI <- rbind(densCI, dens.temp)
  
  dens$model <- factor(dens$model)
  levels(dens$model) <-
    c(expression(Sigma %prop% G), expression(Sigma %prop% P))
  densCI$model <- factor(densCI$model)
  levels(densCI$model) <-
    c(expression(Sigma %prop% G), expression(Sigma %prop% P))
}

# CPCOLS <- c("#FFF68F", "#87CEEB")
CPCOLS <- c("#EEE685", "#436EEE")

ggplot(modeldifs, aes(model, loglik)) +
  geom_jitter()


p2.1 <-
  dens %>% subset(., variable == "loglik") %>%
  ggplot(., aes(x, y)) +
  facet_grid(variable ~ model, scales = "free", labeller = label_parsed) +
  geom_area(aes(x = x, y = y, fill = model),
            densCI %>% subset(., variable == "loglik"),
            show.legend = F) +
  geom_line(aes(x = x, y = y), dens %>% subset(., variable == "loglik")) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_bw() +
  scale_fill_manual(values = CPCOLS) +
  xlab("BM-BMkappa") +
  ylab("")
# annotate("label",y=0.016,x=-100,label="BM")+
# annotate("label",y=0.016,x=100,label="BMkappa")
p2.2 <-
  dens %>% subset(., variable == "BIC") %>%
  ggplot(., aes(-x, y)) +
  facet_grid(variable ~ model, scales = "free", labeller = label_parsed) +
  geom_area(aes(x = -x, y = y, fill = model),
            densCI %>% subset(., variable == "BIC"),
            show.legend = F) +
  geom_line(aes(x = -x, y = y), dens %>% subset(., variable == "BIC")) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_bw() +
  scale_fill_manual(values = CPCOLS) +
  xlab("BMkappa-BM") +
  ylab("")
# xlim(-300,300)+
# annotate("label",y=0.0075,x=-230,label="BM")+
# annotate("label",y=0.0075,x=230,label="BMkappa")

p2 <- plot_grid(p2.1, p2.2, nrow = 2, align = "v")
p2

ggsave("results/model_compare_difs.pdf",
       p2,
       width = 7,
       height = 7)

sigma_dist <- modelFits$BM[1:i][I] %>%
  llply(., function(x) {
    tcrossprod(x$modelOptim$Sigma_x[, , 1])
  })
names(sigma_dist) <- paste0("Sigma", 1:length(sigma_dist))
P_dist <- alply(Ps[, , I], 3, function(x)
  x)
names(P_dist) <- paste0("P", 1:sum(I))
G_dist <- alply(Gs[, , I], 3, function(x)
  x)
names(G_dist) <- paste0("G", 1:sum(I))
allMats <- c(G_dist[sample(1:sum(I),100)], P_dist[sample(1:sum(I),100)], sigma_dist[sample(1:sum(I),100)])

simFull <- PCAsimilarity(allMats, parallel = T)

p3 <-
  rbind(
    # data.frame(
    #   comparison = "within Sigmas",
    #   similarity = simFull[201:300, 201:300] %>% c 
    # ),
    data.frame(
      comparison = "with Gs",
      similarity = simFull[201:300, 1:100] %>% c
    ),
    data.frame(
      comparison = "with Ps",
      similarity = simFull[201:300, 101:200] %>% c 
    )
  ) %>%
  subset(., similarity != 0) %>%
  ggplot(., aes(comparison, similarity)) +
  geom_quasirandom(alpha = 1, size = 0.1) +
  geom_violin(scale = 3, alpha = 0.5) +
  xlab("") +
  ylab("PCA similarity") +
  ylim(c(0.9, 1.0))
p3
ggsave("results/matrix_sim.pdf",
       p3,
       width = 5,
       height = 5)

i=1


vector_corr.plot<-foreach(i=1:100, .combine = "rbind") %do% {
  abs(diag(t(eigen(sigma_dist[[i]])$vector) %*% eigen(P_dist[[i]])$vector))
} %>% melt %>%
  ggplot(., aes(as.factor(X2), value))+
  stat_summary(fun = median,
               geom = "pointrange",
               fun.min = min,
               fun.max = max)+
  xlab("Eigenvector Rank")+
  ylab("Vector Correlation")
ggsave("results/vector_corr.pdf",
       vector_corr.plot,
       width = 5,
       height = 5)





p4 <- data.frame(
  P = modelFits$BMkappaP[1:i][I] %>%
    laply(., function(x)
      x$modelOptim$kappa[1]),
  G = modelFits$BMkappaG[1:i][I] %>%
    laply(., function(x)
      x$modelOptim$kappa[1])
) %>%
  melt %>%
  ggplot(., aes(variable, value)) +
  # geom_boxplot()+
  geom_quasirandom(alpha = 1, size = 0.1) +
  geom_violin(scale = 3, alpha = 0.5) +
  scale_x_discrete(labels = c(expression(Sigma %prop% P), expression(Sigma %prop%
                                                                       G))) +
  ylab(expression(kappa)) +
  xlab("model")
p4
ggsave("results/kappas.pdf", p4, width = 5, height = 5)


p5 <- 
  rbind(
    ldply(modelFits$BM[1:i][I], function(x) {
      mod <- RetrieveBestModel(x)
      C <- mod$Sigma_x[, , 1] %>% tcrossprod
      m1 <- emorph2:::varSum(C[c(1, 4), c(1, 4)])
      m2 <- emorph2:::varSum(C[c(2, 5), c(2, 5)])
      m3 <- emorph2:::varSum(C[c(3, 6), c(3, 6)])
      data.frame(m1, m2, m3, matrix = "Rate matrix")}),
    ldply(modelFits$BMkappaP[1:i][I], function(x) {
      mod <- RetrieveBestModel(x) #%>% PCMApplyTransformation
      C <- mod$Sigma_x[, , 1] %>% tcrossprod
      m1 <- emorph2:::varSum(C[c(1, 4), c(1, 4)])
      m2 <- emorph2:::varSum(C[c(2, 5), c(2, 5)])
      m3 <- emorph2:::varSum(C[c(3, 6), c(3, 6)])
      data.frame(m1, m2, m3, matrix = "P matrix")
    }),
    ldply(modelFits$BMkappaG[1:i][I], function(x) {
      mod <- RetrieveBestModel(x) #%>% PCMApplyTransformation
      C <- mod$Sigma_x[, , 1] %>% tcrossprod
      m1 <- emorph2:::varSum(C[c(1, 4), c(1, 4)])
      m2 <- emorph2:::varSum(C[c(2, 5), c(2, 5)])
      m3 <- emorph2:::varSum(C[c(3, 6), c(3, 6)])
      data.frame(m1, m2, m3, matrix = "G matrix")
    })) %>%
  subset(., m1<0.02) %>%
  melt %>%
  ggplot(., aes(variable, value)) +
  geom_quasirandom(
    aes(color = matrix),
    alpha = 1,
    size = 0.1,
    dodge.width = 0.9) +
  geom_violin(aes(fill = matrix), scale = 3, alpha = 0, show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_color_manual(name="Matrix",values = c("#87CEEB", "#EEE685", "#436EEE")) +
  theme(legend.position = c(0.2,0.8),
        legend.box.background = element_rect(color="gray", size=2))+
  xlab("")+
  ylab("Rate of Evolution/Variance")
  # ylab(expression(sigma ^ 2)) +
  # guides(colour = guide_legend(override.aes = list(size = 1)))+
  # scale_y_continuous(
  #   # Features of the first axis
  #   name = "Rate of evolution",
  #   sec.axis = sec_axis( trans=~.*1, name="Variance")
  # )
p5
ggsave("results/sigma_areas.pdf", p5, width = 5, height = 5)

p5_red <-
  ldply(modelFits$BM[I], function(x) {
    mod <- RetrieveBestModel(x)
    C <- mod$Sigma_x[, , 1] %>% tcrossprod
    m1 <- emorph2:::varSum(C[c(1, 4), c(1, 4)])
    m2 <- emorph2:::varSum(C[c(2, 5), c(2, 5)])
    m3 <- emorph2:::varSum(C[c(3, 6), c(3, 6)])
    data.frame(m1, m2, m3, matrix = "BM")}) %>%
  subset(., m1<0.02) %>%
  melt %>%
  ggplot(., aes(variable, value)) +
  geom_quasirandom(
    aes(),
    alpha = 0.5,
    size = 0.2,
    dodge.width = 0.9) +
  geom_violin(scale = 3, alpha = 0) +
  ylab(expression(sigma ^ 2))

ggsave("results/sigma_areas_red.pdf", p5_red, width = 5, height = 5)

#######################

load("../Proportionality_diagonal/fitDist.Rdata")

i <- sum(!laply(modelFits$BM, is.null))
likBMmodel <- laply(modelFits$BM[1:i], function(X)
  X$logLikOptim[1])
likBMkappaPmodel <-
  laply(modelFits$BMkappaP[1:i], function(X)
    X$logLikOptim[1])
likBMkappaGmodel <-
  laply(modelFits$BMkappaG[1:i], function(X)
    X$logLikOptim[1])
n <- dim(X)[2]

modelCompare <- rbind(
  data.frame(
    model = "BM",
    logLik = likBMmodel,
    k = PCMParamCount(modelFits$BM[[1]]$modelOptim)
  ),
  data.frame(
    model = "BMkappaP",
    logLik = likBMkappaPmodel,
    k = PCMParamCount(modelFits$BMkappaP[[1]]$modelOptim)
  ),
  data.frame(
    model = "BMkappaG",
    logLik = likBMkappaGmodel,
    k = PCMParamCount(modelFits$BMkappaG[[1]]$modelOptim)
  )
) %>%
  mutate(
    .,
    AICc = (2 * k - 2 * logLik) + (2 * k) * (k + 1) / (n - k - 1),
    BIC = (log(n) * k - 2 * logLik)
  )

# I <- likBMmodel>0 & likBMkappaGmodel>0 & likBMkappaPmodel>0
# I <- likBMmodel>500 & likBMkappaGmodel>500 & likBMkappaPmodel>500
I <-
  likBMkappaGmodel > 190 &
  likBMmodel > likBMkappaPmodel &
  likBMmodel > likBMkappaGmodel

d2p2 <- modelCompare %>% select(., -k) %>%
  subset(., c(I, I, I)) %>%
  melt

CPCOLS <- c("#87CEEB", "#EEE685", "#436EEE")

p23<-rbind(data.frame(type="diagonal model",d2p2),
           data.frame(type="covariance model",d2p)) %>%
  subset(., variable=="BIC") %>%
  ggplot(., aes(model, value)) +
  facet_wrap(type~., scales = "free_y") +
  geom_violin(scale = 3,
              alpha = 1,
              show.legend = FALSE) +
  geom_quasirandom(
    aes(color = model),
    alpha = 1,
    size = 0.1,
    show.legend = F
  ) +
  xlab("") +
  ylab("BIC") +
  scale_color_manual(values = CPCOLS) +
  # coord_flip()+
  scale_x_discrete(labels = c("BM", expression(Sigma %prop% P), expression(Sigma %prop% G)))

ggsave("results/model_compare2.pdf",
       p23,
       width = 8,
       height = 5)


