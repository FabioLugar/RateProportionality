library(ape)
library(geiger)
library(phytools)
library(nlme)
library(evomap)
library(caper)
library(geomorph)
library(geiger)
library(diversitree)
library(treeplyr)

## 
# data <- list()
# data$ratios <- read.table("./Data_log_ratios.txt",header=T,sep="\t",row.names=1)
# data$md <- read.table("./Data_log_md.txt",header=T,sep="\t",row.names=1)
# data$bl <- read.table("./Data_log_bl-tal.txt",header=T,sep="\t",row.names=1)
# bl-tal -> dw; bl-tri -> mw; md -> l

Gdata <- read.csv("./hluskoG.csv")
geoCV <- function(cv){
  log(cv^2+1)
}
geoCVse <- function(n, cv, se){
  .cv <- rnorm(n, cv, se)
  log(.cv^2+1)
}
Gdata <- mutate(Gdata, Var.residual = (1-total.c2)*Var, G.raw = total.h2*Var, CV=sqrt(Var)/Mean, CV.residual=sqrt(Var.residual)/Mean, 
                P.raw=CV^2, P.raw.residual=CV.residual^2, P.ln = geoCV(CV), P.ln.residual = geoCV(CV.residual), G.ln = geoCV(sqrt(G.raw)/Mean))

Gpooled <- group_by(Gdata, Molar, Trait) %>% summarize(., Var.pooled = sum((n-1)*Var)/(sum(n) - length(n)),  
                                                       Var.pooled.residual = sum((n-1)*Var.residual)/(sum(n) - length(n)),
                                                       G.raw = sum((n-1)*Var*total.h2)/(sum(n) - length(n)), 
                                                       CV.ln=sqrt(sum((n-1)*(Var/Mean^2)/(sum(n)-length(n)))),
                                                       CV.residual.ln=sqrt(sum((n-1)*(Var.residual/Mean^2)/(sum(n)-length(n)))),
                                                       G.ln = geoCV(sqrt(sum((n-1)*(Var*total.h2/Mean^2)/(sum(n) - length(n))))),
) %>% mutate(., P.raw=CV.ln^2, P.raw.residual=CV.residual.ln^2, P.ln = geoCV(CV.ln), P.ln.residual = geoCV(CV.residual.ln)) %>% arrange(., Trait)

h2samples <- lapply(1:nrow(Gdata), function(x) rnorm(100, Gdata$resid.h2[x], Gdata$resid.h2SE[x])*(1-Gdata$total.c2[x]))
G.lnsamples <- lapply(1:nrow(Gdata), function(x) Gdata$Var[x]*h2samples[[x]]/Gdata$Mean[x]^2)
G.lnSamp <- list()
for(i in 1:6){
  G.lnSamp[[i]] <- c(G.lnsamples[[i]], G.lnsamples[[6+i]])
}


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

Gcor
