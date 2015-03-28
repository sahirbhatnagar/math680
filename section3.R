##################################
# R source code file for section 2 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

# rm(list=ls())

# setwd("~/Biostats PhD/MATH680/Assignment 2/")
# setwd("~/git_repositories/math680/")


# ---- importfunctions ----
# only import functions and data to test this Rscript..
# comment out when compiling Rnw because these have been previously called
#source("functions.R")
#source("trans_choles_data.R")
# B <- 4000
# samples <- replicate(B,DT[sample(1:nrow(DT),replace=T),],simplify=F)

# Compliance for observation 1
#obs <- -2.32316


## ---- table-3 ----
source("functions.R")
jj <- fit.all(data=DT,pred=obs,bootsamples = samples, B=4000)

#plot conf.intervals
p1 <- data.frame(type="Standard", Interval=paste0("(",round(jj$l.stand[1],2), ", ", round(jj$u.stand[1],2),")"), 
                 Length=jj$length[1], "Center point"=jj$center[1])
p2 <- data.frame(type="Percentile", Interval=paste0("(",round(jj$l.quant[1],2), ", ", round(jj$u.quant[1],2),")"), 
                 Length=jj$length.1[1], "Center point"=jj$center.1[1])
p3 <- data.frame(type="Smoothed", Interval=paste0("(",round(jj$l.smooth[1],2), ", ", round(jj$u.smooth[1],2),")"), 
                 Length=jj$length.2[1], "Center point"=jj$center.2[1])

p.conf <- rbind(p1,p2,p3)

 
## ---- figure-2 ----

# JJ <- fit.all(data=DT, pred=DT[,x], bootsamples = samples, B=4000, single=FALSE)
# se.bar <- summary(lm(y ~ x + x2 + x3, data=DT))$sigma

JJ <- read.table("figure2data", 
                 quote="\"", stringsAsFactors=FALSE)
colnames(JJ) <- c("compliance", "muhat", "mutilde", "sdhat", "sdtilde", "l.stand", "u.stand", "center", "length", "coverage", 
                  "sdbar", "l.quant", "u.quant", "center.1", "length.1", "coverage.1", "l.smooth", "u.smooth", 
                  "center.2", "length.2", "coverage.2")

X <- cbind(1, as.matrix(DT[, list(x, x2, x3)]))
c.vec <- function(c) c(1, c, c^2, c^3)

SEbar <- foreach(i = 1:nrow(JJ), .combine=c) %do% {
    sqrt(JJ$sdbar[i]^2 * t(c.vec(JJ$compliance[i])) %*% 
                  solve(t(X) %*% X) %*% 
                  c.vec(JJ$compliance[i]))
}

plot(x=JJ$compliance, y=JJ$sdtilde/JJ$sdhat, ylim=c(0,1.6), 
     col='red', pch=19, cex=0.6, xlab="adjusted compliance",
     ylab="stdev ratio")
lines(lowess(data.frame(x=JJ$compliance, y=JJ$sdtilde/JJ$sdhat)),
      col='red')
points(x=JJ$compliance, y=JJ$sdtilde/SEbar, 
     col='blue', pch=19, cex=0.6)
abline(h=c(1,0), lty=c(1,2), lwd=c(2,1))
abline(v=0, lty=2)
