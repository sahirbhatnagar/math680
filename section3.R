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

JJ <- fit.all(data=DT, pred=DT[,x], bootsamples = samples, B=4000, single=FALSE)
se.bar <- summary(lm(y ~ x + x2 + x3, data=DT))$sigma

plot(x=DT[,x], y=JJ$sdtilde/JJ$sdhat, ylim=c(0,1), 
     col='blue', pch=19, cex=0.6)

ggplot(as.data.frame(x=DT[,x], y=JJ$sdtilde/JJ$sdhat), aes(x,y)) + geom_point(colour="red", fill="white") +
    geom_point(data = as.data.frame(x=DT[,x], y=JJ$sdhat/se.bar), aes(x,y) colour="blue", fill="white") + 
    geom_vline(xintercept = 0, colour="black", linetype = "longdash", size=1.5) + 
    geom_hline(yintercept = 1, colour = "black") + geom_hline(yintercept = 0, colour="black", linetype = "longdash") +
    xlab("adjusted compliance") + ylab("stdev ratio")