##################################
# R source code file for section 1 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

#rm(list=ls())

# setwd("~/Biostats PhD/MATH680/Assignment 2/")
# setwd("~/git_repositories/math680/")

## ---- import-data ----
source("trans_choles_data.R")

## ---- hist-compliance ----
hist(DT$x, xlab="transformed compliance", main="Histogram of transformed compliance")

# Trying to guess 11 subject used for Figure 1 and 5 ----------------------------
#DT[x %in% subjects.comp]

## ---- figure-1 ----
subjects.comp <- c(-2.25093, -1.36986, -0.82647, -0.54556,-0.22347,-0.00764,0.24704, 
                   0.53672, 0.74327, 1.27808, 1.97051)

fit <- lm(y ~ x + x2 + x3, data=DT)
new.comp <- with(DT, seq(min(x),max(x), length.out=1000))
new.val <- predict(object = fit, newdata = data.frame(x=new.comp, x2=new.comp^2, x3=new.comp^3))

plot(y ~ x, data=DT, xlab="compliance", ylab="cholesterol decrease", ylim=c(-55,max(DT$y)))
lines(x=new.comp, y=new.val, col="green", lwd=4)
abline(h=0, lty=2)
calibrate::textxy(subjects.comp, rep(-45,length(subjects.comp)), 
                  labs = seq_along(1:length(subjects.comp)), cex = 0.8)
for (k in seq_along(1:length(subjects.comp))) abline(v=subjects.comp[k], lty=2)
points(-2.2509,11.50, pch=19, col="red")
arrows(-2.1,35,-2.25,15, angle=15)
calibrate::textxy(-2.0,36,labs = 1, cex=1)