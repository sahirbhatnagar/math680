##################################
# R source code file for section 1 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

rm(list=ls())

# import transformed data -------------------------------------------------

# setwd("~/Biostats PhD/MATH680/Assignment 2/")
setwd("~/Dropbox/Winter 2015/MATH680/Assignments/A2/")
source("trans_choles_data.R")

# Histogram of Compliance values to show its N(0,1) -----------------------

hist(DT$x, xlab="transformed compliance", main="Histogram of transformed compliance")

# Trying to guess 11 subject used for Figure 1 and 5 ----------------------------

subjects.comp <- c(-2.32316, -1.37117, -0.82660, -0.54559,-0.22347,-0.00764,0.24705, 
                  0.53688 ,0.74428, 1.28329, 2.06079) 

DT[x %in% subjects.comp]

# Reproduce Figure 1 ------------------------------------------------------

fit <- lm(y ~ x + x2 + x3, data=DT)
new.comp <- with(DT, seq(min(x),max(x), length.out=1000))
new.val <- predict(object = fit, newdata = data.frame(x=new.comp, x2=new.comp^2, x3=new.comp^3))

plot(y ~ x, data=DT, xlab="compliance", ylab="cholesterol decrease", ylim=c(-55,max(DT$y)))
lines(x=new.comp, y=new.val, col="green", lwd=4)
abline(h=0, lty=2)
calibrate::textxy(subjects.comp, rep(-45,length(subjects.comp)), 
                  labs = seq_along(1:length(subjects.comp)), cex = 0.8)
for (k in seq_along(1:length(subjects.comp))) abline(v=subjects.comp[k], lty=2)
points(-2.32316,11.50, pch=19, col="red")
arrows(-2,40,-2.25,15, angle=60)



