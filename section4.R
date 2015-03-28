##################################
# R source code file for section 4 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 27
# Updated: 
# NOTE: 
##################################

# rm(list=ls())

# setwd("~/Biostats PhD/MATH680/Assignment 2/")
# setwd("~/git_repositories/math680/")


# ---- importfunctions ----
# only import functions and data to test this Rscript..
# comment out when compiling Rnw because these have been previously called
# library(MASS)
#source("functions.R")
#source("trans_choles_data.R")

## ---- figure-5 ----
subjects.comp <- c(-2.25093, -1.36986, -0.82647, -0.54556,-0.22347,-0.00764,0.24704, 
                   0.53672, 0.74327, 1.27808, 1.97051) 
id.unord <- match(subjects.comp, round(DT[,x], 5))
id.ord <- match(subjects.comp, round(sort(DT[,x]), 5))

#X <- cbind(1, as.matrix(DT[,-c(1,ncol(DT)), with=FALSE]))

X <- model.matrix(~x+x2+x3+x4+x5+x6,data=DT )

#Sigma depends on compliance

sigma.fn <- function(c){
    23.7 + 5.49*c - 2.25*c^2 - 1.03*c^3
}

Sigma.mat <- diag(sigma.fn(DT[,x]))

mu.bold <- X %*% coef(lm(y~., DT[,-ncol(DT), with=FALSE]))

G <- t(X) %*% solve(Sigma.mat) %*% X
beta.hat <- solve(G) %*% t(X) %*% solve(Sigma.mat) %*% as.matrix(DT[,y])
V <- solve(G)

#Simulation

Ystar.i <- mvrnorm(n = 100, mu = mu.bold, Sigma = Sigma.mat)
muHatStar.i <- apply(Ystar.i, 1, function(row) X %*% coef(lm.fit(y=row, x=X)))

Ystarstar.ij <- foreach(i = 1:100, .packages="MASS") %dopar% {
    mvrnorm(n = 1000, mu = muHatStar.i[,i], Sigma = Sigma.mat)
}

par.bs <- foreach(i = seq_along(1:length(Ystarstar.ij))) %dopar% {param(Ystarstar.ij[[i]])}

mu.smooth <- sapply(par.bs, function(l) l[,1])
sd.smooth <- sapply(par.bs, function(l) l[,2])

#Figure 5 with selected subjects
plot(x=1, y=1, type='n', ylim=c(0,2), xlim=c(min(DT[,x]), max(DT[,x])),
     xlab="adjusted compliance, subjects 1 through 11", ylab="standard deviation estimates")
for(i in 1:100) points(x=DT[,x][id.unord], y=sd.smooth[id.unord,i], 
                       pch='-', cex=0.7)
lines(x=DT[,x][order(DT[,x])][id.ord], 
      y=apply(sd.smooth, 1, mean)[order(DT[,x])][id.ord], lty=2, 
      col='blue', type='b', pch=19)
lines(x=DT[,x][order(DT[,x])][id.ord], 
      y=apply(mu.smooth, 1, sd)[order(DT[,x])][id.ord],
      col='red', lwd=2)
calibrate::textxy(DT[,x][id.unord], rep(0.05,length(DT[,x][id.unord])), 
                  labs = seq_along(1:length(DT[,x][id.unord])), cex = 0.8)