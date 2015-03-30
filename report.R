##################################
# R source code file for user defined functions used in A2
# 
# Created by Sahir, March 22
# Updated: 
# NOTE: 
##################################

library(plyr)
library(data.table)
library(doParallel)
library(foreach)
library(bootstrap)
library(splines)
library(ggplot2)
library(reshape2)
library(car)
library(grid)
library(calibrate)
library(xtable)
library(MASS)
#library(glmnet)
registerDoParallel(cores = 4)



# K fold cross validation -------------------------------------------------
my.kfold <- function(DF, K=5, REP=5, formulas) {
    X <- DF[,c("x","x2","x3","x4","x5","x6")]
    y <- DF[,"y"]
    n <- nrow(DF)
    CVErr <- 0
    for (iREP in 1:REP) {
        fold <- sample(rep(1:K, length = n))
        SumSqErr <- 0 
        for (k in 1:K) {
            iTest <- fold == k
            Xk <- X[!iTest, , drop = FALSE]
            yk <- y[!iTest]
            
            #data frame with k observations removed
            Xyk <- data.frame(as.data.frame(Xk), y = yk)
            
            #fit all possible models for the Xyk dataset
            fits <- lapply(formulas, function(i) lm(as.formula(i),data=Xyk))    
            
            # calculate the predicted responses for each of the models
            yHat <- lapply(fits, function(i) predict(i, newdata = X[iTest, ,drop = FALSE], type = "response"))
            
            #sum of the residuals for the k subset for all 32 models
            SumSqErr <- SumSqErr + sapply(yHat, function(i) sum((y[iTest]-i)^2))
        }
        #CVk for one repitition, for 32 models 
        CVErr <- CVErr + SumSqErr/n
    }
    CVErr/REP
}

# this function works for the compliance dataset ----
fit.best <- function(j, method="AIC", K, REP, predict=-2.25093){
    
    obs <- predict
    DF <- j
    n <- nrow(DF)
    y <- DF[,.(y)]
    x <- DF[,.(x,x2,x3,x4,x5,x6)]
    
    # 6 candidate models
    formulas <- list("y~x", "y~x+x2", "y~x+x2+x3", "y~x+x2+x3+x4", 
                     "y~x+x2+x3+x4+x5", "y~x+x2+x3+x4+x5+x6")
    
    # fit all candidate models
    fits <- lapply(formulas, function(i) lm(as.formula(i),data=DF))
    
    # log likelihoods for the models, produces same result as logLik function 
    # in R base (sapply(fits,logLik))
    log.liks <- sapply(fits, function(i) 
        { (-n/2)*log(sum(i$residuals^2)/n) - n/2 - (n/2)*log(2*pi)} )
    
    # coefficient values
    dfCoefNum <- ldply(fits, function(x) as.data.frame( t(coef(x))))
    
    # number of coefficients per model including intercept
    p <- apply(dfCoefNum, 1, function(i) sum(!is.na(i)))
    
    #selection criterion
    
    if (method=="AIC") {criterion <- -2*log.liks + 2*(p+1)}
    
    if (method=="BIC") {criterion <- -2*log.liks + log(n)*(p+1)}
    
    if (method=="CV")  {criterion <- sapply(fits, function(i) 
        {(1/n)*sum(  ((y-fitted.values(i))^2)/(1-influence(i)$hat)^2)})}
    
    if (method=="CP")  {
        # sigma squared (which here I call sigma_full) full model
        # sigma_full <- sum(residuals(fits[[6]])^2)/(164-6-1)
        # sigma of 22 was in all boostrap replications pg 993 above table 1
        sigma_full <- 22^2
        ssres <- sapply(fits, function(i) { sum(i$residuals^2) }) 
        criterion <- ssres + (2*p)*sigma_full
    }
    
    if (method=="GCV") {
        ssres <- sapply(fits, function(i) { (1/n)*sum(i$residuals^2) }) 
        criterion <- ssres/(1-p/n)^2
    }
    
    
    if (method=="KFOLD"){
        criterion <- my.kfold(DF=DF, K=K, REP=REP, formulas=formulas)
    }
    
    # fitted value for observation 
    muhat <- sapply(1:6,function(i) predict(fits[[i]], 
                    newdata=data.frame(x=obs, x2=obs^2, x3=obs^3,x4=obs^4, 
                                       x5=obs^5,x6=obs^6)))
    
    # create data frame of results
    dat <- !is.na(dfCoefNum)
    dat <- data.frame(dat,criteria=criterion, muhat=muhat)
    
    # which model gets selected
    v.min <- which(dat[,"criteria"]==min(dat[,"criteria"]))
    
    # label the model that gets selected
    # y~ LINEAR
    m1 <- nrow(subset(dat[v.min,], x==T & x2==F & x3==F & x4==F & x5==F & x6==F))>0
    # y~ QUADRATIC
    m2 <- nrow(subset(dat[v.min,], x==T & x2==T & x3==F & x4==F & x5==F & x6==F))>0
    # y~ CUBIC
    m3 <- nrow(subset(dat[v.min,], x==T & x2==T & x3==T & x4==F & x5==F & x6==F))>0
    # y~ QUARTIC
    m4 <- nrow(subset(dat[v.min,], x==T & x2==T & x3==T & x4==T & x5==F & x6==F))>0
    # y~ QUINIC
    m5 <- nrow(subset(dat[v.min,], x==T & x2==T & x3==T & x4==T & x5==T & x6==F))>0
    # y~ SEXTIC
    m6 <- nrow(subset(dat[v.min,], x==T & x2==T & x3==T & x4==T & x5==T & x6==T))>0
    
    # muhat for chosen model
    muhat <- dat[v.min, "muhat"]
    
    l <- list("m1" = m1,"m2" = m2,"m3" = m3,"m4" = m4,"m5" = m5,"m6" = m6, 
              "muhat"=muhat)
    return(l)
}

# this function is used to get the numbers for Table 1, it outputs a data.frame ----
# with the Cp criteria for each model, as well as a predicted value based on the fitted model
fit.once <- function(j, predict=-2.25093, single=TRUE){
    
    obs <- predict
    DF <- j
    n <- nrow(DF)
    y <- DF[,.(y)]
    x <- DF[,.(x,x2,x3,x4,x5,x6)]

    # 6 candidate models
    formulas <- list("y~x", "y~x+x2", "y~x+x2+x3", "y~x+x2+x3+x4", 
                     "y~x+x2+x3+x4+x5", "y~x+x2+x3+x4+x5+x6")
    
    # fit all candidate models
    fits <- lapply(formulas, function(i) lm(as.formula(i),data=DF))
    
    # log likelihoods for the models, produces same result as logLik function 
    # in R base (sapply(fits,logLik))
    log.liks <- sapply(fits, function(i) 
    { (-n/2)*log(sum(i$residuals^2)/n) - n/2 - (n/2)*log(2*pi)} )
    
    # coefficient values
    dfCoefNum <- ldply(fits, function(x) as.data.frame( t(coef(x))))
    
    # number of coefficients per model including intercept
    p <- apply(dfCoefNum, 1, function(i) sum(!is.na(i)))
    
    #selection criterion
    sigma_full <- sapply(fits, function(i) car::sigmaHat(i)^2)
    ssres <- sapply(fits, function(i) { sum(i$residuals^2) }) 
    criterion <- ssres + (2*p)*sigma_full
    
    #SE under given model
    se.model <- sapply(fits, function(i) summary(i)$sigma)
    
    # fitted value for observation 
    if(single){
        muhat <- sapply(1:6,function(i) predict(fits[[i]], 
                                                newdata=data.frame(x=obs, x2=obs^2, x3=obs^3,x4=obs^4, 
                                                                   x5=obs^5,x6=obs^6)))
        # create data frame of results
        dat <- !is.na(dfCoefNum)
        dat <- data.frame(dat,criteria=criterion, muhat=muhat, se.bar=se.model, p=p)
        
    }else{
        dat <- foreach(l = obs) %dopar% {
            muhat <- sapply(1:6,function(i){
                predict(fits[[i]],
                        newdata=data.frame(x=l, x2=l^2, x3=l^3,x4=l^4,
                                           x5=l^5,x6=l^6))
            })
        data.frame(!is.na(dfCoefNum),criteria=criterion, muhat=muhat, se.bar=se.model, p=p)
        }
    }
    
    return(dat)
}

# Evaluate center, length, coverage ----------------------------------------
vals <- function(x, tstar, n.boot){
    #length
    len <- abs(x[2]-x[1])
    #center
    cen <- (x[2]+x[1])/2
    
    coverage.stand <- sapply(tstar, function(i) { i>= x[1] & i<=x[2]  })
    #percent coverage
    cov <- sum(coverage.stand)/n.boot
    
    return(data.frame(center=cen, length=len, coverage=cov))
}

# calculates sdhat, sdsmooth, mutild, musmooth, quantiles CI for Lasso, MCP, SCAD, bic, cp, --------------
# This function is for NON-PARAMETRIC BOOTSTRAP and works for prostate dataset
fit.all <- function(data, pred , bootsamples, single = TRUE, B, penalty) {
    #single: if TRUE then predicts single observation, else predicts fitted values for all j=1,...,n
    #pred: data frame of observation that you want to predict
    #n.lambda: number of lambdas to test for tuning parameter in lasso
    #B: number of bootstrap samples
    # bootsamples: bootstrap sample
    #penalty: SCAD, lasso, MCP, bic, cp
    
    #data=DT; pred=-2.32316;single=TRUE;B=4000;bootsamples=samples
    ############################################################
    x <- data[,c(-1,-8), with=FALSE]
    y <- data[,y]
    n <- nrow(data)
    #cov.names <- c("x","x2","x3","x4","x5","x6")
    
    if (single){
        
        # model selection on original data
        # \hat{\mu} = t(y) , predict one observation
        #t.y <- pred.penalty(data=data, x=x, y=y, pred=pred, penalty=penalty)
        mat.results <- fit.once(data, predict = pred)
        t.y <- mat.results[which.min(mat.results$criteria),]$muhat
        #se.bar <- mat.results[which.min(mat.results$criteria),]$se.bar
        
        # t_i^*(y_i^*), i.e., the predicted values for 1 individual j for each 
        # bootstrap sample i=1,...,B, 
        # need to change this code for different datasets
        tstar.i <- foreach(i = bootsamples, .export = "fit.once",
                           .packages=c("data.table", "plyr"), .combine=c) %dopar% {
            mat.results <- fit.once(i, predict = pred)
            mat.results[which.min(mat.results$criteria),]$muhat
        }
        #tstar.i <- unlist(tstar.i)
        
        # s(y) = B^{-1} \sumn t_i^*(y_i^*) for one j
        # this is the mean for one predicted value across each bootstrap sample, i.e. , the mutilde
        mutilde <- mean(tstar.i)
        # \hat{sd_B}, the standard bootstrap standard error
        sd.hat <- sd(tstar.i)
        
        # Standard interval ----
        stand.int <- t.y+qnorm(c(0.025, 0.975))*sd.hat
        # center, length and coverage
        m <- vals(stand.int, tstar.i, n.boot=B)
        
        # quantile interval
        quantile.int <- quantile(tstar.i, probs=c(0.025,0.975))
        quantile.result <- as.data.frame(vals(quantile.int, tstar.i, n.boot=B))
        quantile.int <- data.frame(quantile.int)
        
        #smoothed interval    
        label <- seq(1,n+1, by=1) 
        #how many times j appears in bootstrap sample i, j=1,..., 164
        y.ij <- foreach(i = bootsamples) %dopar% as.data.frame(table(cut(i$label,label,include.highest = TRUE, right=FALSE)))[,2]
        #data.frame version of above, these are the Y_{ij}^{*}
        y.ij.star <- sapply(y.ij,cbind)
        
        #how many times each observation j=1,...,164, showed up in bootstrap sample i=1,...,4000
        # y.j.star is Y__{\cdot j}^{*}, i.e.,  (1/n.boot) \sum_{i=1}^{n.boot} y_{ij}^*
        y.j.star <- apply(as.matrix(y.ij.star),1,sum)/B
        
        covj <- rep(NA,n)
        
        for (j in 1:n){
            temp = 0
            for (i in 1:B){
                temp <- temp + (y.ij.star[j,i]-y.j.star[j])*(tstar.i[i]-mutilde)
            }
            covj[j] <- temp/B
        }
        
        #\tilde{sd_B} smoothed standard deviation
        sd.tilde <- sqrt(sum(covj^2))
        
        #smoothed interval
        smooth.int <- mutilde+qnorm(c(0.025,0.975))*sd.tilde 
        smooth.result <- vals(smooth.int, tstar.i, n.boot=B)
        
        
        result <- data.frame(muhat=t.y, mutilde=mutilde, sdhat=sd.hat, sdtilde=sd.tilde,  
                             l.stand=stand.int[1], u.stand=stand.int[2], m, sdbar=mat.results$se.bar[3],
                             l.quant=quantile.int[1,], u.quant=quantile.int[2,],quantile.result,
                             l.smooth=smooth.int[1], u.smooth=smooth.int[2], smooth.result)
        
        
    } else { 
               
        # model selection on original data
        # \hat{\mu} = t(y) , predict one observation
        #t.y <- pred.penalty(data=data, x=x, y=y, pred=pred, penalty=penalty)
        mat.results <- fit.once(data, predict = pred, single=FALSE)
        t.y <- sapply(mat.results, function(i) i[which.min(i$criteria),]$muhat)
        #se.bar <- mat.results[which.min(mat.results$criteria),]$se.bar
        
        # t_i^*(y_i^*), i.e., the predicted values for 1 individual j for each 
        # bootstrap sample i=1,...,B, 
        # need to change this code for different datasets
        tstar.i <- foreach(i = bootsamples, .export = "fit.once",
                           .packages=c("data.table", "plyr", "foreach")) %dopar% {
                               mat.results <- fit.once(i, predict = pred, single=FALSE)
                               sapply(mat.results, function(l) l[which.min(l$criteria),]$muhat)
                           }
        tstar.i <- matrix(unlist(tstar.i), byrow=TRUE, nrow=length(bootsamples))
        
        # s(y) = B^{-1} \sumn t_i^*(y_i^*) for one j
        # this is the mean for one predicted value across each bootstrap sample, i.e. , the mutilde
        mutilde <- apply(tstar.i, 2, mean)
        # \hat{sd_B}, the standard bootstrap standard error
        sd.hat <- apply(tstar.i, 2, sd)
        
        # Standard interval ----
        stand.int <- cbind(t.y+qnorm(0.025)*sd.hat, t.y+qnorm(0.975)*sd.hat)
        # center, length and coverage
        m <- foreach(i = 1:nrow(stand.int), .combine = rbind,
                     .export="vals") %dopar% vals(stand.int[i,], tstar.i[,i], n.boot=B)
        
        # quantile interval
        quantile.int <- t(apply(tstar.i, 2, quantile, probs=c(0.025,0.975)))
        quantile.result <- foreach(i = 1:nrow(quantile.int), .combine = rbind,
                                   .export="vals") %dopar% vals(quantile.int[i,], tstar.i[,i], n.boot=B)
        #quantile.int <- data.frame(quantile.int)
        
        #smoothed interval    
        label <- seq(1,n+1, by=1) 
        #how many times j appears in bootstrap sample i, j=1,..., 164
        y.ij <- foreach(i = bootsamples) %dopar% as.data.frame(table(cut(i$label,label,include.highest = TRUE, right=FALSE)))[,2]
        #data.frame version of above, these are the Y_{ij}^{*}
        y.ij.star <- sapply(y.ij,cbind)
        
        #how many times each observation j=1,...,164, showed up in bootstrap sample i=1,...,4000
        # y.j.star is Y__{\cdot j}^{*}, i.e.,  (1/n.boot) \sum_{i=1}^{n.boot} y_{ij}^*
        y.j.star <- apply(as.matrix(y.ij.star),1,sum)/B
        
        sd.tilde <- foreach(l = 1:ncol(tstar.i), .combine = c) %dopar% {

        covj <- rep(NA,n)
        
        for (j in 1:n){
            temp = 0
            for (i in 1:B){
                temp <- temp + (y.ij.star[j,i]-y.j.star[j])*(tstar.i[i,l]-mutilde[l])
            }
            covj[j] <- temp/B
        }
        
        #\tilde{sd_B} smoothed standard deviation
        sqrt(sum(covj^2))
        }
        
        #smoothed interval
        smooth.int <- cbind(mutilde+qnorm(0.025)*sd.tilde, mutilde+qnorm(0.975)*sd.tilde) 
        smooth.result <- foreach(i = 1:nrow(smooth.int), .export="vals",
                                 .combine = rbind) %dopar% vals(smooth.int[i,], tstar.i[,i], n.boot=B)        
        
        result <- data.frame(muhat=t.y, mutilde=mutilde, sdhat=sd.hat, sdtilde=sd.tilde,  
                             l.stand=stand.int[1,], u.stand=stand.int[2,], m,
                             l.quant=quantile.int[1,], u.quant=quantile.int[2,],quantile.result,
                             l.smooth=smooth.int[1,], u.smooth=smooth.int[2,], smooth.result)
            }
    
    return(result)
    
}

#parametric bootstrap----

param <- function(data){
    
    #compute beta
    betaStar.i <- apply(data, 1, function(row){
        solve(G) %*% t(X) %*% solve(Sigma.mat) %*% row
    })
    betaStar.i <- t(betaStar.i)
    
    betaStar.mean <- apply(betaStar.i, 2, mean)
    B.mat <- apply(betaStar.i, 1, function(row) row - betaStar.mean)
    B.mat <- t(B.mat)
    
    #compute smoothed estimate
    muHat.i <- apply(betaStar.i, 1, function(row) X %*% row)
    muTilde.i <- apply(muHat.i, 1, mean)
    
    #empirical covariance
    cov.hat <- foreach(i = 1:164) %dopar% {
        t(B.mat) %*% (muHat.i[i,]-muTilde.i[i])/1000
    }
    
    #empirical bs variance
    V.bar <- t(B.mat) %*% B.mat/1000
    
    sd.tilde <- sapply(cov.hat, function(cov){
        sqrt(t(cov) %*% solve(V.bar) %*% cov)
    })
    
    return(cbind(muTilde.i, sd.tilde))
}

##################################
# R source code file for transforming Cholesterol data 
# as in Efron JASA 2014 paper
# 
# Created by Maxime, March 20
# Updated: 
# NOTE: 
##################################

library(bootstrap)

fit <- lm(y ~ z + I(z^2) + I(z^3), data=cholost)
new.comp <- seq(0,100, length.out=1000)
new.val <- predict(object = fit, newdata = data.frame(z=new.comp))

data <- cholost
data[95,"y"] <- data[95,"y"]+45
data$transf <- qnorm((rank(data$z) - 0.5)/nrow(data))


fit2 <- lm(y ~ transf + I(transf^2) + I(transf^3), data)
new.comp2 <- seq(min(data$transf),max(data$transf), length.out=1000)
new.val2 <- predict(object = fit2, newdata = data.frame(transf=new.comp2))

# get data into form so that it works with Sy's functions
DT <- as.data.table(data[,c("y","transf")])
setnames(DT,"transf","x")

# make the transformation now. This is how my function was written up
DT <- transform(DT, x2=x^2, x3=x^3, x4=x^4, x5=x^5, x6=x^6)

# round the compliances to 5 digits, so that we can find the 11 subjects for Figure 5
DT[,x:=round(x,5)]

# need to add a label to each data point. This is used in the calculation of 
# Efron's smoothed standard error estimate, i.e., we need to know how many
# times each subject appears in the bootstrap sample
DT[,label:=seq_along(1:nrow(DT))]


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

##################################
# R source code file for section 2 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

## ---- bootstrap-samples ----
B <- 4000
samples <- replicate(B,DT[sample(1:nrow(DT),replace=T),],simplify=F)


## ---- Cp-boot ----
source("functions.R")

#obs <- -2.32316
obs <- -2.25093

# Chosen models based on bootstrap replicates
true.cp <- foreach(i = samples, .packages = c("data.table", "plyr")) %dopar% fit.best(i, method="CP", predict=obs)
df.cp <- as.data.table(ldply(true.cp, data.frame))

# how many times each model was chosen 
perc.chosen <- apply(df.cp[,-7,with=FALSE], 2, sum)/B

# percentile interval for bootstrap replicates (also used for Figure 3)
ci.muhat <- data.frame(x=c(quantile(df.cp$muhat,c(0.025,0.975))[1],
               quantile(df.cp$muhat,c(0.025,0.975))[2]),y=c(0,0))


## ---- Cp-original ----
fit.all <- fit.once(DT)
muhat.obs <-fit.all[which.min(fit.all$criteria),"muhat"]


# table 1 ----
tab1 <- data.frame(model=c("Linear","Quadratic","Cubic","Quartic","Quintic",
                           "Sextic"),
                   m=fit.all$p, Cp=round(fit.all$criteria,0)-80000, "Bootstrap"=round(perc.chosen*100))

# table 2 ----
tab2 <- data.frame(rbind(sapply(paste0("m",1:6), function(i) df.cp[which(get(i)),round(mean(muhat),2)]),
sapply(paste0("m",1:6), function(i) df.cp[which(get(i)),round(sd(muhat),2)])))
rownames(tab2) <- c("Mean","St.dev.")

#some bootstrap samples that led to the sextic model
#being selected give a terrible fit: 25, 675, 757, 1567, 2191, 2414, 3716
#same thing for those that led to quintic model: 1538, 2778, 3096, 3315, 3417

#removing them and recomputing table 2 shows that the discrepancy
#between our table 2 and Efron's is due to bootstrap variability
tab21 <- data.frame(rbind(sapply(paste0("m",1:6), 
                                 function(i){
                                     df.cp[-c(25,675,757,1567,2191,2414,
                                              3716,1538,2778,3096,3315,3417)][which(get(i)),
                                                                              round(mean(muhat),2)]}),
                         sapply(paste0("m",1:6), 
                                function(i){
                                    df.cp[-c(25,675,757,1567,2191,2414,
                                             3716,1538,2778,3096,3315,3417)][which(get(i)),
                                                                             round(sd(muhat),2)]})))
rownames(tab21) <- c("Mean","St.dev.")

## ---- figure-3 ----
# histogram of 4000 muhat bootstrap estimates: Figure 3
ggplot(df.cp, aes(muhat)) + geom_histogram(binwidth=1, colour="black", fill="white")+xlim(-31,31)+
    geom_vline(xintercept = muhat.obs, colour="red", linetype = "longdash",size=1.5)+
    annotate("text", label = round(muhat.obs,2), x = 5, y = -7, size = 5, colour = "red")+
    geom_point(data=ci.muhat, aes(x,y), size=4,shape=24, fill="red")+
    xlab("bootstrap estimates for subject 1") + ylab("Frequency")

##################################
# R source code file for section 3 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

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

##################################
# R source code file for section 4 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 27
# Updated: 
# NOTE: 
##################################

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

##################################
# R source code file for discussion of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 27
# Updated: 
# NOTE: 
##################################

## ---- discussion-1 ----
m1 <- data.frame(muhat=df.cp[m1==TRUE,.(muhat)], model="linear")
m2 <- data.frame(muhat=df.cp[m2==TRUE,.(muhat)], model="quadratic")
m3 <- data.frame(muhat=df.cp[m3==TRUE,.(muhat)], model="cubic")
m4 <- data.frame(muhat=df.cp[m4==TRUE,.(muhat)], model="quartic")
m5 <- data.frame(muhat=df.cp[m5==TRUE,.(muhat)], model="quintic")
m6 <- data.frame(muhat=df.cp[m6==TRUE,.(muhat)], model="sextic")

dat2 <- rbind(m1,m3,m5)

#histograms of three most chosen models
ggplot(dat2, aes(muhat, fill=model)) + geom_histogram(alpha = .8, position = 'identity')+
    geom_vline(xintercept = muhat.obs, colour="red", linetype = "longdash")+
    geom_vline(xintercept = muhat.obs, colour="black", size=1)+
    labs(x=expression(paste("fitted value   ",hat(mu)[1])), y="frequency")+
    annotate("text", label = "original estimate = 2.71\n based on Cp chosen model", x = -13, y = 210, 
             size = 4, colour = "red")+
    geom_segment(aes(x = -2.51, y = 210, xend = 2.71, yend = 210),
                 color="red", arrow = arrow(angle=20,length = unit(0.5, "cm")))+
    annotate("text", label = paste0(as.numeric(round(perc.chosen["m1"]*100,0)),"%"), 
             x = -15, y = 75, size = 4, colour = "red")+
    annotate("text", label = paste0(as.numeric(round(perc.chosen["m3"]*100,0)),"%"), 
             x = -5, y = 75, size = 4, colour = "red")+
    annotate("text", label = paste0(as.numeric(round(perc.chosen["m5"]*100,0)),"%"), 
             x = 8, y = 75, size = 4, colour = "red") + xlim(c(-30,30))

prostate.data <- read.table("prostate.txt", header=TRUE)
cov.names <- names(prostate.data)[2:9]

################################################################################
# Prostate Cancer Analysis -----------------------------------------------------
################################################################################

# standardized version
prostate.s <- data.frame(label=seq(1,nrow(prostate.data),by=1), log.psa=log(prostate.data$psa),sapply(prostate.data[,cov.names], scale))

# standardized version, center log.psa
prostate.s.c <- data.frame(label=seq(1,nrow(prostate.data),by=1), log.psa=scale(log(prostate.data$psa), center=TRUE, scale=FALSE) ,sapply(prostate.data[,cov.names], scale))

# unstandardized version
prostate.u <- data.frame(label=seq(1,nrow(prostate.data),by=1), log.psa=log(prostate.data$psa),prostate.data[,2:9])

# unstandardized version, center log.psa
prostate.u.c <- data.frame(label=seq(1,nrow(prostate.data),by=1), log.psa=scale(log(prostate.data$psa), center=TRUE, scale=FALSE) ,prostate.data[,2:9])

#check that transformed have mean 0, variance 1
apply(prostate.u.c, 2, mean)

#best models of each size 
allsubsets <- regsubsets(as.formula(paste("log.psa ~", paste0(cov.names, collapse="+"), sep = "")), 
                         data = prostate.s, nbest=max(choose(8, 1:8)), really.big=T)
allsubsets.summary<-summary(allsubsets)
allsubsets.summary$which[which.min(allsubsets.summary$cp),]

B <- 4000
samples <- replicate(B,prostate.s[sample(1:nrow(prostate.s),replace=T),],simplify=F)
#predicted value based on the following observation
obs <- prostate.s[95,c(-1,-2)]


#CP
true.cp <- foreach(i = samples) %dopar% fit.best.pros(i, method="CP", predict=obs)
df.cp <- ldply (true.cp, data.frame)
apply(df.cp, 2, sum)/B


#chosen model based on CP is model 7
fit7 <- lm(log.psa~lcavol+lweight+svi+lbph+age+pgg45+lcp,data=prostate.s) 

plot(fit7)

#predict for model 7 (remove gleason)
muhat.obs<-predict(fit7, newdata=prostate.s[95,c(-1,-2,-9)]);muhat.obs

#histogram of 4000 muhat bootstrap estimates
ggplot(df.cp, aes(muhat)) + geom_histogram(binwidth=1/20, colour="red", fill="#00CC00")+
  geom_vline(xintercept = muhat.obs, colour="black", size=2)+
labs(title="Fitted values for subject 95, from B=4000 nonparametric bootsrap replications\n of the Cp chosen model; 60% of the replications greater than the \n original estimate 3.6", x=expression(paste("fitted value   ",hat(mu)[95])))+
  annotate("text", label = "original estimate = 3.6\n based on Cp chosen model", x = 2.9, y = 260, size = 6, colour = "red")+
  geom_segment(aes(x = 3.4, y = 260, xend = 3.6, yend = 260),color="red", arrow = arrow(angle=20,length = unit(0.5, "cm")))

#60% of greater less than 
ecdf(df.cp$muhat)(muhat.obs)

boxplot(df.cp$muhat)
abline(h=muhat.obs, col="red")

#muhat bootstrap based on 4000 samples
mean(df.cp$muhat)
#sd bootstrap
sd(df.cp$muhat)

m1 <-data.frame(muhat=subset(df.cp,m1==TRUE)[,"muhat"], model="m1")
m2 <-data.frame(muhat=subset(df.cp,m2==TRUE)[,"muhat"], model="m2")
m3 <-data.frame(muhat=subset(df.cp,m3==TRUE)[,"muhat"], model="m3")
m4 <-data.frame(muhat=subset(df.cp,m4==TRUE)[,"muhat"], model="m4")
m5 <-data.frame(muhat=subset(df.cp,m5==TRUE)[,"muhat"], model="m5")
m6 <-data.frame(muhat=subset(df.cp,m6==TRUE)[,"muhat"], model="m6")
m7 <-data.frame(muhat=subset(df.cp,m7==TRUE)[,"muhat"], model="m7")
m8 <-data.frame(muhat=subset(df.cp,m8==TRUE)[,"muhat"], model="m8")

dat2 <- rbind(m2,m3,m4,m5,m6,m7,m8)

p <- ggplot(dat2, aes(model, muhat))
p + geom_boxplot() + annotate("text", label = "1%", x = 1, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "18%", x = 2, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "12%", x = 3, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "22%", x = 4, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "15%", x = 5, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "24%", x = 6, y = 4.6, size = 6, colour = "red")+
  annotate("text", label = "8%", x = 7, y = 4.6, size = 6, colour = "red")+
  labs(title="Boxplot of fitted values for Subject 95 for the model chosen by Cp criteria
        based on B=4000 nonparametric bootsrap samples", y=expression(paste("fitted value   ",hat(mu)[95])))+
  annotate("text", label = "Model 7", x = 6, y = 3.4, size = 4, colour = "blue")+
  annotate("text", label = "******", x = 6, y = 3.2, size = 6, colour = "blue")


#density of three most chosen models
dat <- rbind(m3,m5,m7)
ggplot(dat, aes(muhat, fill=model)) + geom_density(alpha=0.3)

#histograms of three most chosen models
ggplot(dat, aes(muhat, fill=model)) + geom_histogram(binwidth=1/20,alpha = .8, position = 'identity')+
  geom_vline(xintercept = muhat.obs, colour="red", linetype = "longdash")+
  geom_vline(xintercept = muhat.obs, colour="black", size=1.5)+
  labs(title="Fitted values for subject 95, from B=4000 nonparametric bootsrap \n replications separated by three most frequently chosen models by Cp", x=expression(paste("fitted value   ",hat(mu)[95])))+
  annotate("text", label = "original estimate = 3.6\n based on Cp chosen model", x = 2.9, y = 100, size = 5, colour = "red")+
  geom_segment(aes(x = 3.4, y = 100, xend = 3.6, yend = 100),color="red", arrow = arrow(angle=20,length = unit(0.5, "cm")))+
  annotate("text", label = "18%", x = 4, y = 50, size = 6, colour = "red")+
  annotate("text", label = "22%", x = 3.75, y = 69, size = 6, colour = "red")+
  annotate("text", label = "24%", x = 3.2, y = 30, size = 6, colour = "red")


# Smooth intervals for prostate data LASSO, SCAD, MCP, bic, cp ----
obs <- prostate.u[95,c(-1,-2)]
#LASSO, need to use unstandardized
pros_lasso <- fit.all(data=prostate.u, x = prostate.u[,c(-1,-2)], y=prostate.u[,2] , 
                       pred=obs, single=TRUE, n.lambda=200, B=4000, penalty="lasso")

#MCP, need to use unstandardized data
pros_MCP <- fit.all(data=prostate.u, x = prostate.u[,c(-1,-2)], y=prostate.u[,2] , 
                    pred=obs, single=TRUE, n.lambda=200, B=4000, penalty="MCP")

#SCAD, need to use unstandardized data
pros_SCAD <- fit.all(data=prostate.u, x = prostate.u[,c(-1,-2)], y=prostate.u[,2] , 
                    pred=obs, single=TRUE, n.lambda=200, B=4000, penalty="SCAD")

#cp, need to use standardized data
pros_cp <- fit.all(data=prostate.s, x = prostate.s[,c(-1,-2)], y=prostate.s[,2] , 
                    pred=obs, single=TRUE, n.lambda=200, B=4000, penalty="cp")

#bic, need to use standardized data
pros_bic <- fit.all(data=prostate.s, x = prostate.s[,c(-1,-2)], y=prostate.s[,2] , 
                    pred=obs, single=TRUE, n.lambda=200, B=4000, penalty="bic")




#histogram of B muhat bootstrap estimates
ggplot(pros_bic, aes(tstar.i)) + geom_histogram(alpha=0.5)+
  geom_vline(xintercept = pros_bic$muhat, colour="red", linetype = "longdash")

ggplot(pros_cp, aes(tstar.i)) + geom_histogram(alpha=0.5)+
  geom_vline(xintercept = pros_cp$muhat, colour="red", linetype = "longdash")

ggplot(df.cp, aes(muhat)) + geom_histogram(binwidth=1/20, colour="red", fill="#00CC00")+
  geom_vline(xintercept = muhat.obs, colour="black", size=2)+
  labs(title="Fitted values for subject 95, from B=4000 nonparametric bootsrap replications\n of the Cp chosen model; 60% of the replications greater than the \n original estimate 3.6", x=expression(paste("fitted value   ",hat(mu)[95])))+
  annotate("text", label = "original estimate = 3.6\n based on Cp chosen model", x = 2.9, y = 260, size = 6, colour = "red")+
  geom_segment(aes(x = 3.4, y = 260, xend = 3.6, yend = 260),color="red", arrow = arrow(angle=20,length = unit(0.5, "cm")))

#all on one graph
dat.hist <- data.frame(MCP=pros_MCP$tstar.i, SCAD=pros_SCAD$tstar.i, LASSO=pros_lasso$tstar.i, BIC=pros_bic$tstar.i, Cp=pros_cp$tstar.i)
dat.hist.m <- melt(dat.hist)

dat.hist2 <- data.frame(MCP=pros_MCP$muhat[1], SCAD=pros_SCAD$muhat[1], LASSO=pros_lasso$X1[1], BIC=pros_bic$muhat[1], Cp=pros_cp$muhat[1])
dat.hist2.m <- melt(dat.hist2)

ggplot(subset(dat.hist.m, variable %in% c("SCAD","LASSO","MCP")), aes(value)) + geom_histogram(binwidth=1/30, colour="red", fill="#00CC00")+
  facet_wrap(~ variable)+  
  labs(title="Fitted values for subject 95, from B=4000 nonparametric bootsrap replications", x=expression(paste("fitted value   ",hat(mu)[95])))+
  geom_vline(data=subset(dat.hist2.m, variable %in% c("SCAD","LASSO","MCP")), aes(xintercept=value,label=variable), colour="black", size=1)


ggplot(subset(dat.hist.m, variable %in% c("BIC","Cp")), aes(value)) + geom_histogram(binwidth=1/5, colour="red", fill="#00CC00")+
  facet_wrap(~ variable)+
  labs(title="Fitted values for subject 95, from B=4000 nonparametric bootsrap replications", x=expression(paste("fitted value   ",hat(mu)[95])))+
  geom_vline(data=subset(dat.hist2.m, variable %in% c("BIC","Cp")), aes(xintercept=value,label=variable), colour="black", size=1)




ecdf(pros_SCAD$tstar.i)(pros_SCAD$X1[1]) 

#plot conf.intervals
p1 <- data.frame(lower=pros_lasso$l.stand[1], upper=pros_lasso$u.stand[1], type="standard", model="LASSO", estimate=pros_lasso$X1[1], length=pros_lasso$length[1], sd = pros_lasso$sdhat[1], coverage = pros_lasso$coverage[1])
p2 <- data.frame(lower=pros_lasso$l.quant[1], upper=pros_lasso$u.quant[1], type="quantile", model="LASSO", estimate=NA, length=pros_lasso$length.1[1], sd = NA, coverage = pros_lasso$coverage.1[1])
p3 <- data.frame(lower=pros_lasso$l.smooth[1],upper=pros_lasso$u.smooth[1],type="smooth", model="LASSO", estimate=pros_lasso$mutilde[1], length=pros_lasso$length.2[1], sd = pros_lasso$sdtilde[1], coverage = pros_lasso$coverage.2[1])

p4 <- data.frame(lower=pros_SCAD$l.stand[1], upper=pros_SCAD$u.stand[1], type="standard", model="SCAD", estimate=pros_SCAD$muhat[1], length=pros_SCAD$length[1], sd = pros_SCAD$sdhat[1], coverage = pros_SCAD$coverage[1])
p5 <- data.frame(lower=pros_SCAD$l.quant[1], upper=pros_SCAD$u.quant[1], type="quantile", model="SCAD", estimate=NA, length=pros_SCAD$length.1[1], sd = NA, coverage = pros_SCAD$coverage.1[1])
p6 <- data.frame(lower=pros_SCAD$l.smooth[1],upper=pros_SCAD$u.smooth[1],type="smooth", model="SCAD", estimate=pros_SCAD$mutilde[1], length=pros_SCAD$length.2[1], sd = pros_SCAD$sdtilde[1], coverage = pros_SCAD$coverage.2[1])

p7 <- data.frame(lower=pros_MCP$l.stand[1], upper=pros_MCP$u.stand[1], type="standard", model="MCP", estimate=pros_MCP$muhat[1], length=pros_MCP$length[1], sd = pros_MCP$sdhat[1], coverage = pros_MCP$coverage[1])
p8 <- data.frame(lower=pros_MCP$l.quant[1], upper=pros_MCP$u.quant[1], type="quantile", model="MCP", estimate=NA, length=pros_MCP$length.1[1], sd = NA, coverage = pros_MCP$coverage.1[1])
p9 <- data.frame(lower=pros_MCP$l.smooth[1],upper=pros_MCP$u.smooth[1],type="smooth", model="MCP", estimate=pros_MCP$mutilde[1], length=pros_MCP$length.2[1], sd = pros_MCP$sdtilde[1], coverage = pros_MCP$coverage.2[1])

p10 <- data.frame(lower=pros_bic$l.stand[1], upper=pros_bic$u.stand[1], type="standard", model="BIC", estimate=pros_bic$muhat[1], length=pros_bic$length[1], sd = pros_bic$sdhat[1], coverage = pros_bic$coverage[1])
p11 <- data.frame(lower=pros_bic$l.quant[1], upper=pros_bic$u.quant[1], type="quantile", model="BIC", estimate=NA, length=pros_bic$length.1[1], sd = NA, coverage = pros_bic$coverage.1[1])
p12 <- data.frame(lower=pros_bic$l.smooth[1],upper=pros_bic$u.smooth[1],type="smooth", model="BIC", estimate=pros_bic$mutilde[1], length=pros_bic$length.2[1], sd = pros_bic$sdtilde[1], coverage = pros_bic$coverage.2[1])

p13 <- data.frame(lower=pros_cp$l.stand[1], upper=pros_cp$u.stand[1], type="standard", model="Cp", estimate=pros_cp$muhat[1], length=pros_cp$length[1], sd = pros_cp$sdhat[1], coverage = pros_cp$coverage[1])
p14 <- data.frame(lower=pros_cp$l.quant[1], upper=pros_cp$u.quant[1], type="quantile", model="Cp", estimate=NA, length=pros_cp$length.1[1], sd = NA, coverage = pros_cp$coverage.1[1])
p15 <- data.frame(lower=pros_cp$l.smooth[1],upper=pros_cp$u.smooth[1],type="smooth", model="Cp", estimate=pros_cp$mutilde[1], length=pros_cp$length.2[1], sd = pros_cp$sdtilde[1], coverage = pros_cp$coverage.2[1])

p.all <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15)

library(xtable)
newdat<- p.all[,c("model","type","estimate","sd","length","coverage")];newdat
print(xtable(newdat), include.rownames=FALSE)


p.conf <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9)

#scad, MCP, lasso
pd <- position_dodge(width=0.3,height=NULL)
se1 <- ggplot(p.conf, aes(model, estimate, ymin = 0, ymax=length, colour = type)) + 
  geom_linerange(position=pd, width=0.3, size=2) +
  scale_y_continuous(limits=c(0, max(p.conf$length))) +
  coord_flip() +theme_bw() + 
  labs(x="penalty", y="length", title="Length of 95% Confidence Intervals, Subject 95 based \n on B=4000 nonparametric bootsrap samples for MCP, SCAD and LASSO 
       penalties")
se1 


p.conf2 <- rbind(p10,p11,p12,p13,p14,p15)

#BIC and CP 
pd <- position_dodge(width=0.3,height=NULL)
se2 <- ggplot(p.conf2, aes(model, estimate, ymin = 0, ymax=length, colour = type)) + 
  geom_linerange(position=pd, width=0.3, size=2) +
  scale_y_continuous(limits=c(0, max(p.conf2$length))) +
  coord_flip() +theme_bw() + 
  labs(x="penalty", y="length", title="Length of 95% Confidence Intervals for fitted value of Subject 95 based \n on B=4000 nonparametric bootsrap samples for Cp and BIC")
se2 



multiplot(se1,se2, cols=2)


