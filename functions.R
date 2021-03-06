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
    
    #     # all possible models
    #     Cols <- names(DF)
    #     Cols <- Cols[! Cols %in% "y" & ! Cols %in% "x0"]
    #     p <- length(Cols)
    #     id <- unlist(lapply(1:p, function(i) combn(1:p,i,simplify=F)), recursive=F)
    #     formulas <- sapply(id,function(i) paste("y~",paste(Cols[i],collapse="+")))
    
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
        
#         #model selection on original data
#         fit.las <- cv.glmnet(x=as.matrix(x),y=y, nfolds=5)
#         
#         # \hat{\mu} = t(y) , predict everyone
#         t.y <- predict(fit.las,as.matrix(pred), s="lambda.min")
#         
#         # Bootstrap samples for simulation 
#         samples.las <- replicate(B,data[sample(1:nrow(data),replace=T),],simplify=F)
#         
#         # t_i^*(y_i^*), i.e., the predicted values for each individual j=1,...,n, for each 
#         # bootstrap sample i=1,...,B, so we have 200 predicted values for each bootstrap sample.
#         tstar.i <- foreach(i = samples.las, .combine='cbind') %dopar% pred.las(i, pred, n.lam=n.lambda)
#         
#         # s(y) = B^{-1} \sumn t_i^*(y_i^*) for each j=1,...,n
#         # this is the mean for each predicted value across each bootstrap sample
#         s.y.j <- apply(tstar.i,1,mean)
#         # \hat{sd_B}, the standard bootstrap standard error
#         sd.hat <- apply(tstar.i,1,sd)
#         
#         # Standard interval ----
#         stand.int <- cbind(t.y+qnorm(0.025)*sd.hat,t.y+qnorm(0.975)*sd.hat)
#         #center, length and coverage
#         
#         m <- lapply(1:n, function(j) vals(stand.int[j,], tstar.i[j,], n.boot=B))
#         result <- ldply (m, data.frame)
#         result <- transform(result.stand, muhat=t.y, mutilde=s.y.j, sdhat=sd.hat, lower=stand.int[,1], upper=stand.int[,2])
         
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