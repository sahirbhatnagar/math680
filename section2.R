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

#source("functions.R")
#source("trans_choles_data.R")

## ---- bootstrap-samples ----
B <- 4000
samples <- replicate(B,DT[sample(1:nrow(DT),replace=T),],simplify=F)


## ---- Cp-boot ----
source("functions.R")

obs <- -2.32316

# Chosen models based on bootstrap replicates
true.cp <- foreach(i = samples) %dopar% fit.best(i, method="CP", predict=obs)
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
                   m=fit.all$p, Cp=round(fit.all$criteria,0), "Bootrap"=perc.chosen)

# table 2 ----
tab2 <- data.frame(rbind(sapply(paste0("m",1:6), function(i) df.cp[which(get(i)),round(mean(muhat),2)]),
sapply(paste0("m",1:6), function(i) df.cp[which(get(i)),round(sd(muhat),2)])))
rownames(tab2) <- c("Mean","St.dev.")


## ---- figure-3 ----
# histogram of 4000 muhat bootstrap estimates: Figure 3
ggplot(df.cp, aes(muhat)) + geom_histogram(binwidth=1, colour="black", fill="white")+xlim(-31,31)+
    geom_vline(xintercept = muhat.obs, colour="red", linetype = "longdash",size=1.5)+
    annotate("text", label = round(muhat.obs,2), x = 5, y = -7, size = 5, colour = "red")+
    geom_point(data=ci.muhat, aes(x,y), size=4,shape=24, fill="red")+
    xlab("bootstrap estimates for subject 1")