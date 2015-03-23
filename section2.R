##################################
# R source code file for section 2 of Efron JASA 2014 paper
# 
# Created by Sahir, Max, March 22
# Updated: 
# NOTE: 
##################################

rm(list=ls())

# import transformed data and functions -------------------------------------------------

# setwd("~/Biostats PhD/MATH680/Assignment 2/")
setwd("~/Dropbox/Winter 2015/MATH680/Assignments/A2/")
source("trans_choles_data.R")
source("functions.R")


# Bootstrap samples -------------------------------------------------------

B <- 4000
set.seed(4)
samples <- replicate(B,DT[sample(1:nrow(DT),replace=T),],simplify=F)

# Compliance for observation 1
obs <- -2.32316

# Chosen models based on bootstrap replicates
true.cp <- foreach(i = samples) %dopar% fit.best(i, method="CP", predict=obs)
df.cp <- ldply (true.cp, data.frame)

# how many times each model was chosen 
perc.chosen <- apply(df.cp[,-7], 2, sum)/B

# percentile interval for bootstrap replicates (also used for Figure 3)
ci.muhat <- data.frame(x=c(quantile(df.cp$muhat,c(0.025,0.975))[1],
               quantile(df.cp$muhat,c(0.025,0.975))[2]),y=c(0,0))


# Model selection on original data ----------------------------------------

fit.all <- fit.once(DT)
muhat.obs <-fit.all[which.min(fit.all$criteria),"muhat"]


# histogram of 4000 muhat bootstrap estimates: Figure 3
ggplot(df.cp, aes(muhat)) + geom_histogram(binwidth=1, colour="black", fill="white")+xlim(-31,31)+
    geom_vline(xintercept = muhat.obs, colour="red", linetype = "longdash",size=1.5)+
    annotate("text", label = round(muhat.obs,2), x = 5, y = -7, size = 5, colour = "red")+
    geom_point(data=ci.muhat, aes(x,y), size=4,shape=24, fill="red")+
    xlab("bootstrap estimates for subject 1")