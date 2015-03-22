##################################
# R source code file for transforming Cholesterol data 
# as in Efron JASA 2014 paper
# 
# Created by Maxime, March 20
# Updated: 
# NOTE: 
##################################

library(bootstrap)

# z=compliance, y=improvement

fit <- lm(y ~ z + I(z^2) + I(z^3), data=cholost)
new.comp <- seq(0,100, length.out=1000)
new.val <- predict(object = fit, newdata = data.frame(z=new.comp))
plot(y ~ z, data=cholost, 
     xlab="Compliance", ylab="Improvement")
lines(x=new.comp, y=new.val, col="green")

#transformation 1 - scaled ranks

data <- cholost[order(cholost$z),]

data$transf <- qnorm((1:nrow(data) - 0.5)/nrow(data))

fit2 <- lm(y ~ transf + I(transf^2) + I(transf^3), data)
new.comp2 <- seq(min(data$transf),max(data$transf), length.out=1000)
new.val2 <- predict(object = fit2, newdata = data.frame(transf=new.comp2))
plot(y ~ transf, data, 
     xlab="Compliance", ylab="Improvement")
lines(x=new.comp2, y=new.val2, col="green")

#transformation 2 - accounting for ties

z.tr <- sapply(split(data$transf, data$z),mean)

data$transf2 <- NA
for(i in 1:nrow(data)){
  data$transf2[i] <- z.tr[names(z.tr) == as.character(data$z[i])]
}

# get data into form so that it works with Sy's functions
DT <- as.data.table(data[,c("y","transf2")])
setnames(DT,"transf2","x")

# make the transformation now. This is how my function was written up
DT <- transform(DT, x2=x^2, x3=x^3, x4=x^4, x5=x^5, x6=x^6)

# round the compliances to 5 digits, so that we can find the 11 subjects for Figure 5
DT[,x:=round(x,5)]

# Clean up environment
rm(list=setdiff(ls(),"DT"))

