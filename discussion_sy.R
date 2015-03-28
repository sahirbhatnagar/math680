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


