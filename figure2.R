##########################
# Figure 2 code - Grinux #
##########################

i <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

source("trans_choles_data.R")

B <- 4000
samples <- replicate(B,DT[sample(1:nrow(DT),replace=T),],simplify=F)

obs <- DT[i,x]

source("functions.R")
jj <- fit.all(data=DT,pred=obs, bootsamples = samples, B)

write.table(c(obs, jj), file=paste0("results", i, ".txt"), 
            row.names=FALSE, col.names=FALSE)