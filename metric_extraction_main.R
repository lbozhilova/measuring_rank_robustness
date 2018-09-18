###########################################################
######## Measuring rank robustness in scored PINs #########
################# Metric extraction: main #################
###########################################################

# Created: 17/09/18
# Last edited: 17/09/18
# Edit history:
### 17/09/18: Copied from original script. 

# Description:
# Main code for metric extraction. Due to computational cost and having to share server time,
# this was done across a number of OPIG servers, and results were saved in batches and then combined. 
# Here we illustrate with SYN-PVX, although the application to other networks is the same. The separate batches
# are called "pvx_syn_A.RData", "pvx_syn_B.RData", "pvx_syn_C.RData". The dataframes in these were then combined,
# and saved under "F_pvx_syn_metrics.RData". The remaining networks were treated analogously. 

# If you want to use this for your own network, I recommend either removing the NatConDiff calculation 
# in <metric_extraction_aux.R>, or using a very big computer and limiting yourself to networks on 
# ~6,000 nodes or less.

#--- Packages ---#
library("igraph")
library("parallel")
library("Matrix")
library("Brobdingnag")
#----------------#

#----- Metric extraction -----#
rm(list=ls()); gc()
load("pvx_syn_data.RData")
source("metric_extraction_aux.R")

ths <- seq(150, 990, 10)
csize = 40

pvx.syn.standard.df <- node.metrics.standard(pvx.syn.net, ths[1:30], csize)
pvx.syn.loo.df <- node.metrics.loo(pvx.syn.net, ths[1:30], csize, periphery = FALSE)
save(pvx.syn.standard.df, pvx.syn.loo.df, file="pvx_syn_A.RData")
rm(pvx.syn.standard.df, pvx.syn.loo.df); gc()

pvx.syn.standard.df <- node.metrics.standard(pvx.syn.net, ths[31:60], csize)
pvx.syn.loo.df <- node.metrics.loo(pvx.syn.net, ths[31:60], csize, periphery = FALSE)
save(pvx.syn.standard.df, pvx.syn.loo.df, file="pvx_syn_B.RData")
rm(pvx.syn.standard.df, pvx.syn.loo.df); gc()

pvx.syn.standard.df <- node.metrics.standard(pvx.syn.net, ths[61:85], csize)
pvx.syn.loo.df <- node.metrics.loo(pvx.syn.net, ths[61:85], csize, periphery = FALSE)
save(pvx.syn.standard.df, pvx.syn.loo.df, file="pvx_syn_C.RData")
rm(pvx.syn.standard.df, pvx.syn.loo.df); gc()

rm(list=ls()); gc()

#----- Merging data -----#
load("pvx_syn_A.RData")
m1 <- merge(pvx.syn.standard.df, pvx.syn.loo.df, by=c("name", "threshold")) # merge node metric types
load("pvx_syn_B.RData")
m2 <- merge(pvx.syn.standard.df, pvx.syn.loo.df, by=c("name", "threshold")) # note these DFs are now different
load("pvx_syn_C.RData") 
m3 <- merge(pvx.syn.standard.df, pvx.syn.loo.df, by=c("name", "threshold")) # note these DFs are now different

pvx.syn.metrics.df <- rbind(m1, m2, m3)
save(pvx.syn.metrics.df, file="F_pvx_syn_metrics.RData")