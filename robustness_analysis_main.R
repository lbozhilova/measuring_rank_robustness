#############################################################
######### Measuring rank robustness in scored PINs ##########
################# Robustness Analysis: main #################
#############################################################

# Created: 17/09/18
# Last edited: 17/09/18
# Edit history:
### 17/09/18: Copied from original scripts. 

# Description:
# Robustness analysis for the SYN-GNP network. Other networks, and their associated .RData files, follow the
# same general framework.

source("robustness_analysis_aux.R")

#----- Obtain threshold and overall ranks -----#
load("F_pvx_syn_metrics.RData")
pvx.syn.ranks.df <- rank.metrics(pvx.syn.metrics.df)
pvx.syn.overall.ranks.df <- overall.ranks(pvx.syn.ranks.df, lower.b=600, upper.b=900)

#----- Rank continuity -----#
pvx.syn.continuity.df <- continuity.df(pvx.syn.ranks.df)
pvx.syn.continuity.scores <- continuity.score(pvx.syn.continuity.df, lower.b=600, upper.b=900, alpha=0.90)

#----- Rank identifiability -----#
pvx.syn.identifiability.df <- identifiability.df(pvx.syn.ranks.df, pvx.syn.overall.ranks.df, n=100, alpha=1.5)
pvx.syn.identifiability.scores <- identifiability.score(pvx.syn.identifiability.df, lower.b=600, upper.b=900)

#----- Rank instability -----#
pvx.syn.rank.range.df <- rank.range.df(pvx.syn.ranks.df, lower.b=600, upper.b=900)
pvx.syn.instability.scores <- instability.score(pvx.syn.rank.range.df, pvx.syn.overall.ranks.df, k=0.01)

save(pvx.syn.metrics.df, pvx.syn.ranks.df, pvx.syn.overall.ranks.df, pvx.syn.continuity.df,
     pvx.syn.continuity.scores, pvx.syn.identifiability.df, pvx.syn.identifiability.scores,
     pvx.syn.rank.range.df, pvx.syn.instability.scores, file= "F_pvx_syn_metrics.RData")
