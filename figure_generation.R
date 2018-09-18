#############################################################
######### Measuring rank robustness in scored PINs ##########
################ Figure and table generation ################
#############################################################

# Created: 17/09/18
# Last edited: 17/09/18
# Edit history:
### 17/09/18: Copied from original scripts. 

# Description:
# Figure and table generation for robustness analysis figures.
# Key data can be found in <continuityPlots.RData>, <identifiabilityPlots.RData>,
# and <imstabilityPlots.RData>.

source("plotPreamble.R")

load("F_string_pvivax_metrics.RData")
load("F_string_ecoli_metrics.RData")
load("F_string_yeast_metrics.RData")
load("F_hpred_yeast_metrics.RData")
load("F_gnp_syn_metrics.RData")
load("F_pvx_syn_metrics.RData")
chosenMets <- c("degree", "local_clustering", "betweenness", "NatConDiff")

#----- Rank continuity -----#
# Pull data for k=0.01 to plot
plot.continuity.df <- rbind(
  cbind(subset(string.pvivax.continuity.df, k==0.01), Network="PVX"),
  cbind(subset(string.ecoli.continuity.df, k==0.01), Network="ECOLI"),
  cbind(subset(string.yeast.continuity.df, k==0.01), Network="YEAST"),
  cbind(subset(hpred.yeast.continuity.df, k==0.01), Network="HPRED"),
  cbind(subset(gnp.syn.continuity.df, k==0.01), Network="SYN-GNP"),
  cbind(subset(pvx.syn.continuity.df, k==0.01), Network="SYN-PVX")
)
rm(string.pvivax.continuity.df, string.ecoli.continuity.df,pvx.syn.continuity.df,
   string.yeast.continuity.df, hpred.yeast.continuity.df, gnp.syn.continuity.df)

# Figure 2. Metric rank similarity between consecutive thresholds
df.pin <- subset(plot.continuity.df, metric %in% chosenMets & Network %in% c("PVX", "ECOLI", "YEAST", "HPRED"))
df.syn <- subset(plot.continuity.df, metric %in% chosenMets & Network %in% c("SYN-GNP", "SYN-PVX"))
fig2.pin <- ggplot(df.pin, aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
fig2.syn <- ggplot(df.syn, aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba_2, c("SYN-GNP", "SYN-PVX")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
plotDouble(fig2.pin, fig2.syn)

# Figure S1. Standard metrics, PINs
df.pin <- subset(plot.continuity.df, Network %in% c("PVX", "ECOLI", "YEAST", "HPRED"))
df.syn <- subset(plot.continuity.df, Network %in% c("SYN-GNP", "SYN-PVX"))
figs1 <-  ggplot(subset(df.pin, !(metric %in% metrics.loud)), aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
figs1

# Figure S2. LOUD metrics, PINs
figs2 <-  ggplot(subset(df.pin, metric %in% metrics.loud), aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
figs2

# Figure S3. Standard metrics, synthetic networks
figs3 <-  ggplot(subset(df.syn, !(metric %in% metrics.loud)), aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba_2, labels=c("SYN-GNP", "SYN-PVX")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
figs3

# Figure S4. LOUD metrics, synthetic networks
figs4 <-  ggplot(subset(df.syn, metric %in% metrics.loud), aes(x=lb*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba_2, labels=c("SYN-GNP", "SYN-PVX")) + 
  scale_linetype_manual(values=ltps) +
  geom_hline(yintercept=.90)
figs4

# Tables
continuity.score.matrix <- matrix(0, ncol=6, nrow=25)
rownames(continuity.score.matrix) <- names(metricNames)

colnames(continuity.score.matrix) <- c("PVX", "ECOLI", "YEAST", "HPRED","SYN-GNP", "SYN-PVX")
continuity.score.matrix[names(string.pvivax.continuity.scores), "PVX"] <- string.pvivax.continuity.scores
continuity.score.matrix[names(string.ecoli.continuity.scores), "ECOLI"] <- string.ecoli.continuity.scores
continuity.score.matrix[names(string.yeast.continuity.scores), "YEAST"] <- string.yeast.continuity.scores
continuity.score.matrix[names(hpred.yeast.continuity.scores), "HPRED"] <- hpred.yeast.continuity.scores
continuity.score.matrix[names(gnp.syn.continuity.scores), "SYN-GNP"] <- gnp.syn.continuity.scores
continuity.score.matrix[names(pvx.syn.continuity.scores), "SYN-PVX"] <- pvx.syn.continuity.scores
continuity.score.matrix

# Column Continuity in Table EV1
ev1.continuity <- round(apply(continuity.score.matrix[,1:4], 1, mean), digits=2)
ev1.continuity

# Table EV2
ev2 <- cor(continuity.score.matrix, method="spearman")
round(ev2, 2)

# Table S1
s1 <- round(continuity.score.matrix, 2)
s1

# Save and declutter
save(plot.continuity.df, continuity.score.matrix, file="continuityPlots.RData")
rm(continuity.score.matrix, df.pin, df.syn, ev2, fig2.pin, fig2.syn, figs1, figs2, figs3, figs4,
   string.pvivax.continuity.scores, string.ecoli.continuity.scores, string.yeast.continuity.scores,
   hpred.yeast.continuity.scores, gnp.syn.continuity.scores, pvx.syn.continuity.scores, ev1.continuity,
   s1); gc()

#----- Rank identifiability -----#
plot.identifiability.df <- rbind(
  cbind(string.pvivax.identifiability.df, Network="PVX"),
  cbind(string.ecoli.identifiability.df, Network="ECOLI"),
  cbind(string.yeast.identifiability.df, Network="YEAST"),
  cbind(hpred.yeast.identifiability.df, Network="HPRED"),
  cbind(gnp.syn.identifiability.df,  Network="SYN-GNP"),
  cbind(pvx.syn.identifiability.df,  Network="SYN-PVX")
)
rm(string.pvivax.identifiability.df, string.ecoli.identifiability.df,pvx.syn.identifiability.df,
   string.yeast.identifiability.df, hpred.yeast.identifiability.df, gnp.syn.identifiability.df)

# Figure 3: Relaxed similarity for some chosen metrics
df1 <- subset(plot.identifiability.df, metric %in% chosenMets & !(Network %in% c("SYN-GNP", "SYN-PVX")))
df2 <- subset(plot.identifiability.df, Network %in% c("SYN-GNP", "SYN-PVX") & metric %in% chosenMets)
df2$Network <- factor(df2$Network, levels=c("SYN-GNP", "SYN-PVX"))

fig3.pin <- ggplot(df1, aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.15, 0.28), linetype="dotted", colour=col_lyuba[5]) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black")
fig3.pin
fig3.syn <- ggplot(df2, aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black") +
  scale_colour_manual(name="Network", values = col_lyuba_2, labels=c("SYN-GNP", "SYN-PVX")) 
fig3.syn
plotDouble(fig3.pin, fig3.syn)

# Figure S5. Standard metrics, PINs
df1 <- subset(plot.identifiability.df, !(Network %in% c("SYN-GNP", "SYN-PVX")))
df2 <- subset(plot.identifiability.df, Network %in% c("SYN-GNP", "SYN-PVX"))
df2$Network <- factor(df2$Network, levels=c("SYN-GNP", "SYN-PVX"))

figs5 <- ggplot(subset(df1, !(metric %in% metrics.loud)), aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.15, 0.28), linetype="dotted", colour=col_lyuba[5]) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black")
figs5

# Figure S6. LOUD metrics, PINs
figs6 <- ggplot(subset(df1, (metric %in% metrics.loud)), aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  scale_colour_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.15, 0.28), linetype="dotted", colour=col_lyuba[5]) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black")
figs6

# Figure S7. Standard metrics, synthetic networks
figs7 <- ggplot(subset(df2, !(metric %in% metrics.loud)), aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black") +
  scale_colour_manual(name="Network", values = col_lyuba_2, labels=c("SYN-GNP", "SYN-PVX")) 
figs7

# Figure S8. Standard metrics, synthetic networks
figs8 <- ggplot(subset(df2, (metric %in% metrics.loud)), aes(x=threshold*.001, y=score, colour=Network)) + 
  theme_lyuba +
  geom_line(size=0.9) + facet_wrap(~metric, labeller=metric_labeller, ncol=3) +
  scale_x_continuous("Threshold", breaks=c(seq(0.15, 0.9, 0.15), 1), limits=c(0.15,1)) +
  scale_y_continuous("Similarity") +
  geom_hline(yintercept=.90) +
  geom_vline(xintercept=c(0.6, 0.9), linetype="dotted", colour="black") +
  scale_colour_manual(name="Network", values = col_lyuba_2, labels=c("SYN-GNP", "SYN-PVX"))
figs8

# Tables
identifiability.score.matrix <- matrix(0, ncol=6, nrow=25)
rownames(identifiability.score.matrix) <- names(metricNames)
colnames(identifiability.score.matrix) <- c("PVX", "ECOLI", "YEAST", "HPRED","SYN-GNP", "SYN-PVX")
identifiability.score.matrix[names(string.pvivax.identifiability.scores), "PVX"] <- string.pvivax.identifiability.scores
identifiability.score.matrix[names(string.ecoli.identifiability.scores), "ECOLI"] <- string.ecoli.identifiability.scores
identifiability.score.matrix[names(string.yeast.identifiability.scores), "YEAST"] <- string.yeast.identifiability.scores
identifiability.score.matrix[names(hpred.yeast.identifiability.scores), "HPRED"] <- hpred.yeast.identifiability.scores
identifiability.score.matrix[names(gnp.syn.identifiability.scores), "SYN-GNP"] <- gnp.syn.identifiability.scores
identifiability.score.matrix[names(pvx.syn.identifiability.scores), "SYN-PVX"] <- pvx.syn.identifiability.scores
identifiability.score.matrix

# Table EV1, Identifiability column
ev1.identifiability <- round(apply(identifiability.score.matrix[,1:4], 1, mean), digits=2)
ev1.identifiability

# Table EV2
ev3 <- cor(identifiability.score.matrix, method="spearman")
round(ev3, 2)

# Table S2
s2 <- round(identifiability.score.matrix, 2)
s2

save(plot.identifiability.df, identifiability.score.matrix, file="identifiabilityPlots.RData")
rm(df1, df2, fig3.pin, plot.identifiability.df, fig3.syn, figs5, figs6, figs7, figs8,
  string.pvivax.identifiability.scores, string.ecoli.identifiability.scores, string.yeast.identifiability.scores,
  hpred.yeast.identifiability.scores, gnp.syn.identifiability.scores, pvx.syn.identifiability.scores,
  ev1.identifiability, ev3, s2, identifiability.score.matrix)

#----- Rank instability -----#
instability.score.matrix <- matrix(0, ncol=6, nrow=25)
rownames(instability.score.matrix) <- names(metricNames)
colnames(instability.score.matrix) <- c("PVX", "ECOLI", "YEAST", "HPRED","SYN-GNP", "SYN-PVX")
instability.score.matrix[names(string.pvivax.instability.scores), "PVX"] <- string.pvivax.instability.scores
instability.score.matrix[names(string.ecoli.instability.scores), "ECOLI"] <- string.ecoli.instability.scores
instability.score.matrix[names(string.yeast.instability.scores), "YEAST"] <- string.yeast.instability.scores
instability.score.matrix[names(hpred.yeast.instability.scores), "HPRED"] <- hpred.yeast.instability.scores
instability.score.matrix[names(gnp.syn.instability.scores), "SYN-GNP"] <- gnp.syn.instability.scores
instability.score.matrix[names(pvx.syn.instability.scores), "SYN-PVX"] <- pvx.syn.instability.scores
instability.score.matrix

# Figure 4. Rank instability
df <- data.frame(metric=rownames(instability.score.matrix),
                 score=as.vector(instability.score.matrix),
                 Network=rep(colnames(instability.score.matrix), each=25))

df1 <- subset(df, metric %in% chosenMets & Network %in% c("PVX", "ECOLI", "YEAST", "HPRED"))
df1$Network <- factor(df1$Network, levels=c("PVX", "ECOLI", "YEAST", "HPRED"))
df2 <- subset(df, metric %in% chosenMets)
df2$Network <- factor(df2$Network, levels=c("PVX", "ECOLI", "YEAST", "HPRED", "SYN-GNP", "SYN-PVX"))

fig4.pin <- ggplot(df1, aes(x=Network, y=score, fill=Network)) + 
  theme_lyuba + 
  geom_col() + 
  facet_wrap(~metric, labeller=metric_labeller) +
  scale_y_continuous("Instability") +
  scale_fill_manual(name="Network", values = col_lyuba, labels=c("PVX", "ECOLI", "YEAST", "HPRED")) + 
  geom_hline(yintercept = 0.01, linetype="dotted")
fig4.pin
fig4.syn <- ggplot(df2, aes(x=Network, y=score, fill=Network)) + 
  theme_lyuba + 
  geom_col() + 
  facet_wrap(~metric, labeller=metric_labeller) +
  scale_y_continuous("Instability",  breaks=seq(0, 1, 0.25), limits=c(0,1) ) +
  scale_fill_manual(name="Network", values = c(col_lyuba, col_lyuba_2), labels=c("PVX", "ECOLI", "YEAST", "HPRED", "SYN-GNP", "SYN-PVX")) +
  geom_hline(yintercept = 0.01, linetype="dotted") 
fig4.syn
plotDouble(fig4.pin, fig4.syn)

# Table EV1, Instability column
ev1.instability <- round(apply(instability.score.matrix[,1:4], 1, mean), digits=2)
ev1.instability

# Table EV4
ev4 <- cor(instability.score.matrix, method="spearman")
round(ev4, 2)

# Table S3
s3 <- round(instability.score.matrix, 2)
s3

save(df, instability.score.matrix, file="instabilityPlots.RData")
