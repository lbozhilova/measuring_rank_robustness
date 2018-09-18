#########################################
##### Plot basics: packages & theme #####
#########################################

#--- Packages ---#
require("ggplot2")
require("ggthemes")
require("grid")
require("gridExtra")
require("igraph")
#----------------#

#--- Theme ---#
theme_lyuba <- theme_minimal() +
  theme(text = element_text(color = "gray20"),
        legend.position = c("top"), # position the legend in the upper left 
        legend.direction = "horizontal",
        legend.justification = 0.1, # anchor point for legend.position.
        legend.text = element_text(size = 12, color = "gray10"),
        axis.text = element_text(face = "italic"),
        axis.title.x = element_text(vjust = -1), # move title away from axis
        axis.title.y = element_text(vjust = 2), # move away for axis
        axis.ticks.x = element_line(color="gray40", size=0.3),
        axis.ticks.y = element_blank(), # element_blank() is how we remove elements
        axis.line = element_line(color = "gray40", size = 0.3),
        axis.line.y = element_blank(),
        panel.grid.major = element_line(color = "gray50", size = 0.5),
        panel.grid.major.x = element_blank()
  )
col_lyuba <- c("#00274c",  "#ffc300", "#af2a35", "#ab7895")
col_lyuba_2 <- c("coral", "turquoise4")
metricNames <- c("avg_btw"="Avg. betweenness", "avg_closeness"="Avg. closeness", "avg_e_ego_one"="Avg. e_one",
                 "avg_local_c"="Avg. local C", "avg_n_ego_diff"="Avg. n_diff", "avg_n_ego_ratio"="Avg. n_ratio",
                 "avg_n_ego_sqdiff"="Avg. n_sqdiff", "avg_n_ego_two"="Avg. n_two", "avg_path"="Avg. shortest path",
                 "avg_redundancy"="Avg. redundancy", "betweenness"="Betweenness", "closeness"="Closeness",
                 "closeness_MEJN"="Harmonic c-ty", "degree"="Degree", "e_ego_one"="e_one", "global_c"="Global C",
                 "local_clustering"="Local C", "n_ego_diff"="n_diff", "n_ego_ratio"="n_ratio", "n_ego_sqdiff"="n_sqdiff",
                 "n_ego_two"="n_two", "NatConDiff" = "Natural connectivity", "num_con_pairs"="Number of pairs",
                 "pagerank"="PageRank", "redundancy"="Redundancy")
metrics.loud <- c("avg_btw", "avg_closeness", "avg_e_ego_one",
                  "avg_local_c", "avg_n_ego_diff", "avg_n_ego_ratio",
                  "avg_n_ego_sqdiff", "avg_n_ego_two", "avg_path",
                  "avg_redundancy", "global_c", "num_con_pairs",
                  "NatConDiff")
metric_labeller <- labeller(metric=metricNames)
#-------------#

#--- Plot functions ---#
plotDouble <- function(myplot1, myplot2, nc=2){
  myplot1 <- arrangeGrob(myplot1, top = textGrob("A", x = unit(0, "npc")
                                                 , y   = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot2 <- arrangeGrob(myplot2, top = textGrob("B", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  grid.arrange(myplot1, myplot2, ncol = nc)
}
plotQuad <- function(myplot1, myplot2, myplot3, myplot4, nc=2){
  myplot1 <- arrangeGrob(myplot1, top = textGrob("A", x = unit(0, "npc")
                                                 , y   = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot2 <- arrangeGrob(myplot2, top = textGrob("B", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot3 <- arrangeGrob(myplot3, top = textGrob("C", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot4 <- arrangeGrob(myplot4, top = textGrob("D", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  grid.arrange(myplot1, myplot2, myplot3, myplot4, ncol = nc)
}
plotFive <- function(myplot1, myplot2, myplot3, myplot4, myplot5, nc=2){
  myplot1 <- arrangeGrob(myplot1, top = textGrob("A", x = unit(0, "npc")
                                                 , y   = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot2 <- arrangeGrob(myplot2, top = textGrob("B", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot3 <- arrangeGrob(myplot3, top = textGrob("C", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot4 <- arrangeGrob(myplot4, top = textGrob("D", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  myplot5 <- arrangeGrob(myplot5, top = textGrob("E", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
  grid.arrange(myplot1, myplot2, myplot3, myplot4, myplot5, ncol = nc)
}
#----------------------#