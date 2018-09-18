################################################################
########### Measuring rank robustness in scored PINs ###########
##### Data pre-processing and synthetic network generation #####
################################################################

# Created: 17/09/18
# Last edited: 17/09/18
# Edit history:
### 17/09/18: Copied from original scripts. 

# Description: Script for data pre-processing (for the three STRING and one HitPredict networks), and synthetic
# network generation (for both SYN-GNP and SYN-PVX).
# Each network section is self-contained.
# Where files are written to or accessed, the exact path name has been removed.

#----- Packages -----#
require("stringr")
require("igraph")
#--------------------#

#----- YEAST network -----#
# STRING v.10.5, accessed on 19/02/18
string.yeast <- read.table("4932.protein.links.v10.5.txt.gz", header=TRUE)
df <- data.frame(protein1=str_replace(string.yeast$protein1, "4932.", ""),
                 protein2=str_replace(string.yeast$protein2, "4932.", ""),
                 combined_score=string.yeast$combined_score)
g <- graph_from_data_frame(df, directed=FALSE)
# Simplify and check edges have the same weight in both directions.
g <- simplify(g, remove.multiple = TRUE, remove.loops=TRUE, edge.attr.comb = "concat")
x <- E(g)$combined_score
foo <- function(z){
  if(z[1]==z[2]) return(z[1])
  else stop("Edge score at odds.")
}
x <- unlist(lapply(x, foo)) # foo() doesn't throw any errors
g <- delete_edge_attr(g, "combined_score")
g <- set.edge.attribute(g, "combined_score", index=E(g), value=x)
string.yeast.df <- as_data_frame(g, "edges")
colnames(string.yeast.df) <- c("protein1", "protein2", "combined_score")
string.yeast.net <- g
save(string.yeast.df, string.yeast.net, file="string_yeast_data.RData")
rm(list=ls()); gc()

#----- ECOLI network -----#
# STRING v.10.5, accessed on 19/02/18
string.ecoli <- read.table("511145.protein.links.v10.5.txt.gz", header=TRUE)
df <- data.frame(protein1=str_replace(string.ecoli$protein1, "511145.", ""),
                 protein2=str_replace(string.ecoli$protein2, "511145.", ""),
                 combined_score=string.ecoli$combined_score)
g <- graph_from_data_frame(df, directed=FALSE)
# Simplify and check edges have the same weight in both directions.
g <- simplify(g, remove.multiple = TRUE, remove.loops=TRUE, edge.attr.comb = "concat")
x <- E(g)$combined_score
foo <- function(z){
  if(z[1]==z[2]) return(z[1])
  else stop("Edge score at odds.")
}
x <- unlist(lapply(x, foo))
g <- delete_edge_attr(g, "combined_score")
g <- set.edge.attribute(g, "combined_score", index=E(g), value=x)
string.ecoli.df <- as_data_frame(g, "edges")
colnames(string.ecoli.df) <- c("protein1", "protein2", "combined_score")
string.ecoli.net <- g
save(string.ecoli.df, string.ecoli.net, file="string_ecoli_data.RData")
rm(list=ls()); gc()

#----- PVX network -----#
# STRING v.10.5, accessed on 12/03/18
string.pvivax <- read.table("5855.protein.links.v10.5.txt.gz", header=TRUE)
df <- data.frame(protein1=str_replace(string.pvivax$protein1, "5855.", ""),
                 protein2=str_replace(string.pvivax$protein2, "5855.", ""),
                 combined_score=string.pvivax$combined_score)
g <- graph_from_data_frame(df, directed=FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops=TRUE, edge.attr.comb = "concat")
x <- E(g)$combined_score
foo <- function(z){
  if(z[1]==z[2]) return(z[1])
  else stop("Edge score at odds")
}
x <- unlist(lapply(x, foo))
g <- delete_edge_attr(g, "combined_score")
g <- set.edge.attribute(g, "combined_score", index=E(g), value=x)
string.pvivax.df <- as_data_frame(g, "edges")
colnames(string.pvivax.df) <- c("protein1", "protein2", "combined_score")
string.pvivax.net <- g
save(string.pvivax.df, string.pvivax.net, file="string_pvivax_data.RData")
rm(list=ls()); gc()

#----- HPRED network -----#
# HitPredict, current release as of 19/02/18. 
# Manually removed preamble from text file and Ensembl column names.
hpred.yeast <- read.table("hitpredict_S_cerevisiae_interactions.txt", header=TRUE, fill=TRUE) # 116083
hpred.yeast <- subset(hpred.yeast, Taxonomy==559292) # 115005
df <- hpred.yeast[,c(1, 2, 12)]
colnames(df) <- c("protein1", "protein2", "combined_score")
df[,3] <- as.numeric(as.character(df[,3]))
x <- df[,3] <- round(df[,3]*1000)
g <- graph_from_data_frame(df, directed=FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops=TRUE, edge.attr.comb = "concat")
hpred.yeast.df <- as_data_frame(g, "edges")
colnames(hpred.yeast.df) <- c("protein1", "protein2", "combined_score")
hpred.yeast.net <- g
save(hpred.yeast.df, hpred.yeast.net, file="hpred_yeast_data.RData")
rm(list=ls()); gc()

#----- SYN-GNP network -----#
# Create a GNP network on 500 nodes with edge probability p=0.06, and assign random scores
# from the PVX network
load("string_pvivax_data.RData")
set.seed(73)
gnp.syn.net <- erdos.renyi.game(n=500, p.or.m=0.06)
gnp.syn.net <- set.vertex.attribute(gnp.syn.net, "name", value=V(string.pvivax.net)$name[1:500])
gnp.syn.net <- set.edge.attribute(gnp.syn.net, "combined_score", 
                                  value=sample(scores, size=ecount(gnp.syn.net), replace=TRUE))
save(gnp.syn.net, file="gnp_syn_data.RData")
rm(list=ls()); gc()

#----- SYN-PVX network -----#
# Take 1000 nodes from the PVX network and shuffle the scores.
set.seed(363)
load("string_pvivax_data.RData")
N <- vcount(string.pvivax.net)
ids <- sample.int(N, size=1000, replace=F)
pvx.syn.net <- induced.subgraph(string.pvivax.net, vids=ids)
scores <- E(pvx.syn.net)$combined_score
E(pvx.syn.net)$combined_score <- sample(scores)
save(pvx.syn.net, file="pvx_syn_data.RData")
rm(list=ls()); gc()