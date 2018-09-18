################################################################
########## Measuring rank robustness in scored PINs ############
################# Metric extraction: auxiliary #################
################################################################

# Created: 17/09/18
# Last edited: 17/09/18
# Edit history:
### 17/09/18: Copied from original scripts. 

# Description:
# Metric extraction functions. Each function takes a network and a set of thresholds as input and 
# outputs a dataframe of node metrics. Optional parameters include csize (Int) for the number of 
# clusters to use, and periphery (Bool) for whether to do the evaluation on nodes of degree 0 or 1.

# node.metrics.standard() calculates the standard node metrics
# node.metrics.loo() calculates LOUD metrics. Since this can be much slower (remove NatConDiff for a
# speed-up), intermediate results are saved under "loo.metrics.th.NUM.RData", where NUM=current threshold.

#---- Packages -----#
require("igraph")
require("Matrix")
require("Brobdingnag")
require("parallel")
#-------------------#

node.metrics.standard <- function(g, ths, csize=6){
  foo <- function(th){
    g <- delete.edges(g, E(g)[combined_score<th])
    degree <- degree(g, v=vids) 
    local_clustering <- transitivity(g, type="local", vids=vids) 
    local_clustering[is.nan(local_clustering)] <- 0
    redundancy <- local_clustering*(degree-1) 
    e_ego_one <- degree + local_clustering*choose(degree, 2)
    n_ego_two <- ego_size(g, order=2, nodes=vids)
    n_ego_diff <- n_ego_two - degree - 1
    n_ego_sqdiff <- (degree+1)^2-n_ego_two
    n_ego_ratio <- (degree+1)/n_ego_two
    pagerank <- try(page_rank(g, vids=vids, directed=FALSE, damping=.85), silent=T) #
    if (class(pagerank)=="try-error")
      pagerank <- rep(NA, length(vids))
    else 
      pagerank <- unname(pagerank[[1]][vids])
    distance_matrix <- distances(g, v=vids)
    closeness <- sapply(1:vcount(g), function(vid) { 
      x <- distance_matrix[vid,-vid]; 
      x[is.infinite(x)] <- vcount(g); 
      return(1/mean(x))})
    closeness_MEJN <- 1/(vcount(g)-1)*sapply(1:vcount(g), function(vid) {sum(1/distance_matrix[vid,-vid])}) 
    betweenness <- betweenness(g, v=vids, directed=FALSE) 
    df <- data.frame(name=vids, threshold=th, degree, local_clustering, redundancy, e_ego_one,
                     n_ego_two, n_ego_diff, n_ego_sqdiff, n_ego_ratio, pagerank, closeness,
                     closeness_MEJN, betweenness)
    rownames(df) <- NULL
    return(df)
  }
  vids <- V(g)$name
  cl <- makeCluster(csize, "FORK")
  df_list <- parLapply(cl, ths, foo)
  stopCluster(cl)
  df <- do.call("rbind", df_list)
  return(df)
}

node.metrics.loo <- function(g, ths, csize=6, periphery=TRUE){
  foo <- function(th){
    bar <- function(vid){
      Q <- get.adjacency(g)
      Q[vid,] <- Q[,vid] <- 0
      g <- graph_from_adjacency_matrix(Q)
      adj_eigs <- eigen(Q, symmetric=TRUE, only.values=TRUE)$values
      NatConDiff <- log(sum(as.brob(exp(adj_eigs))))
      degs <- degree(g)
      local.clustering <- transitivity(g, "local")
      local.clustering[is.na(local.clustering)] <- 0
      avg_local_c <- mean(local.clustering)
      global_c <- transitivity(g, "global")
      avg_redundancy <- mean(local.clustering*(degs-1))
      avg_e_ego_one <- mean(degs + local.clustering*choose(degs, 2))
      nodes.ego.two <- ego_size(g, order=2)
      avg_n_ego_two <- mean(nodes.ego.two)
      avg_n_ego_diff <- mean(nodes.ego.two - degs -1)
      avg_n_ego_sqdiff <- mean((degs+1)^2-nodes.ego.two)
      avg_n_ego_ratio <- mean((degs+1)/nodes.ego.two)
      csize <- components(g)$csize
      num_con_pairs <- sum(choose(csize, 2))
      avg_path <- mean_distance(g, directed = FALSE, unconnected = TRUE)
      avg_btw <- num_con_pairs*(avg_path-1)
      avg_closeness <- 1/avg_path
      df <- data.frame(name=vid, NatConDiff, avg_local_c, global_c, avg_redundancy,
                       avg_e_ego_one, avg_n_ego_two, avg_n_ego_diff, avg_n_ego_sqdiff, avg_n_ego_ratio,
                       num_con_pairs, avg_path, avg_btw, avg_closeness)
      return(df)
    }
    g <- delete.edges(g, E(g)[combined_score<th])
    if(!periphery)
      vids <- V(g)[degree(g)>1]$name
    else
      vids <- V(g)$name
    cl <- makeCluster(csize, "FORK")
    df_list <- parLapply(cl, vids, bar)
    stopCluster(cl)
    xdf <- do.call("rbind", df_list)
    xdf <- cbind(xdf, threshold=th)
    df <- cbind(threshold=xdf[,15], xdf[,-15])
    if(!periphery){
      vids <- setdiff(V(g)$name, vids)
      if (length(vids)>0)
        df <- rbind(df, data.frame(threshold=th, name=vids, NatConDiff=NA, avg_local_c=NA, global_c=NA, avg_redundancy=NA,
                                   avg_e_ego_one=NA, avg_n_ego_two=NA, avg_n_ego_diff=NA, avg_n_ego_sqdiff=NA, avg_n_ego_ratio=NA,
                                   num_con_pairs=NA, avg_path=NA, avg_btw=NA, avg_closeness=NA))
    }
    save(df, file=paste0("loo.metrics.th.", th, ".RData"))
    return(df)
  }
  df_list <- lapply(ths, foo)
  df <- do.call("rbind", df_list)
  df[,-(1:2)] <- -df[,-(1:2)]
  return(df)
}
