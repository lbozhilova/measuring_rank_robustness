##################################################################
########### Measuring rank robustness in scored PINs #############
################# Robustness Analysis: auxiliary #################
##################################################################

# Created: 17/09/18
# Last edited: 19/09/18
# Edit history:
### 17/09/18: Copied from original scripts. 
### 19/09/18: Added NaN handling in continuity.score()

# Description:
# Functions for rank robustness analysis.
# rank.metrics() takes a data frame of metrics and calculates node ranks.
# overall.ranks() calculates overall ranks; the interval bounds should be set to lower.b=600,
#   upper.b=900 for the STRING and synthetic networks, and to lower.b=150, upper.b=280 for the 
#   HPRED network
# tr.sim() calculates the rank similarity between two sequences, rA, and rB, at some value k
# continuity.df() calculates rank similarity between consecutive thresholds (divided by a gap of th_gap), 
#   for given values of k, for all metrics in a rank data frame
# continuity.score() calculates rank continuity. Interval bounds should be set as with overall.ranks().
#   The score cut-off value is given by alpha (in this paper, alpha=0.90)
# identifiability.df() calculates the relaxed similarity (alpha=1.5) between the top n 
#   (n=100, or 20 in the SYN-GNP case) overall ranks, and regular threshold ranks.
# identifiability.score() calculates rank identifiability. Threshold range bounds need to be set as above.
# rank.range.df() calculates the rank ranges for each protein by each metric. Threshold range bounds 
#   need to be set as above.
# instability.score() calculates rank instability scores.

require("parallel")

#----- Threshold ranks -----#
rank.metrics <- function(metric.df){
  foo <- function(th, mdf){
    mdf <- subset(mdf, threshold==th)
    colid <- which(!(names(mdf) %in% c("name", "threshold")))
    rdf <- apply(mdf[,colid], 2, function(m) rank(m, na.last=FALSE, ties.method = "random"))
    rdf <- data.frame(name=mdf$name, threshold=mdf$threshold, rdf)
    return(rdf)
  }
  ths <- sort(unique(metric.df$threshold))
  rdf_list <- lapply(ths, function(th) foo(th, metric.df))
  rdf <- do.call("rbind", rdf_list)
  rownames(rdf) <- NULL
  return(rdf)
}

#----- Overall ranks ------#
overall.ranks <- function(rank.df, lower.b=150, upper.b=990){
  rank.df <- subset(rank.df, threshold>=lower.b & threshold<=upper.b)
  metrics <- colnames(rank.df); metrics <- metrics[!(metrics %in% c("name", "threshold"))]
  # Calculate mean ranks for each node
  aux.protname <- function(protname){
    rdf <- subset(rank.df, name==protname)[,metrics]
    meanR <- apply(rdf, 2, mean)
    meanR <- data.frame(name=protname, t(meanR))
    return(meanR)
  }
  protnames <- levels(rank.df$name)
  df <- lapply(protnames, aux.protname)
  df <- do.call("rbind", df)
  colnames(df) <- c("name", metrics)
  df <- as.data.frame(df)
  df[,metrics] <- apply(df[,metrics], 2, function(x) rank(x, na.last=TRUE, ties.method="random"))
  return(df)
}

#----- Trajanovski Rank Similarity -----#
tr.sim <- function(rA, rB, k){
  k <- k[k!=0]
  N <- length(rA)
  j <- floor(N*k)
  olap <- sapply(j, function(x) sum((rA>(N-x)) & (rB>(N-x))) )
  return(data.frame(k, j, olap, score=olap/j))
}

#----- Continuity (data frame) -----#
continuity.df <- function(rank.df, k=seq(0, 0.05, 0.001), th_gap=10){
  # Initialise
  lb <- seq(150, 990, 10); ub <- lb+th_gap; lb <- lb[ub<=990]; ub <- ub[ub<=990]
  metrics <- colnames(rank.df); metrics <- metrics[!(metrics %in% c("name", "threshold"))]
  input.df <- data.frame("lb"=lb, "ub"=ub, "metric"=rep(metrics, each=length(lb)))
  foo <- function(input){
    olap.df <- tr.sim(subset(rank.df, threshold==input$lb)[,as.character(input$metric)], subset(rank.df, threshold==input$ub)[,as.character(input$metric)], k)
    olap.df <- cbind("lb"=input$lb, "ub"=input$ub, "metric"=as.character(input$metric), olap.df)
    return(olap.df)
  }
  cl <- makeCluster(6, "FORK")
  olap_list <- parLapply(cl, 1:nrow(input.df), function(x) foo(input.df[x,]))
  stopCluster(cl)
  return(do.call("rbind", olap_list))
}

#----- Continuity (score) -----#
continuity.score <- function(cont.df, lower.b=150, upper.b=990, alpha=0.9){
  cont.df <- subset(cont.df, lb>=lower.b & ub <=upper.b & !(is.nan(score)))
  metrics <- levels(cont.df$metric)
  aux <- function(mc){
    c.df <- subset(cont.df, metric==mc)
    return(sum(c.df$score>=alpha)/nrow(c.df))
  }
  scores <- sapply(metrics, aux)
  names(scores) <- metrics
  return(scores)
}

#----- Identifiability (data frame) -----#
identifiability.df <- function(ranks.df, overall.ranks.df, n=100, alpha=1.5){
  N <- length(levels(overall.ranks.df$name))
  mcs <- colnames(overall.ranks.df)[-c(1, 27)] 
  bar <- function(mc){
    pnames <- as.character(overall.ranks.df[overall.ranks.df[,mc]>(N-n),"name"])
    ranks.df <- subset(ranks.df, name %in% pnames)
    ranks.df <- ranks.df[ranks.df[,mc]>(N-alpha*n),]
    if (nrow(ranks.df)==0){
      df <- data.frame(metric=mc, threshold=seq(150, 990, 10), score=0)
      return(df)
    }
    foo <- function(th){
      ranks.df <- subset(ranks.df, threshold==th)
      return(nrow(ranks.df)/n)
    }
    df <- data.frame(metric=mc, threshold=seq(150, 990, 10), score=sapply(seq(150, 990, 10), foo))
    return(df)
  }
  df <- do.call("rbind", lapply(mcs, bar))
  return(df)
}

#----- Identifiability (score) -----#
identifiability.score <- function(id.df, lower.b=150, upper.b=990){
  mcs <- levels(id.df$metric)
  id.df <- subset(id.df, threshold>=lower.b & threshold<=upper.b)
  scs <- sapply(mcs, function(m) min(subset(id.df, metric==m)$score))
  names(scs) <- mcs
  return(scs)
}

#----- Rank ranges -----#
rank.range.df  <- function(rank.df, lower.b=150, upper.b=990){
  rank.df <- subset(rank.df, threshold>=lower.b & threshold<=upper.b)
  aux.protname <- function(protname){
    rdf <- subset(rank.df, name==protname)
    rdf <- apply(rdf[,-(1:2)], 2, range)
    rdf <- rdf[2,]-rdf[1,]
    rdf <- data.frame(name=protname, metric=names(rdf), rr=unname(rdf))
    return(rdf)
  }
  protnames <- levels(rank.df$name)
  df <- lapply(protnames, aux.protname)
  df <- do.call("rbind", df)
  df <- cbind(df, rr_ratio=df$rr/length(protnames))
  return(df)
}

#----- Instability score -----#
instability.score <- function(rank.range.df, overall.ranks.df, k=0.01){
  mcs <- levels(rank.range.df$metric)
  N <- length(levels(overall.ranks.df$name))
  rranges <- numeric(length(mcs)); names(rranges) <- mcs
  for (mc in mcs){
    idxs <- which(overall.ranks.df[,mc]>(N*(1-k)))
    prots <- as.character(overall.ranks.df[idxs,]$name)
    rranges[mc] <- mean(subset(rank.range.df, name %in% prots & metric==mc)$rr_ratio)
  }
  return(rranges)
}
