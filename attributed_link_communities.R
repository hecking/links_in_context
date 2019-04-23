require(igraph)
require(dplyr)
require(tidyr)
require(Matrix)
require(magrittr)


getLinkSimilarityNetFeat <- function(anchorNode, endpoint1, endpoint2, egos, features) {

  relevantAttributes1 <- which((features[anchorNode,] * features[endpoint1,]) == 1)
  relevantAttributes2 <- which((features[anchorNode,] * features[endpoint2,]) == 1)
  commonRelevantAttributes <- intersect(relevantAttributes1, relevantAttributes2)
  if (length(commonRelevantAttributes) > 0) {

    relevantAttributes <- union(relevantAttributes1, relevantAttributes2)

    nb1 <- egos[[endpoint1]]$name
    nb2 <- egos[[endpoint2]]$name
    commonNbs <- intersect(nb1, nb2)

    sum(features[commonNbs,commonRelevantAttributes]) /
      (length(relevantAttributes) * length(union(nb1, nb2)))

  } else {
    0
  }
}

getLinkSimilarityNetwork <- function(anchorNode, endpoint1, endpoint2, egos, features=NULL, ...) {

  nb1 <- egos[[endpoint1]]
  nb2 <- egos[[endpoint2]]

  length(intersect(nb1, nb2)) / length(union(nb1, nb2))
}

getLinkSimilarityFeature <- function(anchorNode, endpoint1, endpoint2, g=NULL, features, ...) {

  relevantAttributes1 <- which((features[anchorNode,] * features[endpoint1,]) == 1)
  relevantAttributes2 <- which((features[anchorNode,] * features[endpoint2,]) == 1)
  ifelse(((length(relevantAttributes1) > 0) && (length(relevantAttributes2) > 0)),
    length(intersect(relevantAttributes1, relevantAttributes2)) / length(union(relevantAttributes1, relevantAttributes2)),
    0
  )
}

getLinkSimilarity <- function(g, features, simFun, ...) {

  endpoints <- igraph::as_data_frame(g)
  endpoints$eId <- 1:nrow(endpoints)
  endpoints <- bind_rows(endpoints,
                         data_frame(
                           from=endpoints$to,
                           to=endpoints$from,
                           eId=endpoints$eId)) # similarity calculation for both anchor points.
  nbs <- ego(g, 1)
  names(nbs) <- V(g)$name

  endpoints %>%
    group_by(from) %>%
    do({

      if (nrow(.) > 1) {

        edgePairs <- combn(1:nrow(.), 2)

        sim <- sapply(1:ncol(edgePairs), function(i) {

          simFun(.$from[1], .$to[edgePairs[1,i]], .$to[edgePairs[2,i]], nbs, features, ...)
        })


        data_frame(e1=.$eId[edgePairs[1,]], e2=.$eId[edgePairs[2,]], similarity=sim)
      } else {

        data_frame()
      }
    }) %>% ungroup %>% select(-from) # remove grouping from identifier.
}

partitionDensity <- function(g, clusters, ...) {

  sapply(clusters, function(cl) {

    if (length(cl) > 1) {
      sg <- subgraph.edges(g, cl)

      ecount(sg) * ((ecount(sg) - (vcount(sg) - 1))) /
        ((vcount(sg) - 2) * (vcount(sg) - 1))
    } else {
      0
    }
  }) %>% sum * 2 / ecount(g)
}

attributedPartitionDensity <- function(g, clusters, features, numFeatures) {

  m <- ((features[V(g)$name,] %*% t(features[V(g)$name,])) * as_adj(g)) %>% sum / 2

   # number of pairwise shared attributes
  sapply(clusters, function(cl) {

    if (length(cl) > 1) {
      sg <- subgraph.edges(g, cl)
      mc <- ((features[V(sg)$name,] %*% t(features[V(sg)$name,])) * as_adj(sg)) %>% sum / 2

      # mc * ((mc - (vcount(sg) - 1))) /
      #   ((numFeatures * vcount(sg) - 2) * (vcount(sg) - 1))
      mc * mc /
        ((numFeatures * vcount(sg)) * (vcount(sg) - 1))
    } else {
      0
    }
  }) %>% sum * 2 / m
}

# Parameters: Dendogram object, graph object, feature matrix (one hot coding), number of different feature type.
cutLinkDendogram <- function(dendogram, g, features, numFeatures) {

  cuts <- seq(0.2, max(dendogram$height), 0.05) 
  
  #cutree(dendogram, k=cuts) %>%
  cutree(dendogram, h = cuts) %>%
    dplyr::as_data_frame() %>%
    summarise_all(function(linkMembership) {
      lapply(unique(linkMembership), function(linkCluster) {

        which(linkMembership == linkCluster)
      }) %>% list
    }) %>%
    gather(key="cut", value="link_clusters") %>%
    rowwise %>%
    do({
      # Change from 10/12/2018
      
      data_frame(cut=.$cut, link_clusters=list(.$link_clusters),
                 num_clusters=length(.$link_clusters),
                 pd=partitionDensity(g, .$link_clusters),
                 apd=attributedPartitionDensity(g, .$link_clusters, features, numFeatures))
    })
}

removeNoContextEdges <- function(g, features) {

  edges <- igraph::as_data_frame(g)
  noContextEdges <- which(rowSums(features[edges$from,] * features[edges$to,]) == 0)
  print(paste("Delete", length(noContextEdges), "without context"))

  delete_edges(g, noContextEdges)
}

attributedLinkClustering <- function(g, features, simFun, ...) {

  m <- Matrix(0, nrow=ecount(g), ncol=ecount(g))
  rownames(m) <- 1:nrow(m)
  colnames(m) <- 1:ncol(m)
  linkSim <- getLinkSimilarity(g, features, simFun)

  m[cbind(linkSim$e1, linkSim$e2)] <- linkSim$similarity
  m[cbind(linkSim$e2, linkSim$e1)] <- linkSim$similarity
  rs <- rowSums(m) > 0
  print(paste("remove", length(which(!rs)), "edges", sep=" "))
  m <- m[rs,rs]
  (1 - as.dist(m)) %>%
     #hclust(method="ward.D2")
    hclust(method="single")
}


