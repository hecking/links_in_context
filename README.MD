# Links in Context

R implementation of the Links in Context approach as an extension of the Link Communities method (Ahn et al., 2010) for attributed networks.

Hecking, T., & Hoppe, H. U. (2018). Links in Context: Detecting and Describing the Nested Structure of Communities in Node-Attributed Networks. In International Conference on Complex Networks and their Applications (pp. 195-206). Springer, Cham.

Ahn, Y. Y., Bagrow, J. P., & Lehmann, S. (2010). Link communities reveal multiscale complexity in networks. nature, 466(7307), 761.

**Example usage:**
Lazega Lawyers dataset [[https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm]]
```R
source("attributed_link_communities.R")
source("read_sienna.R")

atts <- readAttributesSienna("ELattr.dat") 
gFriend <- readGraphSienna("ELfriend.dat") %>% removeNoContextEdges(atts)
```

Run link clustering
```R
# Possible link similarity measures:
## getLinkSimilarityNetFeat: Standart LiC method (network and feature based) 
## getLinkSimilarityFeature: Feature overlap only
## getLinkSimilarityFeature: Network only (original Link Communities method)

linkDendogram <- attributedLinkClustering(gFriend, atts, getLinkSimilarityNetFeat)
plot(linkDendogram)
```

Find best cut of dendogram and transform link clusters to node clusters
```R
linkClusters <- cutLinkDendogram(linkDendFriends, gFriend, atts, 7) # 7 = number of attribute types
nodeClusters <- linkClustersToNodeClusters(gFriend, lcFriends$link_clusters[[which.max(lcFriends$apd)]])
```
