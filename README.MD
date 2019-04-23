# Links in Context

R implementation of the Links in Context approaoch.

Hecking, T., & Hoppe, H. U. (2018, December). Links in Context: Detecting and Describing the Nested Structure of Communities in Node-Attributed Networks. In International Conference on Complex Networks and their Applications (pp. 195-206). Springer, Cham.

**Example usage:**
Lazega Lawyers dataset [[
```R
source("attributed_link_communities.R")
source("read_sienna.R")

atts <- readAttributesSienna("ELattr.dat") 
gFriend <- readGraphSienna("ELfriend.dat") %>% removeNoContextEdges(atts)
```

Run link clustering
```R
linkDendogram <- attributedLinkClustering(gFriend, atts, getLinkSimilarityNetFeat)
plot(linkDendogram)
```

Find best cut of dendogram and transform link clusters to node clusters
```R
linkClusters <- cutLinkDendogram(linkDendFriends, gFriend, atts, 7) # 7 = number of attribute types
nodeClusters <- linkClustersToNodeClusters(gFriend, lcFriends$link_clusters[[which.max(lcFriends$apd)]])
```
