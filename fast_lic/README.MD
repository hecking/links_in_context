# Fast LiC

This is a faster implementation of different versions of link clustering in Julia (https://julialang.org/).
It runs almost twice as fast as the R implementation but post-processing is much less user friendly.

### Example usage:

```Julia
include("fb_network_reader.jl")
include("atttributed_link_communties.jl")

g, atts = readAttributedNetwork(107)

result = attributedLinkClustering(g, atts, getLinkSimilarityNetFeat)
```
