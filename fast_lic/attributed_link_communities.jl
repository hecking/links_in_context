using LightGraphs
using DataFrames
using DataFrames: groupby
using GraphDataFrameBridge
using CSV
using DelimitedFiles
using NamedArrays
using SparseArrays
using Clustering
using Lazy

function getLinkSimilarityNetFeat(anchorNode::Int64, endpoint1::Int64, endpoint2::Int64,
  g, features::Array{Int8})::Float64

  commonRelevantAttributes = findall(
    i -> features[endpoint1, i] == 1 && features[endpoint2, i] == 1,
    1:size(features, 2)
  )

  if (length(commonRelevantAttributes) > 0)

    relevantAttributes = reduce((x, i) ->
      ((features[endpoint1, i] == 1) || features[endpoint2, i] == 1) ? x + 1 : x,
      1:size(features, 2),
      init = 0
    )

    @fastmath sum(features[common_neighbors(g, endpoint1, endpoint2),
      commonRelevantAttributes]) /
      (length(relevantAttributes) *
        length(union(neighbors(g, endpoint1), neighbors(g, endpoint2))))

  else
    0
  end
end

function getLinkSimilarityNetwork(anchorNode::Int64, endpoint1::Int64,
  endpoint2::Int64, g, features=nothing)::Float64

   @fastmath length(common_neighbors(g, endpoint1, endpoint2)) /
    length(union(neighbors(g, endpoint1), neighbors(g, endpoint2)))
end

function getLinkDistance(g, features, simFun)

  simMatrix = spzeros(ne(g), ne(g))
  #simMatrix = [1.0 for i in 1:ne(g), j in 1:ne(g)]

  edgeDf = @> DataFrame(end1 = src.(LightGraphs.edges(g)),
   end2 = dst.(LightGraphs.edges(g)), id = 1:ne(g))
print("run 1")
  Threads.@threads for subframe in [i for i in groupby(edgeDf, :end1)]
    for i in 1:(nrow(subframe) - 1), j in (i + 1):nrow(subframe)
      dist = -simFun(subframe[1, :end1],
        subframe[i, :end2], subframe[j, :end2], g, features)

        #@inbounds simMatrix[subframe.id[i], subframe.id[j]] = dist

        @inbounds simMatrix[subframe.id[j], subframe.id[i]] = dist
    end
  end
print("run 2")
  Threads.@threads for subframe in [i for i in groupby(edgeDf, :end2)]

    for i in 1:(nrow(subframe) - 1), j in (i + 1):nrow(subframe)

        #@inbounds simMatrix[subframe.id[i], subframe.id[j]] = dist

        @inbounds simMatrix[subframe.id[j], subframe.id[i]] = -simFun(subframe[1, :end2],
          subframe[i, :end1], subframe[j, :end1], g, features)
    end
  end

  return simMatrix
end


function partitionDensity(g, clusters, atts...)

  f = function(cl)

    sg, vmap = induced_subgraph(g,
    view(collect(edges(g)), cl))

    @fastmath ne(sg) * ((ne(sg) - (nv(sg) - 1))) /
      ((nv(sg) - 2) * (nv(sg) - 1))
  end

  sum(
    [f(cl) for cl in clusters if length(cl) > 1]
  ) * 2 / ne(g)
end

function attributedPartitionDensity(g, clusters, features)

  numFeatures = sum(view(features, 1, :)) # number of features

  m = sum(map(e -> sum(view(atts, e.src,:) .* view(atts, e.dst, :)), edges(g)))

  f = function(cl)
    sg, vmap = induced_subgraph(g,
    view(collect(edges(g)), cl))
    mc = sum(map(e -> sum(view(atts, e.src,:) .* view(atts, e.dst, :)), edges(sg)))

    @fastmath mc * mc /
      ((numFeatures * nv(sg)) * (nv(sg) - 1))
  end

  sum(
    [f(cl) for cl in clusters if length(cl) > 1]
  ) * 2 / m
end

function getCut(dendogram, cutpoint)

  mapping = cutree(dendogram, h=cutpoint)

  map(x -> findall(y -> y == x, mapping), unique(mapping))
end

function cutLinkDendogram(dendogram, g, features)

  cuts = -0.8:0.05:0
  linkClusters = []
  by(DataFrame(cut = cuts), :cut, :cut => function(cp)
    linkClusters = getCut(dendogram, cp[1])
    (num_clusters = length(linkClusters),
     clusters = linkClusters,
     pd = partitionDensity(g, linkClusters),
     apd = attributedPartitionDensity(g, linkClusters, features)
    )
  end)
end

function linkClustersToNodeClusters(g, linkClusters)

  map(lc -> induced_subgraph(g, lc)[2], linkClusters)
end

function attributedLinkClustering(g, features, simFun, atts...)

  @> getLinkDistance(g, features, simFun) hclust cutLinkDendogram(g, atts)
end
