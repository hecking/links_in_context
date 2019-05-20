using LightGraphs
using DataFrames
using GraphDataFrameBridge
using CSV
using DelimitedFiles
using NamedArrays
using Lazy
using MetaGraphs

function readNetwork(path)

    @> path CSV.read(header=["end1", "end2"], delim=" ", types=[String, String]) MetaGraph(:end1, :end2)
end

function readAttributes(attFile, egoFile, nodes = nothing)
    atts = @> attFile readdlm(' ', Int16) NamedArray

    egoAtts = @> egoFile readdlm vec

    setnames!(atts, Array(string.(Int.(atts[:,1]))), 1)

    atts = atts[:,2:end]

    if (!isnothing(nodes))
        atts = atts[nodes,:]
    end

    atts[:,egoAtts .== 1]

    return atts
end

function readAttributedNetwork(circleId)
    g = readNetwork(string(circleId, ".edges"))
    atts = readAttributes(string(circleId, ".feat"),
        string(circleId, ".egofeat"),
        [g[n,:name] for n in vertices(g)])

    return (g, atts)
end
