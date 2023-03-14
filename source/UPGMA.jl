import Pkg
Pkg.add("MCPhyloTree")
using MCPhyloTree

include("../test/test_DistanceMatrix.jl")
records
dict_records = createDictionary(records).dict_records

export upgma, guildTreeInstruction

#Cluster::GeneralNode
cluster = upgma(dm,test_leaf_names)

#create an instruction, what is read from guild tree
function guildTreeInstruction(cluster::GeneralNode)
    guildTree = []
    function addInstruction(cluster::GeneralNode)
        if length(cluster.children) == 0
            return cluster.name
        end
        
        for x in cluster.children
            pair = [addInstruction(y) for y in cluster.children]
            println("pair: $pair")
            push!(guildTree, pair)
            return pair
        end
    end
    addInstruction(cluster)
    guildTree
end
 
#Using
guildTree = guildTreeInstruction(cluster)
#=
4-element Vector{Any}:
 ["A", "D"]
 Any[["A", "D"], "C"]
 ["B", "E"]
 Vector[Any[["A", "D"], "C"], ["B", "E"]]
=#

