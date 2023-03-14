import Pkg
Pkg.add("MCPhyloTree")
using MCPhyloTree


dm = [  0.0    4.0  -9.0    4.0   11.0
        4.0    0.0  14.0  -11.0  -20.0
        -9.0   14.0   0.0   10.0   28.0
        4.0  -11.0  10.0    0.0   57.0
        11.0  -20.0  28.0   57.0    0.0]
leaf_names = ["A","B","C","D","E"]

export upgma, guildTreeInstruction

#Cluster::GeneralNode
cluster = upgma(dm,leaf_names)

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
 ["A", "C"]
 Any[["A", "C"], "D"]
 ["B", "E"]
 Vector[Any[["A", "C"], "D"], ["B", "E"]]
=#

