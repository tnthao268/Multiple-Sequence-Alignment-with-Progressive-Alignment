module UPGMA
import Pkg
Pkg.add("MCPhyloTree")

m = [0 17 21 31 23; 17 0 30 34 21; 21 30 0 28 39; 31 34 28 0 43; 23 21 39 43 0]
m = convert(Array{Float64}, m)
leaf_names = ["A" , "B", "C", "D" , "E"]

using MCPhyloTree


# change name of a cluster into combination of its children' name (e.g "AB")

function children(cluster::GeneralNode)
    i = ""
    for j in cluster.children
        if length(i) > 1 && length(j.name) >= 1 #  if the name is the combination of more than 1 leaf
            i = "(" * i * ")" * "," * "(" * j.name * ")"  # add also the brackets in between the clusters, eg: "AB(E)"
       
        else
            i = i * j.name # otherwise no brackets: eg: "AB"
        end
    end
    return i
end

# test: children(find_by_name(tree,"Cluster_1"))

# return tree with changed cluster's name which are combination of its children' name (e.g "AB")

function change_cluster_name(tree::GeneralNode)
    post_list = post_order(tree)
    for node in post_list
        if length(node.children) > 0
            node.name = children(node)
        end
    end
    return tree
end


# change_cluster_name(tree) 

# return raw list of strings which are cluster's name

function cluster_list(tree::GeneralNode)
    cluster_list = []
    post_list = post_order(tree) # traverse the tree and return a list of nodes
    for node in post_list
        if length(node.children) > 1 # only adds nodes which have more than 0 children into the cluster list
            push!(cluster_list, node.name)
        end
    end
    return cluster_list
end
#= e.g 4-element Vector{Any}:
 "CD"
 "AB"
 "AB(E)"
 "CD(AB(E))"
 =#

# still try, NOT CORRECT YET
# method to split string into meaningful name of sequences to be aligned and store them in a list 
function split_name_sequences(cluster_list::Vector{Any})
    list = []
    for name in cluster_list
        if length(name) < 3
            push!(list, split(name,""))
        else
            push!(list, split(name,"),(", limit = 2))
        
        end
        
    end

    
    new_lst = []
    for i in list
        l = []
        for j in i 
            push!(l,replace(j, ")" => "", "(" => ""))
        end
        push!(new_lst,l)
    
    end

    # filter comma in a string
    new_lst2 = []

    for i in new_lst
        l = []
        for j in i 
            push!(l,replace(j, "," => ""))
        end
        push!(new_lst2,l)
    
    end
    

    return unique(new_lst2) # to avoid duplicate elements in a list 
end

#= eg:4-element Vector{Any}:
 Any["C", "D"]
 Any["A", "B"]
 Any["AB", "E"]
 Any["CD", "ABE"]
 =#



split_name_sequences(cluster_list(change_cluster_name(upgma(m, leaf_names))))

end