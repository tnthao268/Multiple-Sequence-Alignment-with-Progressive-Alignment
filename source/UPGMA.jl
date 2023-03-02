import Pkg
Pkg.add("MCPhyloTree")

m = [0 17 21 31 23; 17 0 30 34 21; 21 30 0 28 39; 31 34 28 0 43; 23 21 39 43 0]
m = convert(Array{Float64}, m)
leaf_names = ["A" , "B", "C", "D" , "E"]

using MCPhyloTree
tree = upgma(m, leaf_names)
node_a = find_by_name(tree,"A")

# function to return all names of tree's leaves
function return_leaves(tree::GeneralNode)
    leafnames = []
        for x in get_leaves(tree)
            push!(leafnames,x.name)
        end
    return leafnames
end

# function to return an array which contains leaves of the tree based on name, return GeneralNode
function find_leaves(leafname::Vector)
    leaves = []
    for l in leafname
        push!(leaves, find_by_name(tree,l))
    end
    return leaves
end

# find_leaves(return_leaves(tree))

# find cluster to do pairwise alignment

function findCluster(tree::GeneralNode)
    # create new 2 dimension array which contains pair of alignment
    saved_cluster = []
    # for every leaf of a tree (every sequence in distance matrix)
    for m in find_leaves(return_leaves(tree))
        for n in find_leaves(return_leaves(tree)) 
            # if 2 sequences are grouped together (with same length)
            if path_length(m.mother,m) == path_length(n.mother,n) && m != n
                push!(saved_cluster,[m,n]) # add these two in the array
            end
        end
    end
    return saved_cluster
end

#c = findCluster(tree)
# c[1][1].name
# c[1][2].name

# function to get rid of duplicate pairs of sequences 
function cluster_name(c::Vector{Any})
    name = []
    # index of pair in c
    for i in 1:length(c)
        # add the name of sequences in pair to name array
        push!(name,[c[i][1].name,c[i][2].name])
    end
    # iterate through name array
    for i in name
        # sort elements in the alphabetical order
        sort!(i) # sort array in-place (sort() does not work)
    end
    return collect(Set(name)) # duplicate elements are ruled out
end
        
# cluster_name(findCluster(tree))

# function to return next node for alignment with paired alignment
function next_node_align(aligned_node::GeneralNode)
    for i in find_leaves(return_leaves(tree))
        if aligned_node.mother.mother == i.mother
                return i.name
        end
    end
end


# Add new alignment pairs in set
function automatic_add_alignment_pair(list::Vector)
    # find node by name (first element of every pair in set)
    for j in list
        node = find_by_name(tree,j[1])
        if next_node_align(node) !== nothing
            push!(list,[next_node_align(node), node.name])

        end

    end
    return list
end

automatic_add_alignment_pair(cluster_name(findCluster(tree)))