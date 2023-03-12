include(dirname(@__FILE__()) * "/../source/UPGMA.jl")

m = [0 17 21 31 23; 17 0 30 34 21; 21 30 0 28 39; 31 34 28 0 43; 23 21 39 43 0]
m = convert(Array{Float64}, m)
leaf_names = ["A" , "B", "C", "D" , "E"]

tree = upgma(m, leaf_names)

list = split_name_sequences(cluster_list(change_cluster_name(upgma(m, leaf_names))))
@test list == [["C", "D"],["A", "B"], ["AB", "E"],["CD", "ABE"]]

