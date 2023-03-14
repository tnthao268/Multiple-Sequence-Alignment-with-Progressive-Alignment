include("../test/test_DistanceMatrix.jl")
include("../source/UPGMA2.jl")
using Test

records
dm = createDistanceMatrix(records)
dict = createDictionary(records)
dict_records = dict.dict_records
leaf_names = dict.leaf_names

cluster = upgma(dm,leaf_names)
guildTree = guildTreeInstruction(cluster)
test_guildTree = [["A", "D"], [["A", "D"], "C"], ["B", "E"],[[["A", "D"], "C"], ["B", "E"]]]
@test guildTree == test_guildTree