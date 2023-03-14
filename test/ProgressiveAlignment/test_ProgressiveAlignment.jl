include("../../source/ProgressiveAlignment/ProgressiveAlignment.jl")
using Test

#Test createGroup
@test createGroup([["A","B"],"C"]) == "ABC" 
@test createGroup("A") == "A"

#Test simpleGuildTree
include("../test_UPGMA.jl")
test_guildTree
simple = simpleGuildTree(test_guildTree)
test_simple = [["A","C"],["AC","D"],["B","E"],["ACD","BE"]]
@test simple == test_simple

#Run progressiveAlignment
include("../test_DistanceMatrix.jl")
@time result = progressiveAlignment(simple,test_dict_records)
#= 0.588170 seconds (606.81 k allocations: 30.859 MiB, 67.52% compilation time)
Dict{String, Record} with 5 entries:
  "B" => Record("Seq2", "-ATCACGATGAACC--")
  "A" => Record("Seq1", "-TCAGGA-TGAA---C")
  "C" => Record("Seq3", "ATCAGGAATGAATC-C")
  "D" => Record("Seq4", "-TCACGATTGAATCGC")
  "E" => Record("Seq5", "-TCAGGAATGAATCGC")
=#
