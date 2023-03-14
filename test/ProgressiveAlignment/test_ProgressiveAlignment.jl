include("../../source/ProgressiveAlignment/ProgressiveAlignment.jl")
using Test

#Test createGroup
@test createGroup([["A","B"],"C"]) == "ABC" 
@test createGroup("A") == "A"

#Test simpleGuildTree
include("../test_UPGMA2.jl")
test_guildTree
simple = simpleGuildTree(test_guildTree)
test_simple = [["A","D"],["AD","C"],["B","E"],["ADC","BE"]]
@test simple == test_simple

#Run progressiveAlignment
include("../test_DistanceMatrix.jl")
@time result = progressiveAlignment(simple,dict_records)
#= 0.675598 seconds (553.43 k allocations: 28.117 MiB, 54.01% compilation time)
Dict{String, Record} with 5 entries:
  "B" => Record("Seq2", "ATCACGA-TGAA-C-C")
  "A" => Record("Seq1", "-TCAGGAT-GAA---C")
  "C" => Record("Seq3", "ATCAGGAATGAATC-C")
  "D" => Record("Seq4", "-TCACGATTGAATCGC")
  "E" => Record("Seq5", "-TCAGGAATGAATCGC")
=#
