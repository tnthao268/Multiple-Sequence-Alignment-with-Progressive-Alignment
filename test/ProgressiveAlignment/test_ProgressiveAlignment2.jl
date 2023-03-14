include("../../source/ProgressiveAlignment/ProgressiveAlignment2.jl")
using Test

#Test createGroup
@test createGroup([["A","B"],"C"]) == "ABC" 
@test createGroup("A") == "A"

#Test nestedInstruction
include("../test_UPGMA2.jl")
nestedinstruction = nestedInstruction(test_guildTree)
test_nestedinstruction = [[["A", "D"], "C"],["B", "E"]]
@test nestedinstruction == test_nestedinstruction

#Run progressiveAlignment2
include("../test_DistanceMatrix.jl")
@time result2 = progressiveAlignment2(nestedinstruction, dict_records)

#=  0.270992 seconds (59.14 k allocations: 2.740 MiB, 24.42% compilation time: 100% of which was recompilation)
Dict{String, Record} with 5 entries:
"B" => Record("Seq2", "ATCACGA-TGAA-C-C")
"A" => Record("Seq1", "-TCAGGAT-GAA---C")
"C" => Record("Seq3", "ATCAGGAATGAATC-C")
"D" => Record("Seq4", "-TCACGATTGAATCGC")
"E" => Record("Seq5", "-TCAGGAATGAATCGC")
=#

#Compare 2 methods progressiveAlignment
include("test_ProgressiveAlignment.jl")
for leaf in keys(result)
    @test result[leaf].sequence == result2[leaf].sequence 
end