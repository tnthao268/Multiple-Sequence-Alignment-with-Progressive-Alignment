include("../source/DistanceMatrix.jl")
using Test

s1 = Record("Seq1","TCAGGATGAAC")
s2 = Record("Seq2","ATCACGATGAACC")
s3 = Record("Seq3","ATCAGGAATGAATCC")
s4 = Record("Seq4","TCACGATTGAATCGC")
s5 = Record("Seq5","TCAGGAATGAATCGC")
records = [s1,s2,s3,s4,s5]

#test createDistanceMatrix
distanceMatrix = [  0 34 36 32 41;
                    34 0 44 36 27;
                    36 44 0 40 58;
                    32 36 40 0 57;
                    41 27 58 57 0]

dm = createDistanceMatrix(records)

@test dm == convert(Matrix{Float64},distanceMatrix)

#test createDictionary
test_leaf_names = ["A","B","C","D","E"]
test_dict_records = Dict(zip(test_leaf_names,records))
dict = createDictionary(records)
@test dict.dict_records == test_dict_records
@test dict.leaf_names == test_leaf_names