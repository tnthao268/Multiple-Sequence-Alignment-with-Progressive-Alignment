include("../source/DistanceMatrix.jl")
using Test

s1 = Record("Seq1","TCAGGATGAAC")
s2 = Record("Seq2","ATCACGATGAACC")
s3 = Record("Seq3","ATCAGGAATGAATCC")
s4 = Record("Seq4","TCACGATTGAATCGC")
s5 = Record("Seq5","TCAGGAATGAATCGC")
records = [s1,s2,s3,s4,s5]

#test createDistanceMatrix
distanceMatrix = [  0.0    4.0  -9.0    4.0   11.0
                    4.0    0.0  14.0  -11.0  -20.0
                    -9.0   14.0   0.0   10.0   28.0
                    4.0  -11.0  10.0    0.0   57.0
                    11.0  -20.0  28.0   57.0    0.0]

dm = createDistanceMatrix(records)

@test dm == distanceMatrix

#test createDictionary
test_leaf_names = ["A","B","C","D","E"]
test_dict_records = Dict(zip(test_leaf_names,records))
dict = createDictionary(records)
@test dict.dict_records == test_dict_records
@test dict.leaf_names == test_leaf_names