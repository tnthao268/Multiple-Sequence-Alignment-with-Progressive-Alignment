include("../DistanceMatrix.jl")
using Test

s1 = FastaRecord("Seq1","TCAGGATGAAC")
s2 = FastaRecord("Seq2","ATCACGATGAACC")
s3 = FastaRecord("Seq3","ATCAGGAATGAATCC")
s4 = FastaRecord("Seq4","TCACGATTGAATCGC")
s5 = FastaRecord("Seq5","TCAGGAATGAATCGC")
records = [s1,s2,s3,s4,s5]

distanceMatrix = [  0 34 36 32 41;
                    34 0 44 36 27;
                    36 44 0 40 58;
                    32 36 40 0 57;
                    41 27 58 57 0]

dm = createDistanceMatrix(records)

@test dm == convert(Matrix{Float64},distanceMatrix)

