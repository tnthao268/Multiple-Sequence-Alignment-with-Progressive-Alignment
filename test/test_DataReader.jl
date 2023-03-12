include(dirname(@__FILE__()) * "/../source/DataReader.jl")


using Test


fasta_record  = Record("test string1", "ACGT")
fasta_record2 = Record("test string1", "ACGT")

println("testing FastaRecord...")
@test fasta_record.description == "test string1"
@test fasta_record.sequence == "ACGT"
@test fasta_record2.description == "test string1"
@test fasta_record2.sequence == "ACGT"


# @test fasta_record == fasta_record2

sequences = readSequences(dirname(@__FILE__()) * "/../data/sequence.fasta")
# print(sequences)
@test length(sequences) > 0
@test length(sequences) == 1
@test sequences[1].description == "NM_001103547.4 Drosophila melanogaster Cyclic-AMP response element binding protein B (CrebB), transcript variant I, mRNA"


sequences1 = readSequences(dirname(@__FILE__()) * "/../data/input_test_sequences.faa")
print(sequences1)
@test length(sequences1) > 0
@test length(sequences1) == 10
@test sequences1[1].description == "1"
@test sequences1[1].sequence == "LWYVRMHMYR"
@test sequences1[10].description == "10"
@test sequences1[10].sequence == "QNHHVGNGPA"






