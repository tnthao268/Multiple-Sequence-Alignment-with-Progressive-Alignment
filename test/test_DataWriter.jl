include(dirname(@__FILE__()) * "/../source/DataWriter.jl")
include(dirname(@__FILE__()) * "/../source/DataReader.jl")


fasta_record  = Record("test string1", "ACGT")
fasta_record2 = Record("test string2", "ACGT")


str = split_sequence(fasta_record.sequence,2)
@test str == ["AC" , "GT"]

record  = Record("NM_001103547.4 Drosophila melanogaster Cyclic-AMP response", "ACGTATGCTCGTAGCTGACGTAGCTGACGTAGCTGACGTAGCTGACGACG")
record2 = Record("NM_001103547.5 Drosophila melanogaster Cyclic-AMP response", "ACGTATGCTCGTAGGTGACGTAGCTGACGTAGCTGACGTAGCTGACGACG")

dict_seq = Dict("A" => record, "B" => record2)

writeSequences("try_seq1.txt",dict_seq)

