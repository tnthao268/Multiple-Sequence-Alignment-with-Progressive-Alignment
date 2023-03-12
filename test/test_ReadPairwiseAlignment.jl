include("../BioAlignment_Versuch.jl")
include("../ReadPairwiseAlignment.jl")
using Test

#--------------
#Test check_DNA

s1 = "ATGCC"
@test check_DNA(s1)
s2 = "ATTGD"
@test !check_DNA(s2)

#------------------
#Test ReadAlignment : test with the result from BioAlignment_Versuch

pair = [Record("A",seq), Record("B",ref)] #Variables seq and ref are sequences in BioAlignment_Versuch
aln = readDNAAlignment(pair)

#Compare aligned sequences in ReadAlignment and in BioAlignment_Versuch
aligned_seqs = aln.aln_pair
@test aln_seq == aligned_seqs[1].sequence
@test aln_ref == aligned_seqs[2].sequence
@test s == aln.score