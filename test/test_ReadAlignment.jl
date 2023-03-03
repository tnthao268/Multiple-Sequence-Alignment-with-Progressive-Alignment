include("../BioAlignment_Versuch.jl")
using Test

#--------------
#Test check_DNA

s1 = "ATGCC"
@test check_DNA(s1)
s2 = "ATTGD"
@test !check_DNA(s2)

#------------------
#Test ReadAlignment

dict_pair = Dict("A" => seq, "B" => ref) #seq, ref are sequences in BioAlignment_Versuch
aln = readDNAAlignment(dict_pair)

#Compare aligned sequences in ReadAlignment and in BioAlignment_Versuch
aligned_seq = collect(values(aln.dict_alnpair))
@test aln_seq in aligned_seq
@test aln_ref in aligned_seq