include("../../source/ProgressiveAligment/MSA.jl")
using Test

#Test getPairScore
@test EDNAFULL['A','A'] == getPairScore("-A","AA")
@test gap_open == getPairScore("AA","T-")
@test gap_open == getPairScore("-T","A-")
@test gap_open == getPairScore("-T","--")
@test gap_open == getPairScore("--","--")
@test getPairScore("A-","--") == 0
@test getPairScore("T-","A-") == 0
@test gap_extend == getPairScore("--","TT")

#Test getSoP_at_a_pos
seqs1 = ["ATAT";
         "ATT";
         "-TTA"]
seqs2 = ["CAA";
         "--AA"]
sop = getSoP_at_a_pos(seqs1,seqs2,1,1)
#(['A','A','-'],['C','-']) --> (AC + A-)*2 + -C + -- 
test_sop = (-4 + gap_open)*2 + gap_open + 0
@test sop == test_sop
sop2 = getSoP_at_a_pos(seqs1,seqs2,2,2)
test_sop2 = -4*3 + -5 + -1*2
@test sop2 == test_sop2

#Test setTraceback
@test setTraceback(5,4,3) == "dia"
@test setTraceback(4,5,3) == "ver"

#Test msa_globalAlignment
aln_seqs1 = Record[Record("seq1","AA-GC"), Record("seq2","AATGC")]
aln_seqs2 = Record[Record("seq3","ACTC")]
msa_globalAlignment(aln_seqs1,aln_seqs2)
#= expected result (short form)
seq1 AA-GC
seq2 AATGC
seq3 ACT-C
=#

#result after alignment
msa.aln_seqs1
#=
2-element Vector{Record}:
 Record("seq1", "AA-GC")
 Record("seq2", "AATGC")
 =#
 msa.aln_seqs2
 #=
 1-element Vector{Record}:
 Record("seq3", "ACT-C")
 =#