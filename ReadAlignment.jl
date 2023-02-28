using BioAlignments
include("BioAlignment_Versuch.jl")

struct ReadAlignment
    score::Int64
    aln_seq1::Vector{Char}
    aln_seq2::Vector{Char}
end

# check if sequence is DNA
function check_DNA(seq::String)
    nu = "ATGC"
    for i in seq
        if i ∉ nu; return false; end
    end
    return true
end

#=  aligns two sequences and save the aligned sequences
    returns a ReadAlignment(score,aln_seq1,aln_seq2)
=#
function readDNAAlignment(seq1::String, seq2::String)
    scoremodel = AffineGapScoreModel(EDNAFULL,gap_open = -10, gap_extend = -1)
    
    @assert check_DNA(seq1) "seq1 isn't DNA"
    @assert check_DNA(seq2) "seq2 isn't DNA"

    res = pairalign(GlobalAlignment(),seq1,seq2,scoremodel)
    aln = alignment(res)

    aln_seq1 = [x for (x,_) in aln]
    aln_seq2 = [y for (_,y) in aln]

    ReadAlignment(score(res),aln_seq1,aln_seq2)
end
