#export check_DNA; readDNAAlignment
using Pkg
Pkg.add("BioAlignments")
using BioAlignments
include("DataReader.jl")

#--------------
#ReadAlignment
struct ReadAlignment
    score::Int64
    aln_pair::Vector{Record}
end

#---------------------------
# check if sequence is DNA
function check_DNA(seq::String)
    nu = "ATGC"
    for i in seq
        if i ∉ nu; return false; end
    end
    return true
end

#--------------------------------------------------------
#=  aligns two sequences and save the aligned sequences
    returns a ReadAlignment, which contains the score of the alignment and 2 aligned sequences
=#
function readPairwiseAlignment(pair::Vector{Record})
    
    #check, if there are 2 sequences for pairwise alignment and if they are DNA seuqences
    @assert length(pair) == 2 "there are more or less than two sequences to do a pairwise alignment"
    foreach(pair) do s
        @assert check_DNA(s.sequence) "$(s.description) is not a (aligned) DNA sequence"
    end

    #pairwise alignment
    scoremodel = AffineGapScoreModel(EDNAFULL,gap_open = -5, gap_extend = -1)
    res = pairalign(GlobalAlignment(),pair[1].sequence,pair[2].sequence,scoremodel)
    aln = alignment(res)

    #sequences in String and create list of AlignmentRecord(description, aligned sequence)
    aln_seq1 = reduce(*,[x for (x,_) in aln])
    aln_seq2 = reduce(*,[y for (_,y) in aln])
    aln_records = [Record(pair[1].description, aln_seq1),
                   Record(pair[2].description, aln_seq2)]
    
    println("score: $(score(res))\n$aln")
    ReadAlignment(score(res),aln_records)
end